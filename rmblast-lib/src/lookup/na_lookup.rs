//! Nucleotide lookup table and subject scanning.
//!
//! Mirrors NCBI's adaptive table-width strategy:
//!
//! - `lut_width` is chosen from the word size and approximate query length
//!   following NCBI's `BlastChooseNaLookupTable` thresholds.  For a 1 MB
//!   query with word_size=14 this gives lut_width=12 instead of 14, keeping
//!   the backbone at 4^12×4 = 67 MB rather than 4^14×24 = 6.4 GB.
//!
//! - The table uses a compact linked-list layout matching NCBI's MBLookupTable:
//!     backbone[hash] = 1 + first_query_position  (0 → empty)
//!     next_pos[i]    = previous backbone value for the same hash
//!   so traversal is: q_off = backbone[hash]-1; while q_off valid: next = next_pos[q_off+1]-1
//!
//! - Subject is scanned at stride `scan_step = word_size − lut_width + 1`.
//!   Any word_size-mer match generates at least one seed (on the correct diagonal).
//!
//! Zero-initializing `backbone` with `vec![0u32; N]` on Linux mmap-maps the
//! virtual address space but the OS supplies zero pages lazily — physical memory
//! is only committed for the ~query_len entries actually written during build.

use crate::encoding::NUCL_SENTINEL;

/// A word hit: positions in query and subject where a lut-word matched.
#[derive(Debug, Clone, Copy)]
pub struct WordHit {
    pub query_offset: u32,
    pub subject_offset: u32,
}

/// Compact nucleotide lookup table (MBLookupTable-style).
pub struct NaLookupTable {
    pub word_size: usize,
    /// Key width stored in the backbone: min ≤ word_size, chosen by query size.
    pub lut_width: usize,
    /// Subject scan stride = word_size − lut_width + 1.
    pub scan_step: usize,
    /// backbone[hash] = 1-indexed first query position with that hash (0 = empty).
    backbone: Vec<u32>,
    /// Linked list continuation: next_pos[1-indexed q_off] = previous backbone value.
    next_pos: Vec<u32>,
}

/// Choose lut_width following NCBI's BlastChooseNaLookupTable thresholds.
/// `approx_entries` ≈ total indexable query positions (query_len).
fn choose_lut_width(word_size: usize, approx_entries: usize) -> usize {
    match word_size {
        4..=6 => word_size,
        7 => if approx_entries < 250 { 6 } else { 7 },
        8 => if approx_entries < 8_500 { 7 } else { 8 },
        9 => {
            if approx_entries < 1_250 { 7 }
            else if approx_entries < 21_000 { 8 }
            else { 9 }
        }
        10 => {
            if approx_entries < 1_250 { 7 }
            else if approx_entries < 8_500 { 8 }
            else if approx_entries < 18_000 { 9 }
            else { 10 }
        }
        11 => {
            if approx_entries < 12_000 { 8 }
            else if approx_entries < 180_000 { 10 }
            else { 11 }
        }
        12 => {
            if approx_entries < 8_500 { 8 }
            else if approx_entries < 18_000 { 9 }
            else if approx_entries < 60_000 { 10 }
            else if approx_entries < 900_000 { 11 }
            else { 12 }
        }
        _ => {
            // word_size ≥ 13: NCBI default case
            if approx_entries < 8_500 { 8 }
            else if approx_entries < 300_000 { 11 }
            else { 12 }
        }
    }
}

impl NaLookupTable {
    /// Build a lookup table from a single BLASTNA-encoded query (forward strand only).
    ///
    /// `query` must include leading/trailing sentinels at index 0 and len-1.
    /// DUST masking should be applied before calling: masked bases (≥ 4) are skipped
    /// and their k-mers are not indexed, matching NCBI's lookup_segments behavior.
    pub fn build(query: &[u8], word_size: usize) -> Self {
        let n = query.len().saturating_sub(2);
        Self::build_from_slice(query, n, word_size)
    }

    /// Build a lookup table from a single BLASTNA-encoded query with an explicit
    /// `approx_entries` override for lut_width selection.
    ///
    /// Used by Combined mode: NCBI's BlastChooseNaLookupTable uses `2 × query_len`
    /// when both forward and RC contexts are counted, giving a different (often
    /// larger) lut_width than the single-strand case.  This constructor lets the
    /// caller pass `2 × fwd_n` to reproduce that selection while still indexing
    /// only forward-query k-mers.
    pub fn build_with_approx_entries(query: &[u8], approx_entries: usize, word_size: usize) -> Self {
        Self::build_from_slice(query, approx_entries, word_size)
    }

    /// Build a combined LUT from DUST-masked forward and RC query buffers.
    ///
    /// Mirrors NCBI's `BlastMBLookupTableNew` with both strand contexts in one table.
    /// `fwd_scan` and `rc_scan` are BLASTNA-encoded (sentinel-wrapped) DUST-masked
    /// buffers; pass the original query's `encode_iupac` output with `dust_mask`
    /// applied.
    ///
    /// The combined buffer layout is:
    ///   [sentinel][fwd_bases][sentinel][rc_bases][sentinel]
    ///
    /// `approx_entries = 2 × fwd_len` matches NCBI's `EstimateNumTableEntries` which
    /// sums positions from both context-0 and context-1 `lookup_segments`.
    ///
    /// ## Known efficiency flaw (shared with NCBI)
    /// Context-1 (RC-query) k-mers generate ~half the raw hits.  For non-palindromic
    /// sequences every context-1 hit fails the exact-match extension step (wrong query
    /// positions vs subject), so they are filtered before reaching the diagonal
    /// tracker.  The wasted work is a faithful reproduction of NCBI behavior.
    ///
    /// ## Known edge-case bug (shared with NCBI)
    /// When `word_size == lut_width` no exact extension is performed, so context-1
    /// hits would NOT be filtered and would produce alignments at wrong coordinates.
    pub fn build_combined(fwd_scan: &[u8], rc_scan: &[u8], word_size: usize) -> Self {
        let fwd_n = fwd_scan.len().saturating_sub(2);
        let rc_n  = rc_scan.len().saturating_sub(2);
        // Combined buffer: [sentinel][fwd_bases][sentinel][rc_bases][sentinel]
        // The inner sentinel resets valid_run during building, keeping the two
        // contexts independent while sharing a single backbone and next_pos array.
        let mut combined = Vec::with_capacity(3 + fwd_n + rc_n);
        combined.push(NUCL_SENTINEL);
        combined.extend_from_slice(&fwd_scan[1..1 + fwd_n]);
        combined.push(NUCL_SENTINEL);
        if rc_n > 0 {
            combined.extend_from_slice(&rc_scan[1..1 + rc_n]);
        }
        combined.push(NUCL_SENTINEL);
        // Use 2×fwd_n as approx_entries (both strand contexts counted, matching NCBI).
        Self::build_from_slice(&combined, 2 * fwd_n, word_size)
    }

    /// Core table builder.  `approx_entries` overrides the default n-derived estimate
    /// so callers can pass 2×fwd_len for the combined case.
    fn build_from_slice(query: &[u8], approx_entries: usize, word_size: usize) -> Self {
        assert!(word_size >= 4 && word_size <= 14, "word_size must be 4-14");

        let n = query.len().saturating_sub(2); // real base count (may span two contexts)
        let lut_width = choose_lut_width(word_size, approx_entries);
        let scan_step = word_size - lut_width + 1;
        let backbone_size = 1usize << (2 * lut_width);

        // Zero-initialized backbone: on Linux the allocator uses mmap for large
        // sizes, so the OS supplies zero pages lazily.  Physical memory is only
        // committed for the ~n entries actually written.
        let mut backbone = vec![0u32; backbone_size];
        // next_pos is indexed by 1-based position into the real-bases region.
        let mut next_pos = vec![0u32; n + 1];

        if n < lut_width {
            return NaLookupTable { word_size, lut_width, scan_step, backbone, next_pos };
        }

        let bases = &query[1..n + 1];
        let mask = backbone_size - 1;
        let mut word: usize = 0;
        let mut valid_run = 0usize;

        for (i, &b) in bases.iter().enumerate() {
            if b >= 4 {
                valid_run = 0;
                word = 0;
                continue;
            }
            word = ((word << 2) | (b as usize)) & mask;
            valid_run += 1;
            if valid_run >= lut_width {
                // q_off is 0-based start of the lut_width-mer within the real-bases region
                let q_off = i + 1 - lut_width;
                let slot = q_off + 1; // 1-based into next_pos
                next_pos[slot] = backbone[word];
                backbone[word] = slot as u32;
            }
        }
        NaLookupTable { word_size, lut_width, scan_step, backbone, next_pos }
    }

    /// Scan `subject` for all query lut-word matches, appending to `hits`.
    /// Subject includes leading/trailing sentinels.
    ///
    /// `scan_start` is the first position (in real-bases space) to scan.  For
    /// plus-strand scans this is always 0.  For minus-strand (RC-subject) scans
    /// it must be `(subject_real_len - lut_width) % scan_step` so that the
    /// covered subject positions map to the same plus-strand positions that NCBI
    /// covers via context-1 hits in its combined-LUT plus-strand scan.
    pub fn scan_subject(&self, subject: &[u8], hits: &mut Vec<WordHit>) {
        self.scan_subject_from(subject, 0, hits);
    }

    /// Like `scan_subject` but starts at `scan_start` (real-bases offset).
    pub fn scan_subject_from(&self, subject: &[u8], scan_start: usize, hits: &mut Vec<WordHit>) {
        let lw = self.lut_width;
        let step = self.scan_step;
        let n = subject.len().saturating_sub(2);
        if n < lw {
            return;
        }

        let bases = &subject[1..n + 1];
        let mask = self.backbone.len() - 1;

        if step == 1 {
            // Streaming hash — no ambiguous-base restart in this path.
            // scan_start is ignored for step=1 (already covers every position).
            let mut word: usize = 0;
            let mut valid_run: usize = 0;
            for (i, &b) in bases.iter().enumerate() {
                if b >= 4 {
                    valid_run = 0;
                    word = 0;
                    continue;
                }
                word = ((word << 2) | (b as usize)) & mask;
                valid_run += 1;
                if valid_run >= lw {
                    let s_off = (i + 1 - lw) as u32;
                    self.emit_hits(word, s_off, hits);
                }
            }
        } else {
            // Strided scan: compute hash fresh at each scan position.
            let mut sp = scan_start;
            while sp + lw <= n {
                let mut hash: usize = 0;
                let mut valid = true;
                for k in 0..lw {
                    let b = bases[sp + k];
                    if b >= 4 {
                        valid = false;
                        break;
                    }
                    hash = (hash << 2) | (b as usize);
                }
                if valid {
                    self.emit_hits(hash, sp as u32, hits);
                }
                sp += step;
            }
        }
    }

    #[inline(always)]
    fn emit_hits(&self, hash: usize, s_off: u32, hits: &mut Vec<WordHit>) {
        let mut slot = self.backbone[hash] as usize;
        while slot != 0 {
            hits.push(WordHit {
                query_offset: (slot - 1) as u32,
                subject_offset: s_off,
            });
            slot = self.next_pos[slot] as usize;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encoding::encode_iupac;

    #[test]
    fn test_exact_match() {
        let query = encode_iupac(b"ACGTACGT");  // 8 bases
        let subject = encode_iupac(b"TTACGTACGTTT");
        let tbl = NaLookupTable::build(&query, 8);
        let mut hits = Vec::new();
        tbl.scan_subject(&subject, &mut hits);
        assert!(!hits.is_empty());
        let found = hits.iter().any(|h| h.query_offset == 0 && h.subject_offset == 2);
        assert!(found, "expected hit at q=0, s=2; got {:?}", hits);
    }

    #[test]
    fn test_no_match() {
        let query = encode_iupac(b"AAAAAAAA");
        let subject = encode_iupac(b"CCCCCCCCCCCC");
        let tbl = NaLookupTable::build(&query, 8);
        let mut hits = Vec::new();
        tbl.scan_subject(&subject, &mut hits);
        assert!(hits.is_empty());
    }

    #[test]
    fn test_lut_width_selection_large_query() {
        // 1 MB query with word_size=14 should select lut_width=12, scan_step=3
        let lw = choose_lut_width(14, 1_000_000);
        assert_eq!(lw, 12);
        assert_eq!(14 - lw + 1, 3); // scan_step
    }

    #[test]
    fn test_lut_width_selection_small_query() {
        // Small query (Alu, 311 bp) with word_size=14 → lut_width=8, scan_step=7
        let lw = choose_lut_width(14, 311);
        assert_eq!(lw, 8);
        assert_eq!(14 - lw + 1, 7); // scan_step
    }

    #[test]
    fn test_stride_finds_match() {
        // word_size=14 with small query (311 bases) → lut_width=8, scan_step=7
        // A 14-base match must be found despite stride-7 scanning.
        let query_seq = b"ACGTACGTACGTAC"; // 14 bases
        let query = encode_iupac(query_seq);
        // Subject: 10 Ns then the query then 10 Ns
        let mut subj_seq: Vec<u8> = b"NNNNNNNNNN".to_vec();
        subj_seq.extend_from_slice(query_seq);
        subj_seq.extend_from_slice(b"NNNNNNNNNN");
        let subject = encode_iupac(&subj_seq);

        // For test to work with small query (14 bases) → lut_width=8, but our
        // query is tiny (14 bases, approx_entries=14), so lut_width=8, scan_step=7.
        let tbl = NaLookupTable::build(&query, 14);
        let mut hits = Vec::new();
        tbl.scan_subject(&subject, &mut hits);
        assert!(!hits.is_empty(), "stride scan should find the 14-mer match");
        // All hits should be on the correct diagonal (s_off - q_off = 10)
        let correct_diag = hits.iter().any(|h| {
            h.subject_offset as i64 - h.query_offset as i64 == 10
        });
        assert!(correct_diag, "at least one hit on diagonal 10; hits={:?}", hits);
    }
}
