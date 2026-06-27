// blast/nalookup.rs — faithful port of NCBI nucleotide lookup table construction.
//
// Primary sources (NCBI BLAST 2.17.0+):
//   blast_lookup.h   → ComputeTableIndex, PV_ARRAY_BTS, PV_ARRAY_MASK,
//                      PV_ARRAY_TYPE, PV_SET, PV_TEST
//   blast_lookup.c   → BlastLookupAddWordHit, BlastLookupIndexQueryExactMatches
//   blast_nalookup.h → NaLookupBackboneCell, BlastNaLookupTable, NA_HITS_PER_CELL
//   blast_nalookup.c → s_BlastNaLookupFinalize, BlastNaLookupTableNew

// ── Constants ────────────────────────────────────────────────────────────────

/// BITS_PER_NUC (blast_nalookup.c line 43): bits per nucleotide in BLASTNA.
pub const BITS_PER_NUC: i32 = 2;

/// NA_HITS_PER_CELL (blast_nalookup.h line 109): max inline hits per cell.
pub const NA_HITS_PER_CELL: usize = 3;

/// PV_ARRAY_BTS (blast_lookup.h line 43): bits to shift from index → pv word.
pub const PV_ARRAY_BTS: u32 = 5;

/// PV_ARRAY_MASK (blast_lookup.h line 44): bit-position mask within a pv word.
pub const PV_ARRAY_MASK: i32 = 31;

// ── PV helpers ───────────────────────────────────────────────────────────────

/// PV_SET (blast_lookup.h line 49): mark backbone index as occupied.
#[inline]
pub fn pv_set(pv: &mut [u32], index: i32) {
    pv[(index >> PV_ARRAY_BTS) as usize] |= 1u32 << (index & PV_ARRAY_MASK);
}

/// PV_TEST (blast_lookup.h line 55): true if backbone index is occupied.
#[inline]
pub fn pv_test(pv: &[u32], index: i32) -> bool {
    pv[(index >> PV_ARRAY_BTS) as usize] & (1u32 << (index & PV_ARRAY_MASK)) != 0
}

// ── NaLookupBackboneCell (blast_nalookup.h lines 113–127) ────────────────────

/// Union payload: inline entries or overflow cursor.
///
/// C uses a `union { Int4 overflow_cursor; Int4 entries[NA_HITS_PER_CELL]; }`.
/// When `num_used <= NA_HITS_PER_CELL`, entries are stored inline; otherwise
/// `overflow_cursor` is the index into the overflow array.
#[derive(Clone, Debug)]
pub enum NaLookupPayload {
    Entries([i32; NA_HITS_PER_CELL]),
    OverflowCursor(i32),
}

impl Default for NaLookupPayload {
    fn default() -> Self {
        NaLookupPayload::Entries([0; NA_HITS_PER_CELL])
    }
}

/// NaLookupBackboneCell (blast_nalookup.h lines 113–127).
#[derive(Clone, Debug, Default)]
pub struct NaLookupBackboneCell {
    pub num_used: i32,
    pub payload:  NaLookupPayload,
}

// ── BlastNaLookupTable (blast_nalookup.h lines 131–156) ──────────────────────

/// BlastNaLookupTable.
pub struct BlastNaLookupTable {
    pub mask:            i32,
    pub word_length:     i32,
    pub lut_word_length: i32,
    pub scan_step:       i32,
    pub backbone_size:   i32,
    pub longest_chain:   i32,
    pub thick_backbone:  Vec<NaLookupBackboneCell>,
    pub overflow:        Vec<i32>,
    pub overflow_size:   i32,
    /// Presence vector: bit i set ↔ thick_backbone[i] has hits.
    pub pv:              Vec<u32>,
}

// ── ComputeTableIndex (blast_lookup.h lines 96–108) ──────────────────────────

/// ComputeTableIndex: rolling hash over `lut_word_length` BLASTNA bytes.
/// `word` must be exactly `lut_word_length` bytes long.
#[inline]
pub fn compute_table_index(word: &[u8]) -> i32 {
    let mut index: i32 = 0;
    for &b in word {
        index = (index << BITS_PER_NUC) | (b as i32);
    }
    index
}

// ── BlastLookupAddWordHit (blast_lookup.c lines 33–77) ───────────────────────

/// Add one query offset to the thin backbone.
///
/// C thin_backbone: `Int4 **` where each non-null chain is
///   `[capacity, num_hits, hit0, hit1, ...]`.
/// Rust representation: `Vec<Vec<i32>>` where:
///   - empty `Vec`  → null (no hits)
///   - chain[0]     → num_hits
///   - chain[1..]   → hit offsets
fn blast_lookup_add_word_hit(
    thin_backbone:   &mut Vec<Vec<i32>>,
    word:            &[u8],
    query_offset:    i32,
) {
    let index = compute_table_index(word) as usize;
    let chain = &mut thin_backbone[index];
    if chain.is_empty() {
        chain.push(1);            // chain[0] = num_hits
        chain.push(query_offset); // chain[1] = first hit
    } else {
        chain[0] += 1;
        chain.push(query_offset);
    }
}

// ── BlastLookupIndexQueryExactMatches (blast_lookup.c lines 79–132) ──────────

/// Index query positions into the thin backbone.
///
/// `query_sequence` — BLASTNA slice (0-based, **no** sentinel byte).
/// `locations`      — `(left, right)` inclusive base ranges.
///
/// Mirrors `void BlastLookupIndexQueryExactMatches(...)`.
pub fn blast_lookup_index_query_exact_matches(
    thin_backbone:   &mut Vec<Vec<i32>>,
    word_length:     i32,
    lut_word_length: i32,
    query_sequence:  &[u8],
    locations:       &[(i32, i32)],
) {
    // invalid_mask: 0xff << BITS_PER_NUC = 0xFC.
    // Any BLASTNA value > 3 is ambiguous (N, R, …).
    let invalid_mask: u8 = 0xff << BITS_PER_NUC as u8;

    for &(from, to) in locations {
        if word_length > to - from + 1 {
            continue;
        }

        // word_target: first query position at which a full lut_word can end.
        // C: `word_target = seq + lut_word_length` (pointer past end of first word).
        let mut word_target: i32 = lut_word_length;
        let mut i: i32 = 0;

        // Walk from the start of this location to `to` (inclusive).
        while from + i <= to {
            let b = query_sequence[(from + i) as usize];

            if i >= word_target {
                let word_start = (from + i - lut_word_length) as usize;
                blast_lookup_add_word_hit(
                    thin_backbone,
                    &query_sequence[word_start..word_start + lut_word_length as usize],
                    from + i - lut_word_length,
                );
            }

            if b & invalid_mask != 0 {
                // Skip past any word that would include this ambiguous base.
                word_target = i + lut_word_length + 1;
            }

            i += 1;
        }

        // Handle the last word (C reads `*seq` before this check; we skip the
        // out-of-range load and just test `i >= word_target`).
        if i >= word_target {
            let word_start = (from + i - lut_word_length) as usize;
            blast_lookup_add_word_hit(
                thin_backbone,
                &query_sequence[word_start..word_start + lut_word_length as usize],
                from + i - lut_word_length,
            );
        }
    }
}

// ── s_BlastNaLookupFinalize (blast_nalookup.c lines 442–546) ─────────────────

/// Convert thin backbone into the compact thick backbone + overflow + pv.
/// Consumes (clears) each non-empty thin_backbone chain.
fn s_blast_na_lookup_finalize(
    thin_backbone: &mut Vec<Vec<i32>>,
    lookup:        &mut BlastNaLookupTable,
) {
    lookup.thick_backbone =
        vec![NaLookupBackboneCell::default(); lookup.backbone_size as usize];

    let pv_words = (lookup.backbone_size >> PV_ARRAY_BTS as i32) as usize + 1;
    lookup.pv = vec![0u32; pv_words];

    // First pass: count overflow slots needed and find longest chain.
    let mut overflow_cells_needed: i32 = 0;
    let mut longest_chain:         i32 = 0;
    for chain in thin_backbone.iter() {
        if chain.is_empty() {
            continue;
        }
        let num_hits = chain[0];
        if num_hits > NA_HITS_PER_CELL as i32 {
            overflow_cells_needed += num_hits;
        }
        longest_chain = longest_chain.max(num_hits);
    }
    lookup.longest_chain = longest_chain;

    if overflow_cells_needed > 0 {
        lookup.overflow = vec![0i32; overflow_cells_needed as usize];
    }

    let mut overflow_cursor: i32 = 0;

    // Second pass: fill thick backbone and overflow.
    for i in 0..lookup.backbone_size as usize {
        if thin_backbone[i].is_empty() {
            continue;
        }

        let num_hits = thin_backbone[i][0];
        lookup.thick_backbone[i].num_used = num_hits;
        pv_set(&mut lookup.pv, i as i32);

        if num_hits <= NA_HITS_PER_CELL as i32 {
            let mut entries = [0i32; NA_HITS_PER_CELL];
            for j in 0..num_hits as usize {
                entries[j] = thin_backbone[i][j + 1];
            }
            lookup.thick_backbone[i].payload = NaLookupPayload::Entries(entries);
        } else {
            lookup.thick_backbone[i].payload =
                NaLookupPayload::OverflowCursor(overflow_cursor);
            for j in 0..num_hits as usize {
                lookup.overflow[overflow_cursor as usize + j] = thin_backbone[i][j + 1];
            }
            overflow_cursor += num_hits;
        }

        thin_backbone[i].clear();
    }

    lookup.overflow_size = overflow_cursor;
}

// ── BlastChooseNaLookupTable (blast_nalookup.c) ──────────────────────────────

/// Choose lut_word_length following NCBI's BlastChooseNaLookupTable thresholds.
///
/// `approx_entries` ≈ total indexable query positions (sum across all contexts:
/// `query_len` for single-strand, `2 × query_len` for combined two-strand).
pub fn choose_lut_width(word_size: usize, approx_entries: usize) -> usize {
    match word_size {
        4..=6 => word_size,
        7 => if approx_entries < 250 { 6 } else { 7 },
        8 => if approx_entries < 8_500 { 7 } else { 8 },
        9 => {
            if      approx_entries < 1_250  { 7 }
            else if approx_entries < 21_000 { 8 }
            else                            { 9 }
        }
        10 => {
            if      approx_entries < 1_250  { 7 }
            else if approx_entries < 8_500  { 8 }
            else if approx_entries < 18_000 { 9 }
            else                            { 10 }
        }
        11 => {
            if      approx_entries < 12_000  { 8 }
            else if approx_entries < 180_000 { 10 }
            else                             { 11 }
        }
        12 => {
            if      approx_entries < 8_500   { 8 }
            else if approx_entries < 18_000  { 9 }
            else if approx_entries < 60_000  { 10 }
            else if approx_entries < 900_000 { 11 }
            else                             { 12 }
        }
        _ => {
            // word_size ≥ 13
            if      approx_entries < 8_500   { 8 }
            else if approx_entries < 300_000 { 11 }
            else                             { 12 }
        }
    }
}

// ── BlastNaLookupTableNew (blast_nalookup.c lines 548–584) ───────────────────

/// Build a BlastNaLookupTable from the query.
///
/// `query_sequence` — BLASTNA slice (0-based, no sentinel).
/// `locations`      — `(left, right)` inclusive base ranges to index.
/// `word_length`    — full word length (e.g. 11).
/// `lut_width`      — lookup table word length (e.g. 8).
pub fn blast_na_lookup_table_new(
    query_sequence: &[u8],
    locations:      &[(i32, i32)],
    word_length:    i32,
    lut_width:      i32,
) -> BlastNaLookupTable {
    let backbone_size = 1i32 << (BITS_PER_NUC * lut_width);

    let mut lookup = BlastNaLookupTable {
        mask:            backbone_size - 1,
        word_length,
        lut_word_length: lut_width,
        scan_step:       word_length - lut_width + 1,
        backbone_size,
        longest_chain:   0,
        thick_backbone:  Vec::new(),
        overflow:        Vec::new(),
        overflow_size:   0,
        pv:              Vec::new(),
    };

    let mut thin_backbone: Vec<Vec<i32>> = vec![Vec::new(); backbone_size as usize];

    blast_lookup_index_query_exact_matches(
        &mut thin_backbone,
        word_length,
        lut_width,
        query_sequence,
        locations,
    );

    s_blast_na_lookup_finalize(&mut thin_backbone, &mut lookup);

    lookup
}
