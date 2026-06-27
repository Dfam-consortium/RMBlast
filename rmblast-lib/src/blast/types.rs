// blast_types.rs — faithful port of NCBI BLAST core data structures.
//
// Primary sources (NCBI BLAST 2.17.0+):
//   blast_def.h        → BlastOffsetPair, BLAST_SequenceBlk (SeqBlk here)
//   blast_hits.h       → BlastSeg, BlastHSP, BlastHSPList
//   blast_extend.h     → BlastUngappedData, BlastInitHSP, BlastInitHitList
//   blast_gapalign.h   → BlastGapDP, BlastGapAlignStruct
//   gapinfo.h          → GapEditScript, GapPrelimEditBlock, GapStateArrayStruct,
//                        EGapAlignOpType

// ── blast_def.h ──────────────────────────────────────────────────────────────

/// COMPRESSION_RATIO: 4 nucleotide bases packed into 1 byte.
pub const COMPRESSION_RATIO: i32 = 4;

/// Analogous to BLAST_SequenceBlk.
///
/// NCBI maintains two forms of a sequence in the same struct:
///   sequence       – primary form (BLASTNA 1-byte/base for query; packed
///                    NCBI2NA 4-bases/byte for database subjects in prelim phase)
///   sequence_start – one position before `sequence`; holds the sentinel byte
///
/// For the query, a third form is also maintained:
///   compressed_nuc_seq       – 4-to-1 packed NCBI2NA built from the BLASTNA
///                              query; used by the ungapped extension scanner
///   compressed_nuc_seq_start – allocation root (3 bytes before compressed_nuc_seq)
///
/// In this Rust port we store the BLASTNA form and the packed form as separate
/// Vecs so that lifetime / ownership is explicit.  The `sequence()` helper
/// returns a slice starting at base-index 0 (i.e., past the sentinel byte),
/// matching C's `seq_blk->sequence`.
pub struct SeqBlk {
    /// BLASTNA (1-byte-per-base) storage.
    /// Layout: data[0] = sentinel (0), data[1..=length] = bases.
    /// `sequence` in NCBI C code points to data[1].
    pub blastna: Vec<u8>,

    /// Packed NCBI2NA storage for this sequence (built by
    /// `blast_compress_blastna_sequence`).
    ///
    /// Layout mirrors NCBI's compressed_nuc_seq:
    ///   packed[0..2] = pre-sequence bytes (right-justified partial packs)
    ///   packed[3..]  = packed bases; packed[3+k] covers bases k..k+3 with
    ///                  base k in bits 7:6 (MSB) and base k+3 in bits 1:0.
    ///
    /// `compressed_nuc_seq` in C points to packed[3].
    pub packed: Vec<u8>,

    /// Sequence length in bases (not including sentinel).
    pub length: i32,
}

impl SeqBlk {
    /// Slice of BLASTNA bases, analogous to `seq_blk->sequence` in NCBI.
    /// Index 0 corresponds to the first base.
    #[inline]
    pub fn sequence(&self) -> &[u8] {
        &self.blastna[1..] // skip sentinel
    }

    /// Mutable BLASTNA bases.
    #[inline]
    pub fn sequence_mut(&mut self) -> &mut [u8] {
        &mut self.blastna[1..]
    }

    /// `compressed_nuc_seq` pointer equivalent: slice starting at the
    /// packed base-0 byte (index 3 in `self.packed`).
    #[inline]
    pub fn compressed_nuc_seq(&self) -> &[u8] {
        &self.packed[3..]
    }
}

/// BlastOffsetPair (blast_def.h).
/// Holds the query/subject offsets of an initial word match.
#[derive(Clone, Copy, Debug, Default)]
pub struct BlastOffsetPair {
    pub q_off: u32,
    pub s_off: u32,
}

// ── blast_extend.h ───────────────────────────────────────────────────────────

/// BlastUngappedData (blast_extend.h).
#[derive(Clone, Copy, Debug, Default)]
pub struct BlastUngappedData {
    pub q_start: i32,
    pub s_start: i32,
    pub length:  i32,
    pub score:   i32,
}

/// BlastInitHSP (blast_extend.h).
#[derive(Clone, Copy, Debug, Default)]
pub struct BlastInitHsp {
    pub offsets:       BlastOffsetPair,
    pub ungapped_data: BlastUngappedData,
    /// Whether ungapped_data is valid (C uses a nullable pointer).
    pub has_ungapped:  bool,
}

/// BlastInitHitList (blast_extend.h) — all initial HSPs for one subject.
///
/// C uses a growable array with explicit `total`/`allocated` fields; Rust Vec
/// provides the same semantics automatically.
#[derive(Clone, Debug, Default)]
pub struct BlastInitHitList {
    pub init_hsp_array: Vec<BlastInitHsp>,
}

impl BlastInitHitList {
    pub fn new() -> Self { Self::default() }
    pub fn reset(&mut self) { self.init_hsp_array.clear(); }
    pub fn total(&self) -> i32 { self.init_hsp_array.len() as i32 }
}

/// DiagStruct (blast_extend.h) — per-diagonal last-hit bookkeeping.
#[derive(Clone, Copy, Debug, Default)]
pub struct DiagStruct {
    pub last_hit: i32,  // offset of the last hit (+ diag_table.offset)
    pub flag:     bool, // TRUE if a hit was saved on this diagonal
}

/// BLAST_DiagTable (blast_extend.h) — diagonal table for two-hit tracking.
#[derive(Debug)]
pub struct BlastDiagTable {
    /// Per-diagonal last_hit / flag entries (length = diag_array_length).
    pub hit_level_array:  Vec<DiagStruct>,
    /// Length of the most recent hit on each diagonal (two-hit mode only).
    pub hit_len_array:    Vec<u8>,
    /// Smallest power of 2 ≥ query_length + subject_length.
    pub diag_array_length: i32,
    /// diag_array_length − 1 (mask for modular indexing).
    pub diag_mask:        i32,
    /// Running offset added to subject positions so the array doesn't need
    /// to be zeroed between subjects.
    pub offset:           i32,
    /// Window size for the two-hit requirement (0 → one-hit mode).
    pub window:           i32,
    pub multiple_hits:    bool,
    pub actual_window:    i32,
}

impl BlastDiagTable {
    /// Allocate a fresh diagonal table for the given query length.
    /// Mirrors BlastExtendWordNew (blast_extend.c line 110).
    pub fn new(query_length: i32, multiple_hits: bool, window_size: i32) -> Self {
        // Round up to the next power of 2.
        let mut len = 1i32;
        while len < query_length { len <<= 1; }
        len <<= 1; // one extra power-of-2 margin (matches NCBI)
        let n = len as usize;
        Self {
            hit_level_array:   vec![DiagStruct::default(); n],
            hit_len_array:     if multiple_hits { vec![0u8; n] } else { Vec::new() },
            diag_array_length: len,
            diag_mask:         len - 1,
            offset:            0,
            window:            window_size,
            multiple_hits,
            actual_window:     0,
        }
    }

    /// Update the offset after scanning a subject of the given length.
    /// Mirrors Blast_ExtendWordExit (blast_extend.c line 162).
    pub fn exit_subject(&mut self, subject_length: i32) {
        self.offset += subject_length + 1;
        if self.offset > i32::MAX / 2 {
            self.hit_level_array.iter_mut().for_each(|d| d.last_hit = 0);
            self.hit_len_array.iter_mut().for_each(|x| *x = 0);
            self.offset = 0;
        }
    }
}

// ── gapinfo.h ────────────────────────────────────────────────────────────────

/// EGapAlignOpType (gapinfo.h).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum GapAlignOp {
    Del     = 0, // gap in query (deletion from query)
    Sub     = 3, // substitution / match
    Ins     = 6, // gap in subject (insertion into query)
}

/// GapEditScript (gapinfo.h) — run-length encoded alignment path.
#[derive(Clone, Debug, Default)]
pub struct GapEditScript {
    pub ops: Vec<(GapAlignOp, i32)>, // (op_type, count) pairs
}

/// GapPrelimEditScript (gapinfo.h) — single op entry before compression.
#[derive(Clone, Copy, Debug)]
pub struct GapPrelimEditScript {
    pub op_type: GapAlignOp,
    pub num:     i32,
}

/// GapPrelimEditBlock (gapinfo.h) — dynamically grown array of prelim ops.
#[derive(Clone, Debug, Default)]
pub struct GapPrelimEditBlock {
    pub edit_ops: Vec<GapPrelimEditScript>,
    pub last_op:  Option<GapAlignOp>,
}

impl GapPrelimEditBlock {
    pub fn new() -> Self { Self::default() }

    /// Analogous to GapPrelimEditBlockReset.
    pub fn reset(&mut self) {
        self.edit_ops.clear();
        self.last_op = None;
    }

    /// Analogous to GapPrelimEditBlockAdd.
    pub fn add(&mut self, op_type: GapAlignOp, num: i32) {
        if let Some(last) = self.edit_ops.last_mut() {
            if last.op_type == op_type {
                last.num += num;
                return;
            }
        }
        self.edit_ops.push(GapPrelimEditScript { op_type, num });
        self.last_op = Some(op_type);
    }
}

// ── blast_hits.h ─────────────────────────────────────────────────────────────

/// BlastSeg (blast_hits.h) — one sequence segment within an HSP.
#[derive(Clone, Copy, Debug, Default)]
pub struct BlastSeg {
    pub frame:        i16,
    pub offset:       i32, // start of hsp (0-based inclusive)
    pub end:          i32, // end of hsp (0-based exclusive)
    pub gapped_start: i32, // where the gapped extension started
}

/// BlastHSP (blast_hits.h).
#[derive(Clone, Debug)]
pub struct BlastHsp {
    pub score:    i32,
    pub query:    BlastSeg,
    pub subject:  BlastSeg,
    pub context:  i32,
    pub gap_info: GapEditScript,
}

// ── blast_gapalign.h ─────────────────────────────────────────────────────────

/// BlastGapDP (blast_gapalign.h) — one DP cell.
#[derive(Clone, Copy, Debug, Default)]
pub struct BlastGapDP {
    /// Score of best path ending in a match at this position.
    pub best:     i32,
    /// Score of best path ending in a gap at this position.
    pub best_gap: i32,
}

/// MININT: used to initialise DP cells that are "unreachable".
/// Mirrors `#define MININT INT4_MIN/2` in blast_gapalign.c line 58.
pub const MININT: i32 = i32::MIN / 2;

/// BlastGapAlignStruct (blast_gapalign.h) — workspace for gapped extension.
///
/// Owns the reusable DP array (dp_mem) and traceback buffers (fwd/rev prelim
/// edit blocks).  All other fields are per-call outputs written by the
/// alignment functions.
pub struct BlastGapAlignStruct {
    /// Reusable DP scratch array.  Length grows on demand.
    pub dp_mem:       Vec<BlastGapDP>,
    /// Forward (right extension) preliminary traceback.
    pub fwd_prelim_tback: GapPrelimEditBlock,
    /// Reverse (left extension) preliminary traceback.
    pub rev_prelim_tback: GapPrelimEditBlock,
    /// X-dropoff currently in use.
    pub gap_x_dropoff: i32,
    // Output fields written by alignment functions:
    pub query_start:   i32,
    pub query_stop:    i32,
    pub subject_start: i32,
    pub subject_stop:  i32,
    pub score:         i32,
}

impl BlastGapAlignStruct {
    /// Analogous to BLAST_GapAlignStructNew.
    pub fn new(x_dropoff: i32, dp_alloc: usize) -> Self {
        Self {
            dp_mem: vec![BlastGapDP::default(); dp_alloc.max(100)],
            fwd_prelim_tback: GapPrelimEditBlock::new(),
            rev_prelim_tback: GapPrelimEditBlock::new(),
            gap_x_dropoff: x_dropoff,
            query_start:   0,
            query_stop:    0,
            subject_start: 0,
            subject_stop:  0,
            score:         0,
        }
    }
}

// ── scoring parameters (subset needed by gapped alignment) ───────────────────

/// Analogous to BlastScoringParameters (blast_parameters.h) — just the
/// penalty fields used by the DP routines.
#[derive(Clone, Copy, Debug)]
pub struct ScoringParams {
    pub gap_open:   i32,
    pub gap_extend: i32,
}
