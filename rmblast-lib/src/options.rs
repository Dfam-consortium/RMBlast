//! Search parameters mirroring the rmblastn CLI options in scope.

/// Gap alignment algorithm to use (always ALIGN_EX for rmblastn with matrix scoring).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GapAlignAlgo {
    AlignEx,
}

/// Threading mode (mirrors NCBI -mt_mode).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MtMode {
    SplitByDb,
    SplitByQueries,
}

/// Lookup-table seeding strategy.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeedMode {
    /// NCBI-faithful combined LUT: both forward and RC query k-mers in one table.
    ///
    /// Matches NCBI's BlastMBLookupTableNew behavior exactly, including:
    /// - lut_width chosen with approx_entries = 2 × query_len (both strands counted).
    /// - DUST applied to the query before LUT construction (masked positions excluded).
    ///
    /// Known efficiency flaw (shared with NCBI): ~half the raw hits are from context-1
    /// (RC-query) k-mers and always fail the exact-match extension step for
    /// non-palindromic sequences, wasting ~2× seeding work.
    ///
    /// Known edge-case bug (shared with NCBI): if word_size == lut_width no exact
    /// extension is needed, so context-1 hits are NOT filtered and produce alignments
    /// at wrong query coordinates for palindromic k-mers.
    Combined,

    /// Efficient single-strand LUT: forward query k-mers only.
    ///
    /// lut_width chosen with approx_entries = query_len (one strand).
    /// DUST applied to the query before LUT construction.
    /// Minus-strand seeding still uses the forward-query LUT against the revcomp
    /// subject (matching NCBI's minus-strand scan semantics without the wasted
    /// context-1 hits).
    SeparateStrands,
}

/// All tunable search parameters used by the alignment engine.
#[derive(Debug, Clone)]
pub struct SearchParams {
    // Scoring
    pub gap_open: i32,
    pub gap_extend: i32,
    pub matrix_name: String,

    // Seeding
    pub word_size: usize,

    // X-dropoff (raw score units)
    pub xdrop_ungap: i32,
    pub xdrop_gap: i32,
    pub xdrop_gap_final: i32,

    // Filtering
    /// Score cutoff for the **preliminary** gapped stage — NOT a floor on the
    /// reported score.  A prelim HSP scoring at or above this is promoted to
    /// traceback, but traceback re-aligns under `xdrop_gap_final` and may return a
    /// lower score, which is then reported as-is.  This matches NCBI rmblastn:
    /// `Blast_TracebackFromHSPList` does not re-test the cutoff (see
    /// PORTING_NOTES.md §9.4).  Measured against rmblastn 2.17.1: ~0.2% of hits fall
    /// below a cutoff of 200, ~2% below a cutoff of 95 (minimum observed score: 1).
    ///
    /// Two later stages *do* re-apply it: `complexity_adjust` (drops adjusted scores
    /// below the cutoff) and the cut-HSP `reevaluate_gapped` path.  Neither closes the
    /// gap for the main path.  Callers needing a hard floor must filter
    /// `AlignResult::hsp.score` themselves.
    pub min_raw_gapped_score: i32,
    /// Ungapped pre-filter keep threshold.  `None` = use the historical fixed
    /// fallback `min_raw_gapped_score / 2`.  `Some(c)` = a Karlin-Altschul–derived
    /// cutoff (see `search::ka_cutoff`), computed per-search from the matrix+gap
    /// ALP params and the average DB sequence length, as NCBI's
    /// `BlastInitialWordParametersNew` does.
    pub ungapped_cutoff: Option<i32>,
    pub complexity_adjust: bool,
    pub dust: bool,
    /// Maximum query coverage by a higher-scoring HSP before this HSP is suppressed.
    /// 80 = RepeatMasker's usual setting (it passes `-mask_level` explicitly);
    /// 101 = effectively disabled, and the DEFAULT, matching NCBI rmblastn's
    /// `-mask_level` default of -1.  Callers that do not set `-mask_level` must
    /// get no masklevel filtering, or self-vs-self searches collapse to their
    /// self-hits alone (a 100%-coverage self-hit dominates every other HSP).
    pub mask_level: u32,

    // Multithreading
    pub num_threads: usize,
    pub mt_mode: MtMode,

    // Seeding strategy
    pub seed_mode: SeedMode,

}

impl Default for SearchParams {
    fn default() -> Self {
        Self {
            gap_open: 4,
            gap_extend: 4,
            matrix_name: String::new(),
            word_size: 8,
            xdrop_ungap: 20,
            xdrop_gap: 30,
            xdrop_gap_final: 100,
            min_raw_gapped_score: 0,
            ungapped_cutoff: None,
            complexity_adjust: false,
            dust: true,
            mask_level: 101,
            num_threads: 1,
            mt_mode: MtMode::SplitByDb,
            seed_mode: SeedMode::Combined,
        }
    }
}
