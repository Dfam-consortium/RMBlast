//! Search parameters mirroring the rmblastn CLI options in scope.

/// Threading mode (mirrors NCBI -mt_mode).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MtMode {
    SplitByDb,
    SplitByQueries,
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

}

impl SearchParams {
    /// The gapped HSP cutoff, NCBI's `hit_params->cutoffs[ctx].cutoff_score`.
    /// With a positive --min-raw-gapped-score it is that value (rmblastn's
    /// matrix-only scoring, BlastHitSavingParametersUpdate).  With 0 NCBI
    /// derives it from E = 10, starting from 1: for a matrix without
    /// Karlin-Altschul statistics that is exactly 1, and for one with them it
    /// is the same E = 10 threshold the --ncbi-compat E-value cull applies.
    /// So it is never below 1: 3.0.12 used 0 and kept score-0 remainders
    /// (identity.matrix printed hits like q 3940-3939 with score 0).
    pub fn gapped_cutoff(&self) -> i32 {
        self.min_raw_gapped_score.max(1)
    }

    /// NCBI's `cutoff_score_max`, which caps the ungapped cutoff: half of a
    /// positive --min-raw-gapped-score (matrix-only scoring), else the
    /// E-derived gapped cutoff, which is at least 1 -- exactly 1 for a matrix
    /// without Karlin-Altschul statistics (the RMBlast fallback in
    /// BlastInitialWordParametersUpdate then uses it as the ungapped cutoff).
    pub fn cutoff_score_max(&self) -> i32 {
        if self.min_raw_gapped_score > 0 { self.min_raw_gapped_score / 2 } else { 1 }
    }
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
        }
    }
}
