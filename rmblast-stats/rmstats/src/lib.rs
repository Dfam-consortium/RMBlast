//! # rmstats
//!
//! Faithful Rust port of the Karlin-Altschul E-value / bit-score statistics
//! of NCBI BLAST 2.17.0 with the RMBlast (-RMH-) patches, for use by the
//! rmblastn Rust port.
//!
//! Sources (all in ncbi-blast-2.17.0+-src/c++/src/algo/blast/core/):
//! * blast_stat.c — KA tables (blastn_values_*, rmblast_*_values),
//!   Blast_KarlinBlkGappedLoadFromTables, Blast_KarlinBlkNuclGappedCalc,
//!   Blast_GetNuclAlphaBeta, BLAST_ComputeLengthAdjustment,
//!   BLAST_KarlinStoE_simple, BlastKarlinEtoS_simple, BLAST_Cutoffs, and the
//!   ungapped machinery (Blast_KarlinLambdaNR, BlastKarlinLtoH,
//!   BlastKarlinLHtoK, score/residue frequencies).
//! * blast_setup.c — Blast_ScoreBlkKbpGappedCalc three-mode hierarchy for
//!   custom matrices, BLAST_CalcEffLengths.
//! * blast_hits.c — Blast_HSPListGetEvalues / Blast_HSPListGetBitScores
//!   sentinel semantics.
//!
//! The high-level entry point is [`RmStats`]: build it once per
//! (matrix, gap costs, query length, database totals) and call
//! [`RmStats::evalue`] / [`RmStats::bit_score`] per HSP raw score.

pub mod eff_lengths;
pub mod evalue;
pub mod karlin;
pub mod length_adjust;
pub mod ncbi_math;
pub mod nucl_tables;
pub mod rmblast_tables;
pub mod tables_data;

#[cfg(feature = "alp-fit")]
pub mod alp;

pub use eff_lengths::{calc_eff_lengths_custom_matrix, EffLengths, SearchSpaceInputs};
pub use evalue::{hsp_bit_score, hsp_evalue, karlin_stoe_simple};
pub use karlin::{kbp_ungapped_calc_blastna, KarlinBlk, BLASTNA_SIZE};
pub use rmblast_tables::{kbp_gapped_calc_custom_matrix, MatrixCliOverrides};

/// Bundled statistics context for one rmblastn search: the gapped Karlin
/// block resolved through the three-mode hierarchy plus the effective search
/// space for the query/database geometry.
///
/// Mirrors the NCBI flow: Blast_ScoreBlkKbpGappedCalc →
/// BLAST_CalcEffLengths → per-HSP Blast_HSPListGetEvalues /
/// Blast_HSPListGetBitScores.
#[derive(Debug, Clone, Copy)]
pub struct RmStats {
    /// Gapped Karlin block (may be the Mode-3 sentinel).
    pub kbp_gap: KarlinBlk,
    /// Effective search space and length adjustment for this query context.
    pub eff: EffLengths,
    /// sbp->round_down: always false on the custom-matrix path.
    pub round_down: bool,
}

impl RmStats {
    /// Set up statistics for a custom-matrix (read_in_matrix) blastn search,
    /// the rmblastn case.
    ///
    /// * `matrix_name` — basename as passed to -matrix (e.g.
    ///   "14p35g.matrix"); matched case-insensitively against the baked
    ///   table.
    /// * `gap_open`, `gap_extend` — the -gapopen/-gapextend costs.
    /// * `cli` — the -matrix_lambda/-matrix_k/-matrix_alpha/-matrix_beta
    ///   overrides (all zero when not supplied).
    /// * `query_length` — length of the query sequence (each blastn strand
    ///   context has this same length and search space).
    /// * `db_length` — total database length (sum of subject lengths), after
    ///   any -dblen override.
    /// * `db_num_seqs` — number of database sequences, after any override.
    /// * `eff_searchsp_override` — the -searchsp override; 0 when unset.
    /// * `kbp_std` — the ungapped Karlin block for this context, needed only
    ///   for the Mode-2 alpha/beta fallback (compute with
    ///   [`kbp_ungapped_calc_blastna`]); pass `None` to mirror an invalid
    ///   ungapped context.
    #[allow(clippy::too_many_arguments)]
    pub fn new_custom_matrix(
        matrix_name: &str,
        gap_open: i32,
        gap_extend: i32,
        cli: &MatrixCliOverrides,
        query_length: i32,
        db_length: i64,
        db_num_seqs: i32,
        eff_searchsp_override: i64,
        kbp_std: Option<&KarlinBlk>,
    ) -> RmStats {
        let kbp_gap = kbp_gapped_calc_custom_matrix(matrix_name, gap_open, gap_extend, cli);
        let inputs = SearchSpaceInputs {
            query_length,
            db_length,
            db_num_seqs,
            eff_searchsp_override,
        };
        let eff = calc_eff_lengths_custom_matrix(
            &inputs,
            &kbp_gap,
            kbp_std,
            cli,
            gap_open,
            gap_extend,
            true, // gapped_calculation
        );
        RmStats { kbp_gap, eff, round_down: false }
    }

    /// True if usable statistics are available (not the Mode-3 sentinel).
    pub fn has_stats(&self) -> bool {
        self.kbp_gap.is_valid()
    }

    /// E-value for a raw HSP score (1.0 when statistics are unavailable,
    /// matching Blast_HSPListGetEvalues).
    pub fn evalue(&self, raw_score: i32) -> f64 {
        hsp_evalue(
            raw_score,
            Some(&self.kbp_gap),
            self.eff.eff_searchsp,
            true,
            self.round_down,
        )
    }

    /// Bit score for a raw HSP score (0.0 when statistics are unavailable,
    /// matching Blast_HSPListGetBitScores).
    pub fn bit_score(&self, raw_score: i32) -> f64 {
        hsp_bit_score(raw_score, Some(&self.kbp_gap))
    }
}
