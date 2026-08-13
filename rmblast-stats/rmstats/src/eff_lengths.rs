//! Effective search space computation, ported from BLAST_CalcEffLengths
//! (blast_setup.c:733-906), restricted to the blastn path.
//!
//! The RMBlast `read_in_matrix` alpha/beta hierarchy (blast_setup.c:829-851)
//! and the Lambda <= 0 sentinel skip (blast_setup.c:873-874) are included.

use crate::karlin::KarlinBlk;
use crate::length_adjust::compute_length_adjustment;
use crate::nucl_tables::get_nucl_alpha_beta;
use crate::rmblast_tables::MatrixCliOverrides;

/// BLAST_REWARD / BLAST_PENALTY (blast_options.h): the fake reward/penalty
/// installed for matrix-only scoring, also used by the Mode-2 alpha/beta
/// fallback.
pub const BLAST_REWARD: i32 = 1;
pub const BLAST_PENALTY: i32 = -3;

/// Inputs describing the search geometry for one query context.
#[derive(Debug, Clone, Copy)]
pub struct SearchSpaceInputs {
    /// Length of this query context (for blastn: the query length; both
    /// strand contexts have the same length and get the same search space).
    pub query_length: i32,
    /// Total database length (user -dblen override if > 0, else the real
    /// total; the caller resolves that choice, mirroring
    /// eff_len_options->db_length vs eff_len_params->real_db_length).
    pub db_length: i64,
    /// Number of database sequences (user -num_seqs override if > 0, else
    /// real count).
    pub db_num_seqs: i32,
    /// User-specified effective search space (-searchsp); 0 when unset.
    /// (s_GetEffectiveSearchSpaceForContext).
    pub eff_searchsp_override: i64,
}

/// Result: per-context effective search space and length adjustment, exactly
/// the two fields BLAST_CalcEffLengths stores on the context.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct EffLengths {
    pub eff_searchsp: i64,
    pub length_adjustment: i32,
}

/// The blastn `read_in_matrix` arm of BLAST_CalcEffLengths for one context.
///
/// * `kbp_gap` — the gapped Karlin block for this context (from
///   `kbp_gapped_calc_custom_matrix`); may be the Mode-3 sentinel.
/// * `kbp_std` — the ungapped Karlin block (kbp_std\[index\]) for the Mode-2
///   alpha/beta fallback; `None` if the ungapped calculation failed (in the C
///   code that would have invalidated the context, so callers should normally
///   pass `Some`).
/// * `cli` — the -matrix_lambda/k/alpha/beta overrides.
/// * `gapped_calculation` — scoring_options->gapped_calculation.
pub fn calc_eff_lengths_custom_matrix(
    inputs: &SearchSpaceInputs,
    kbp_gap: &KarlinBlk,
    kbp_std: Option<&KarlinBlk>,
    cli: &MatrixCliOverrides,
    gap_open: i32,
    gap_extend: i32,
    gapped_calculation: bool,
) -> EffLengths {
    let mut effective_search_space: i64 = inputs.eff_searchsp_override;
    let mut length_adjustment: i32 = 0;

    let query_length = inputs.query_length;
    if query_length > 0 {
        // --- alpha/beta selection (blast_setup.c:829-851) ---
        let (alpha, beta);
        if cli.is_set() {
            // Use user-supplied alpha/beta directly.
            alpha = cli.matrix_alpha;
            beta = cli.matrix_beta; // 0.0 if not supplied
        } else if kbp_gap.h > 0.0 && kbp_gap.lambda > 0.0 {
            // Baked-in table set kbp->H; derive alpha = lambda/H, beta = 0
            // (same approximation BLAST uses for unknown matrices).
            alpha = kbp_gap.lambda / kbp_gap.h;
            beta = 0.0;
        } else {
            // Fall back to fake reward/penalty.
            // NOTE: the C code passes sbp->kbp_std[index] here; if that is
            // NULL the context would have been invalid. With reward/penalty
            // 1/-3 the table lookup itself cannot fail.
            let fallback = KarlinBlk::sentinel();
            let kbp_ungap = kbp_std.unwrap_or(&fallback);
            let (a, b) = get_nucl_alpha_beta(
                BLAST_REWARD,
                BLAST_PENALTY,
                gap_open,
                gap_extend,
                kbp_ungap,
                gapped_calculation,
            )
            .unwrap_or((0.0, 0.0));
            alpha = a;
            beta = b;
        }

        // Mode 3 sentinel: skip the length-adjustment / effective-search-
        // space calculation (blast_setup.c:873-874).
        if !(kbp_gap.lambda <= 0.0) {
            let (la, _converged) = compute_length_adjustment(
                kbp_gap.k,
                kbp_gap.log_k,
                alpha / kbp_gap.lambda,
                beta,
                query_length,
                inputs.db_length,
                inputs.db_num_seqs,
            );
            length_adjustment = la;

            if effective_search_space == 0 {
                let mut effective_db_length: i64 =
                    inputs.db_length - (inputs.db_num_seqs as i64) * (length_adjustment as i64);

                // Just in case effective_db_length < 0
                if effective_db_length <= 0 {
                    effective_db_length = 1;
                }

                effective_search_space =
                    effective_db_length * ((query_length - length_adjustment) as i64);
            }
        }
    }

    EffLengths { eff_searchsp: effective_search_space, length_adjustment }
}
