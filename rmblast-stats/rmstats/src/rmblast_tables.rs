//! RMBlast custom-matrix KA table lookup and the -RMH- three-mode gapped
//! Karlin block hierarchy, ported from the RMBlast-patched blast_stat.c
//! (Blast_KarlinBlkGappedLoadFromTables restricted to the rmblast_* entries)
//! and blast_setup.c (Blast_ScoreBlkKbpGappedCalc, blastn read_in_matrix
//! branch).

use crate::karlin::KarlinBlk;
use crate::tables_data::{RmblastKaEntry, RMBLAST_MATRIX_KA};

/// Status codes of Blast_KarlinBlkGappedLoadFromTables.
#[derive(Debug, PartialEq, Eq)]
pub enum TableLookupError {
    /// status 1: matrix not found
    MatrixNotFound,
    /// status 2: matrix found, but gap open/extend values not supported
    GapCostsNotSupported,
}

/// Find a table entry by matrix name (case-insensitive, as strcasecmp).
pub fn find_matrix_entry(matrix_name: &str) -> Option<&'static RmblastKaEntry> {
    RMBLAST_MATRIX_KA
        .iter()
        .find(|e| e.name.eq_ignore_ascii_case(matrix_name))
}

/// Blast_KarlinBlkGappedLoadFromTables (blast_stat.c), over the RMBlast
/// custom nucleotide matrix entries: exact (case-insensitive) name match,
/// then exact gap-cost match; fills Lambda, K, logK = ln(K), H.
pub fn karlin_blk_gapped_load_from_tables(
    gap_open: i32,
    gap_extend: i32,
    matrix_name: &str,
) -> Result<KarlinBlk, TableLookupError> {
    let entry = find_matrix_entry(matrix_name).ok_or(TableLookupError::MatrixNotFound)?;
    if entry.gap_open == gap_open && entry.gap_extend == gap_extend {
        Ok(KarlinBlk {
            lambda: entry.lambda,
            k: entry.k,
            log_k: entry.k.ln(),
            h: entry.h,
        })
    } else {
        Err(TableLookupError::GapCostsNotSupported)
    }
}

/// The -matrix_lambda/-matrix_k/-matrix_alpha/-matrix_beta command-line
/// overrides (BlastScoringOptions in the RMBlast patch; all default 0.0).
#[derive(Debug, Clone, Copy, Default)]
pub struct MatrixCliOverrides {
    pub matrix_lambda: f64,
    pub matrix_k: f64,
    pub matrix_alpha: f64,
    pub matrix_beta: f64,
}

impl MatrixCliOverrides {
    /// The condition guarding Mode 2 in blast_setup.c: all of lambda, K and
    /// alpha supplied ( > 0).
    pub fn is_set(&self) -> bool {
        self.matrix_lambda > 0.0 && self.matrix_k > 0.0 && self.matrix_alpha > 0.0
    }
}

/// The blastn `read_in_matrix` branch of Blast_ScoreBlkKbpGappedCalc
/// (blast_setup.c:89-126): three-mode hierarchy for the gapped Karlin block
/// of a custom nucleotide matrix.
///
/// Mode 1: baked-in table lookup by matrix name + gap costs.
/// Mode 2: lambda, K and alpha supplied on the command line
///         (H is derived as lambda/alpha).
/// Mode 3: sentinel block (Lambda = K = logK = H = -1); the search proceeds
///         without E-values/bit-scores.
///
/// This never fails — Mode 3 is the terminal fallback, matching the C code's
/// `retval = 0; /* non-fatal */`.
pub fn kbp_gapped_calc_custom_matrix(
    matrix_name: &str,
    gap_open: i32,
    gap_extend: i32,
    cli: &MatrixCliOverrides,
) -> KarlinBlk {
    // Mode 1: Baked-in table lookup by matrix name + gap params. Failure is
    // silent (error_return NULL in the C code).
    if let Ok(kbp) = karlin_blk_gapped_load_from_tables(gap_open, gap_extend, matrix_name) {
        return kbp;
    }
    // Mode 2: lambda, K, and alpha supplied on the command line.
    if cli.is_set() {
        return KarlinBlk {
            lambda: cli.matrix_lambda,
            k: cli.matrix_k,
            log_k: cli.matrix_k.ln(),
            // Derive H = lambda/alpha for internal consistency.
            h: cli.matrix_lambda / cli.matrix_alpha,
        };
    }
    // Mode 3: no ALP params for this matrix+gap combo and no CLI params.
    KarlinBlk::sentinel()
}
