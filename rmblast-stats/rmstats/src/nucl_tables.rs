//! Reward/penalty-based nucleotide KA parameter tables and lookups, ported
//! from blast_stat.c: s_GetNuclValuesArray, s_SplitArrayOf8,
//! s_AdjustGapParametersByGcd, Blast_KarlinBlkNuclGappedCalc,
//! Blast_GetNuclAlphaBeta, s_GetUngappedBeta.

use crate::karlin::KarlinBlk;
use crate::ncbi_math::{blast_gcd, blast_nint};
use crate::tables_data::*;

/// array_of_8 column layout (blast_stat.c):
/// 0: gap open, 1: gap extend, 2: Lambda, 3: K, 4: H, 5: alpha, 6: beta, 7: theta
const K_GAP_OPEN_INDEX: usize = 0;
const K_GAP_EXT_INDEX: usize = 1;
const K_LAMBDA_INDEX: usize = 2;
const K_K_INDEX: usize = 3;
const K_H_INDEX: usize = 4;
const K_ALPHA_INDEX: usize = 5;
const K_BETA_INDEX: usize = 6;

/// Result of s_GetNuclValuesArray.
pub struct NuclValues {
    /// The standard (typically affine) values; owned copies, possibly
    /// gcd-adjusted.
    pub normal: Vec<[f64; 8]>,
    /// The megablast linear (non-affine) values, if this reward/penalty pair
    /// has them.
    pub linear: Option<[f64; 8]>,
    pub gap_open_max: i32,
    pub gap_extend_max: i32,
    pub round_down: bool,
}

/// s_SplitArrayOf8: split off a leading {0,0,...} linear row.
fn split_array_of_8(input: &'static [[f64; 8]]) -> (Vec<[f64; 8]>, Option<[f64; 8]>) {
    if input[0][0] == 0.0 && input[0][1] == 0.0 {
        (input[1..].to_vec(), Some(input[0]))
    } else {
        (input.to_vec(), None)
    }
}

/// s_AdjustGapParametersByGcd: scale gap costs up and Lambda/alpha down by
/// the reward/penalty gcd.
fn adjust_gap_parameters_by_gcd(
    normal: &mut [[f64; 8]],
    linear: &mut Option<[f64; 8]>,
    gap_existence_max: &mut i32,
    gap_extend_max: &mut i32,
    divisor: i32,
) {
    if divisor == 1 {
        return;
    }
    let d = divisor as f64;
    *gap_existence_max *= divisor;
    *gap_extend_max *= divisor;

    for row in normal.iter_mut() {
        row[0] *= d;
        row[1] *= d;
        row[2] /= d;
        row[5] /= d;
    }
    if let Some(row) = linear.as_mut() {
        row[0] *= d;
        row[1] *= d;
        row[2] /= d;
        row[5] /= d;
    }
}

/// s_GetNuclValuesArray (blast_stat.c). Err(()) for an unsupported
/// reward/penalty pair ("Substitution scores %d and %d are not supported").
pub fn get_nucl_values_array(reward: i32, penalty: i32) -> Result<NuclValues, ()> {
    let divisor = blast_gcd(reward, penalty);
    let (mut reward, mut penalty) = (reward, penalty);
    if divisor != 1 {
        reward /= divisor;
        penalty /= divisor;
    }

    let (table, gap_open_max, gap_extend_max, round_down): (&'static [[f64; 8]], i32, i32, bool) =
        match (reward, penalty) {
            (1, -5) => (BLASTN_VALUES_1_5, 3, 3, false),
            (1, -4) => (BLASTN_VALUES_1_4, 2, 2, false),
            (2, -7) => (BLASTN_VALUES_2_7, 4, 4, true),
            (1, -3) => (BLASTN_VALUES_1_3, 2, 2, false),
            (2, -5) => (BLASTN_VALUES_2_5, 4, 4, true),
            (1, -2) => (BLASTN_VALUES_1_2, 2, 2, false),
            (2, -3) => (BLASTN_VALUES_2_3, 6, 4, true),
            (3, -4) => (BLASTN_VALUES_3_4, 6, 3, true),
            (1, -1) => (BLASTN_VALUES_1_1, 4, 2, false),
            (3, -2) => (BLASTN_VALUES_3_2, 5, 5, false),
            (4, -5) => (BLASTN_VALUES_4_5, 12, 8, false),
            (5, -4) => (BLASTN_VALUES_5_4, 25, 10, false),
            _ => return Err(()),
        };

    let (mut normal, mut linear) = split_array_of_8(table);
    let mut gap_open_max = gap_open_max;
    let mut gap_extend_max = gap_extend_max;
    adjust_gap_parameters_by_gcd(
        &mut normal,
        &mut linear,
        &mut gap_open_max,
        &mut gap_extend_max,
        divisor,
    );

    Ok(NuclValues { normal, linear, gap_open_max, gap_extend_max, round_down })
}

/// Errors from Blast_KarlinBlkNuclGappedCalc.
#[derive(Debug, PartialEq, Eq)]
pub enum NuclGappedError {
    /// status -1: "Substitution scores %d and %d are not supported"
    UnsupportedRewardPenalty,
    /// status 1: "Gap existence and extension values %ld and %ld are not
    /// supported for substitution scores %ld and %ld"
    UnsupportedGapCosts,
}

/// Blast_KarlinBlkNuclGappedCalc (blast_stat.c): gapped KA parameters for a
/// reward/penalty scoring system. `kbp_ungap` is used for the infinite gap
/// cost regime. Returns (KarlinBlk, round_down).
pub fn karlin_blk_nucl_gapped_calc(
    gap_open: i32,
    gap_extend: i32,
    reward: i32,
    penalty: i32,
    kbp_ungap: &KarlinBlk,
) -> Result<(KarlinBlk, bool), NuclGappedError> {
    let values = get_nucl_values_array(reward, penalty)
        .map_err(|_| NuclGappedError::UnsupportedRewardPenalty)?;

    if gap_open == 0 && gap_extend == 0 {
        if let Some(linear) = values.linear {
            let k = linear[K_K_INDEX];
            return Ok((
                KarlinBlk {
                    lambda: linear[K_LAMBDA_INDEX],
                    k,
                    log_k: k.ln(),
                    h: linear[K_H_INDEX],
                },
                values.round_down,
            ));
        }
    }

    for row in &values.normal {
        if blast_nint(row[K_GAP_OPEN_INDEX]) == gap_open as i64
            && blast_nint(row[K_GAP_EXT_INDEX]) == gap_extend as i64
        {
            let k = row[K_K_INDEX];
            return Ok((
                KarlinBlk {
                    lambda: row[K_LAMBDA_INDEX],
                    k,
                    log_k: k.ln(),
                    h: row[K_H_INDEX],
                },
                values.round_down,
            ));
        }
    }

    // If gap costs are larger than maximal provided in tables, copy the
    // values from the ungapped Karlin block.
    if gap_open >= values.gap_open_max && gap_extend >= values.gap_extend_max {
        return Ok((*kbp_ungap, values.round_down));
    }

    Err(NuclGappedError::UnsupportedGapCosts)
}

/// s_GetUngappedBeta (blast_stat.c).
fn get_ungapped_beta(reward: i32, penalty: i32) -> f64 {
    if (reward == 1 && penalty == -1) || (reward == 2 && penalty == -3) {
        -2.0
    } else {
        0.0
    }
}

/// Blast_GetNuclAlphaBeta (blast_stat.c). `kbp` is the ungapped Karlin block
/// used for the fallback alpha = Lambda/H. Err(()) mirrors the status return
/// for an unsupported reward/penalty pair (in which case the C caller's
/// alpha/beta retain their previous values).
pub fn get_nucl_alpha_beta(
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    kbp: &KarlinBlk,
    gapped_calculation: bool,
) -> Result<(f64, f64), ()> {
    let values = get_nucl_values_array(reward, penalty)?;

    if gapped_calculation && !values.normal.is_empty() {
        if gap_open == 0 && gap_extend == 0 {
            if let Some(linear) = values.linear {
                return Ok((linear[K_ALPHA_INDEX], linear[K_BETA_INDEX]));
            }
        } else {
            for row in &values.normal {
                if row[K_GAP_OPEN_INDEX] == gap_open as f64
                    && row[K_GAP_EXT_INDEX] == gap_extend as f64
                {
                    return Ok((row[K_ALPHA_INDEX], row[K_BETA_INDEX]));
                }
            }
        }
    }

    // If input values not found in tables, or if this is an ungapped search,
    // return the ungapped values of alpha and beta.
    Ok((kbp.lambda / kbp.h, get_ungapped_beta(reward, penalty)))
}
