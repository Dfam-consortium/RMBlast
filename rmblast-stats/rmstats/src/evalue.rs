//! E-value and bit-score computation, ported from blast_stat.c
//! (BLAST_KarlinStoE_simple, BlastKarlinEtoS_simple, BLAST_GapDecayDivisor,
//! BLAST_Cutoffs) and blast_hits.c (Blast_HSPListGetEvalues,
//! Blast_HSPListGetBitScores) including the -RMH- sentinel guards.

use crate::karlin::KarlinBlk;
use crate::ncbi_math::{blast_powi, NCBIMATH_LN2};

/// BLAST_KarlinStoE_simple (blast_stat.c:4670): E-value from a raw score and
/// the effective search space. Returns -1.0 on invalid parameters.
pub fn karlin_stoe_simple(score: i32, kbp: &KarlinBlk, searchsp: i64) -> f64 {
    if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 {
        return -1.0;
    }
    searchsp as f64 * (-kbp.lambda * score as f64 + kbp.log_k).exp()
}

/// BlastKarlinEtoS_simple (blast_stat.c:4553): score from an expect value.
/// Returns BLAST_SCORE_MIN on invalid parameters.
pub fn karlin_etos_simple(e: f64, kbp: &KarlinBlk, searchsp: i64) -> i32 {
    // Smallest float that might not cause a floating point exception.
    const K_SMALL_FLOAT: f64 = 1.0e-297;

    if kbp.lambda < 0.0 || kbp.k < 0.0 || kbp.h < 0.0 {
        return crate::karlin::BLAST_SCORE_MIN;
    }

    let e = e.max(K_SMALL_FLOAT);
    ((kbp.k * searchsp as f64 / e).ln() / kbp.lambda).ceil() as i32
}

/// BLAST_GapDecayDivisor (blast_stat.c:4592).
pub fn gap_decay_divisor(decayrate: f64, nsegs: u32) -> f64 {
    (1.0 - decayrate) * blast_powi(decayrate, nsegs as i32 - 1)
}

/// BLAST_Cutoffs (blast_stat.c:4603): calculate the cutoff score S and the
/// highest expected score E. `s` and `e` are in/out as in the C interface;
/// returns Err(()) when the Karlin block holds the -1 sentinel.
pub fn blast_cutoffs(
    s: &mut i32,
    e: &mut f64,
    kbp: &KarlinBlk,
    searchsp: i64,
    dodecay: bool,
    gap_decay_rate: f64,
) -> Result<(), ()> {
    if kbp.lambda == -1.0 || kbp.k == -1.0 || kbp.h == -1.0 {
        return Err(());
    }

    let mut es: i32 = 1;
    let esave = *e;
    let mut s_changed = false;

    if *e > 0.0 {
        let mut e_local = *e;
        if dodecay {
            // Invert the adjustment to the e-value that will be applied to
            // compensate for the effect of choosing the best among multiple
            // alignments
            if gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
                e_local *= gap_decay_divisor(gap_decay_rate, 1);
            }
        }
        es = karlin_etos_simple(e_local, kbp, searchsp);
    }
    // Pick the larger cutoff score between the user's choice and that
    // calculated from the value of E.
    if es > *s {
        s_changed = true;
        *s = es;
    }

    // Re-calculate E from the cutoff score, if E going in was too high
    if esave <= 0.0 || !s_changed {
        let mut e_local = karlin_stoe_simple(*s, kbp, searchsp);
        if dodecay && gap_decay_rate > 0.0 && gap_decay_rate < 1.0 {
            // Weight the e-value to compensate for the effect of choosing
            // the best of more than one collection of distinct alignments
            e_local /= gap_decay_divisor(gap_decay_rate, 1);
        }
        *e = e_local;
    }

    Ok(())
}

/// Per-HSP E-value, following Blast_HSPListGetEvalues (blast_hits.c:1816) for
/// the blastn path (no Spouge/Gumbel block, no RPS scaling, gap_decay_rate 0
/// at all blastn call sites so the divisor is 1).
///
/// * The -RMH- sentinel guard: if `kbp` is `None` (invalid context) or
///   Lambda <= 0, the E-value is 1.0.
/// * `round_down` is sbp->round_down (set only by the reward/penalty path;
///   always false for custom matrices): the score is rounded down to an even
///   number for the E-value calculation only.
pub fn hsp_evalue(
    score: i32,
    kbp: Option<&KarlinBlk>,
    eff_searchsp: i64,
    gapped_calculation: bool,
    round_down: bool,
) -> f64 {
    let kbp = match kbp {
        Some(k) if !(k.lambda <= 0.0) => k,
        _ => return 1.0,
    };
    let mut score = score;
    if gapped_calculation && round_down {
        score &= !1;
    }
    karlin_stoe_simple(score, kbp, eff_searchsp)
}

/// Per-HSP bit score, following Blast_HSPListGetBitScores
/// (blast_hits.c:1920). Uses the raw (un-rounded) score. The -RMH- sentinel
/// guard yields 0.0.
pub fn hsp_bit_score(score: i32, kbp: Option<&KarlinBlk>) -> f64 {
    let kbp = match kbp {
        Some(k) if !(k.lambda <= 0.0) => k,
        _ => return 0.0,
    };
    (score as f64 * kbp.lambda - kbp.log_k) / NCBIMATH_LN2
}
