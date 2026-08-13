//! Ungapped Karlin-Altschul parameter computation, ported from blast_stat.c.
//!
//! This covers the blastn (BLASTNA alphabet) path of
//! Blast_ScoreBlkKbpUngappedCalc: query residue composition, score
//! frequencies, and Lambda/H/K estimation (Blast_KarlinLambdaNR,
//! BlastKarlinLtoH, BlastKarlinLHtoK).

use crate::ncbi_math::{blast_expm1, blast_gcd, blast_powi};

/// BLASTNA alphabet size (blast_encoding.h).
pub const BLASTNA_SIZE: usize = 16;
/// blastna code for 'N' (IUPACNA_TO_BLASTNA\['N'\]).
pub const BLASTNA_N: u8 = 14;
/// blastna code for '-' (IUPACNA_TO_BLASTNA\['-'\]).
pub const BLASTNA_GAP: u8 = 15;

/// BLAST_SCORE_MIN / BLAST_SCORE_MAX (blast_stat.h): INT2_MIN / INT2_MAX.
pub const BLAST_SCORE_MIN: i32 = -32768;
pub const BLAST_SCORE_MAX: i32 = 32767;
/// BLAST_SCORE_RANGE_MAX (blast_stat.c).
pub const BLAST_SCORE_RANGE_MAX: i32 = BLAST_SCORE_MAX - BLAST_SCORE_MIN;

const BLAST_KARLIN_LAMBDA0_DEFAULT: f64 = 0.5;
const BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT: f64 = 1.0e-5;
const BLAST_KARLIN_LAMBDA_ITER_DEFAULT: i32 = 17;
const BLAST_KARLIN_K_SUMLIMIT_DEFAULT: f64 = 0.0001;
const BLAST_KARLIN_K_ITER_MAX: usize = 100;

/// Blast_KarlinBlk (blast_stat.h): Lambda, K, logK, H.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KarlinBlk {
    pub lambda: f64,
    pub k: f64,
    pub log_k: f64,
    pub h: f64,
}

impl KarlinBlk {
    /// The -RMH- Mode 3 sentinel block: no valid statistics.
    pub fn sentinel() -> Self {
        KarlinBlk { lambda: -1.0, k: -1.0, log_k: -1.0, h: -1.0 }
    }

    /// True if this block carries usable statistics (Lambda > 0), the test
    /// used by the -RMH- guards in blast_hits.c / blast_setup.c.
    pub fn is_valid(&self) -> bool {
        self.lambda > 0.0
    }
}

/// Blast_ScoreFreq (blast_stat.h): score frequencies of a scoring system.
pub struct ScoreFreq {
    pub score_min: i32,
    pub score_max: i32,
    pub obs_min: i32,
    pub obs_max: i32,
    pub score_avg: f64,
    /// Probability of each score; index `(score - score_min)`.
    pub sprob: Vec<f64>,
}

impl ScoreFreq {
    pub fn sprob(&self, score: i32) -> f64 {
        self.sprob[(score - self.score_min) as usize]
    }
}

/// BlastScoreChk (blast_stat.c).
fn blast_score_chk(lo: i32, hi: i32) -> Result<(), ()> {
    if lo >= 0 || hi <= 0 || lo < BLAST_SCORE_MIN || hi > BLAST_SCORE_MAX {
        return Err(());
    }
    if hi - lo > BLAST_SCORE_RANGE_MAX {
        return Err(());
    }
    Ok(())
}

/// BlastScoreBlkMaxScoreSet (blast_stat.c): lowest/highest matrix score,
/// ignoring sentinel values at or beyond BLAST_SCORE_MIN/MAX.
pub fn matrix_score_range(matrix: &[[i32; BLASTNA_SIZE]; BLASTNA_SIZE]) -> (i32, i32) {
    let mut lo = BLAST_SCORE_MAX;
    let mut hi = BLAST_SCORE_MIN;
    for row in matrix.iter() {
        for &score in row.iter() {
            if score <= BLAST_SCORE_MIN || score >= BLAST_SCORE_MAX {
                continue;
            }
            if lo > score {
                lo = score;
            }
            if hi < score {
                hi = score;
            }
        }
    }
    if lo < BLAST_SCORE_MIN {
        lo = BLAST_SCORE_MIN;
    }
    if hi > BLAST_SCORE_MAX {
        hi = BLAST_SCORE_MAX;
    }
    (lo, hi)
}

/// Blast_ResFreqStdComp for the nucleotide (blastna) alphabet: uniform 0.25
/// for A, C, G, T (blastna codes 0..3), zero elsewhere (nt_prob +
/// Blast_ResFreqNormalize in blast_stat.c).
pub fn res_freq_std_comp() -> [f64; BLASTNA_SIZE] {
    let mut prob = [0.0; BLASTNA_SIZE];
    for p in prob.iter_mut().take(4) {
        *p = 0.25;
    }
    prob
}

/// Blast_ResFreqString for the blastna alphabet: composition of `query`
/// (blastna-encoded), not counting the ambiguous residues 'N' and '-'
/// (BlastResCompStr + Blast_ResFreqResComp in blast_stat.c).
pub fn res_freq_string(query: &[u8]) -> [f64; BLASTNA_SIZE] {
    let mut comp = [0i64; BLASTNA_SIZE];
    for &c in query {
        // "For megablast, check only the first 4 bits of the sequence values"
        comp[(c & 0x0f) as usize] += 1;
    }
    // Don't count ambig. residues.
    comp[BLASTNA_N as usize] = 0;
    comp[BLASTNA_GAP as usize] = 0;

    let sum: i64 = comp.iter().sum();
    let mut prob = [0.0; BLASTNA_SIZE];
    if sum == 0 {
        return prob;
    }
    for (p, &c) in prob.iter_mut().zip(comp.iter()) {
        *p = c as f64 / sum as f64;
    }
    prob
}

/// BlastScoreFreqCalc (blast_stat.c): score frequencies given the matrix and
/// the residue frequencies of the two sequences (rfp1 = query, rfp2 = db).
pub fn score_freq_calc(
    matrix: &[[i32; BLASTNA_SIZE]; BLASTNA_SIZE],
    rfp1: &[f64; BLASTNA_SIZE],
    rfp2: &[f64; BLASTNA_SIZE],
) -> Result<ScoreFreq, ()> {
    let (lo, hi) = matrix_score_range(matrix);
    blast_score_chk(lo, hi)?;

    let range = (hi - lo + 1) as usize;
    let mut sfp = ScoreFreq {
        score_min: lo,
        score_max: hi,
        obs_min: 0,
        obs_max: 0,
        score_avg: 0.0,
        sprob: vec![0.0; range],
    };

    for (index1, row) in matrix.iter().enumerate() {
        for (index2, &score) in row.iter().enumerate() {
            if score >= lo {
                sfp.sprob[(score - lo) as usize] += rfp1[index1] * rfp2[index2];
            }
        }
    }

    let mut score_sum = 0.0;
    let mut obs_min = BLAST_SCORE_MIN;
    let mut obs_max = BLAST_SCORE_MIN;
    for score in lo..=hi {
        if sfp.sprob[(score - lo) as usize] > 0.0 {
            score_sum += sfp.sprob[(score - lo) as usize];
            obs_max = score;
            if obs_min == BLAST_SCORE_MIN {
                obs_min = score;
            }
        }
    }
    sfp.obs_min = obs_min;
    sfp.obs_max = obs_max;

    let mut score_avg = 0.0;
    if score_sum > 0.0001 || score_sum < -0.0001 {
        for score in obs_min..=obs_max {
            let idx = (score - lo) as usize;
            sfp.sprob[idx] /= score_sum;
            score_avg += score as f64 * sfp.sprob[idx];
        }
    }
    sfp.score_avg = score_avg;

    Ok(sfp)
}

/// NlmKarlinLambdaNR (blast_stat.c): safeguarded Newton iteration for lambda,
/// performed in x = exp(-lambda).
/// `probs` is indexed by score (absolute, via closure over the sfp layout).
#[allow(clippy::too_many_arguments)]
fn nlm_karlin_lambda_nr(
    probs: &dyn Fn(i32) -> f64,
    d: i32,
    low: i32,
    high: i32,
    lambda0: f64,
    tolx: f64,
    itmax: i32,
    max_newton: i32,
) -> f64 {
    let mut a = 0.0f64;
    let mut b = 1.0f64;
    let mut f = 4.0f64; // Larger than any possible value of the poly in [0,1]
    let mut is_newton = false;

    debug_assert!(d > 0);

    let x0 = (-lambda0).exp();
    let mut x = if 0.0 < x0 && x0 < 1.0 { x0 } else { 0.5 };

    for k in 0..itmax {
        let fold = f;
        let was_newton = is_newton;
        is_newton = false;

        // Horner's rule for evaluating a polynomial and its derivative
        let mut g = 0.0;
        f = probs(low);
        let mut i = low + d;
        while i < 0 {
            g = x * g + f;
            f = f * x + probs(i);
            i += d;
        }
        g = x * g + f;
        f = f * x + probs(0) - 1.0;
        let mut i = d;
        while i <= high {
            g = x * g + f;
            f = f * x + probs(i);
            i += d;
        }
        // End Horner's rule

        if f > 0.0 {
            a = x;
        } else if f < 0.0 {
            b = x;
        } else {
            break; // x is an exact solution
        }
        if b - a < 2.0 * a * (1.0 - b) * tolx {
            // The midpoint of the interval converged
            x = (a + b) / 2.0;
            break;
        }

        if k >= max_newton || (was_newton && f.abs() > 0.9 * fold.abs()) || g >= 0.0 {
            // bisect
            x = (a + b) / 2.0;
        } else {
            // try a Newton step
            let p = -f / g;
            let y = x + p;
            if y <= a || y >= b {
                x = (a + b) / 2.0;
            } else {
                is_newton = true;
                x = y;
                if p.abs() < tolx * x * (1.0 - x) {
                    break; // Converged
                }
            }
        }
    }
    -x.ln() / d as f64
}

/// Blast_KarlinLambdaNR (blast_stat.c).
pub fn karlin_lambda_nr(sfp: &ScoreFreq, initial_lambda_guess: f64) -> f64 {
    let low = sfp.obs_min;
    let high = sfp.obs_max;
    if sfp.score_avg >= 0.0 {
        // Expected score must be negative
        return -1.0;
    }
    if blast_score_chk(low, high).is_err() {
        return -1.0;
    }

    // Find greatest common divisor of all scores
    let mut d = -low;
    let mut i = 1;
    while i <= high - low && d > 1 {
        if sfp.sprob(i + low) != 0.0 {
            d = blast_gcd(d, i);
        }
        i += 1;
    }

    nlm_karlin_lambda_nr(
        &|score| sfp.sprob(score),
        d,
        low,
        high,
        initial_lambda_guess,
        BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT,
        20,
        20 + BLAST_KARLIN_LAMBDA_ITER_DEFAULT,
    )
}

/// BlastKarlinLtoH (blast_stat.c): H, the relative entropy.
pub fn karlin_l_to_h(sfp: &ScoreFreq, lambda: f64) -> f64 {
    let low = sfp.obs_min;
    let high = sfp.obs_max;

    if lambda < 0.0 {
        return -1.0;
    }
    if blast_score_chk(low, high).is_err() {
        return -1.0;
    }

    let etonlam = (-lambda).exp();
    let mut sum = low as f64 * sfp.sprob(low);
    for score in (low + 1)..=high {
        sum = score as f64 * sfp.sprob(score) + etonlam * sum;
    }

    let scale = blast_powi(etonlam, high);
    if scale > 0.0 {
        lambda * sum / scale
    } else {
        // Underflow of exp( -lambda * high )
        lambda * (lambda * high as f64 + sum.ln()).exp()
    }
}

/// BlastKarlinLHtoK (blast_stat.c): K from lambda and H.
pub fn karlin_lh_to_k(sfp: &ScoreFreq, lambda: f64, h: f64) -> f64 {
    if lambda <= 0.0 || h <= 0.0 {
        return -1.0;
    }
    // Karlin-Altschul theory works only if the expected score is negative
    if sfp.score_avg >= 0.0 {
        return -1.0;
    }

    let low0 = sfp.obs_min;
    let high0 = sfp.obs_max;
    let range0 = high0 - low0;

    // Greatest common divisor ("delta" in Appendix of PNAS 87)
    let mut divisor = -low0;
    let mut i = 1;
    while i <= range0 && divisor > 1 {
        if sfp.sprob(low0 + i) != 0.0 {
            divisor = blast_gcd(divisor, i);
        }
        i += 1;
    }

    let high = high0 / divisor;
    let low = low0 / divisor;
    let lambda = lambda * divisor as f64;
    let range = high - low;

    let mut first_term_closed_form = h / lambda;
    let exp_minus_lambda = (-lambda).exp();

    if low == -1 && high == 1 {
        let k = (sfp.sprob(low * divisor) - sfp.sprob(high * divisor))
            * (sfp.sprob(low * divisor) - sfp.sprob(high * divisor))
            / sfp.sprob(low * divisor);
        return k;
    }

    if low == -1 || high == 1 {
        if high != 1 {
            let score_avg = sfp.score_avg / divisor as f64;
            first_term_closed_form = (score_avg * score_avg) / first_term_closed_form;
        }
        return first_term_closed_form * (1.0 - exp_minus_lambda);
    }

    let sumlimit = BLAST_KARLIN_K_SUMLIMIT_DEFAULT;
    let iterlimit = BLAST_KARLIN_K_ITER_MAX;

    let mut alignment_score_probabilities =
        vec![0.0f64; iterlimit * range as usize + 1];

    let mut outer_sum = 0.0;
    let mut low_alignment_score: i32 = 0;
    let mut high_alignment_score: i32 = 0;
    let mut inner_sum = 1.0f64;
    alignment_score_probabilities[0] = 1.0;

    // probability array reindexed so that `low0` is at index 0
    let prob_at = |i: i32| sfp.sprob(low0 + i);

    let mut iter_counter = 0usize;
    while iter_counter < iterlimit && inner_sum > sumlimit {
        let mut first = range;
        let mut last = range;
        low_alignment_score += low;
        high_alignment_score += high;

        // dynamic program to compute P(i,j)
        let mut p = (high_alignment_score - low_alignment_score) as isize;
        while p >= 0 {
            let ptr1_start = p - first as isize;
            let ptr1_end = p - last as isize;
            let mut sum = 0.0;
            let mut ptr1 = ptr1_start;
            let mut ptr2 = first;
            while ptr1 >= ptr1_end {
                sum += alignment_score_probabilities[ptr1 as usize] * prob_at(ptr2);
                ptr1 -= 1;
                ptr2 += 1;
            }
            inner_sum = sum;
            if first > 0 {
                first -= 1;
            }
            if p <= range as isize {
                last -= 1;
            }
            alignment_score_probabilities[p as usize] = inner_sum;
            p -= 1;
        }

        // Horner's rule
        let mut idx = 0usize;
        inner_sum = alignment_score_probabilities[idx];
        let mut i = low_alignment_score + 1;
        while i < 0 {
            idx += 1;
            inner_sum = alignment_score_probabilities[idx] + inner_sum * exp_minus_lambda;
            i += 1;
        }
        inner_sum *= exp_minus_lambda;

        while i <= high_alignment_score {
            idx += 1;
            inner_sum += alignment_score_probabilities[idx];
            i += 1;
        }

        iter_counter += 1;
        inner_sum /= iter_counter as f64;
        outer_sum += inner_sum;
    }

    -(-2.0 * outer_sum).exp() / (first_term_closed_form * blast_expm1(-lambda))
}

/// Blast_KarlinBlkUngappedCalc (blast_stat.c): Lambda, H, K from score
/// frequencies. On failure the block carries the C error values
/// (Lambda = H = K = -1, logK = HUGE_VAL) and Err is returned.
pub fn karlin_blk_ungapped_calc(sfp: &ScoreFreq) -> Result<KarlinBlk, KarlinBlk> {
    let err = KarlinBlk { lambda: -1.0, h: -1.0, k: -1.0, log_k: f64::INFINITY };

    let lambda = karlin_lambda_nr(sfp, BLAST_KARLIN_LAMBDA0_DEFAULT);
    if lambda < 0.0 {
        return Err(err);
    }
    let h = karlin_l_to_h(sfp, lambda);
    if h < 0.0 {
        return Err(err);
    }
    let k = karlin_lh_to_k(sfp, lambda, h);
    if k < 0.0 {
        return Err(err);
    }
    Ok(KarlinBlk { lambda, k, log_k: k.ln(), h })
}

/// The per-context blastn flow of Blast_ScoreBlkKbpUngappedCalc: compute the
/// ungapped Karlin block for a blastna-encoded query (one strand/context)
/// against the standard uniform nucleotide composition.
pub fn kbp_ungapped_calc_blastna(
    matrix: &[[i32; BLASTNA_SIZE]; BLASTNA_SIZE],
    query: &[u8],
) -> Result<KarlinBlk, KarlinBlk> {
    let stdrfp = res_freq_std_comp();
    let rfp = res_freq_string(query);
    let sfp = score_freq_calc(matrix, &rfp, &stdrfp).map_err(|_| KarlinBlk {
        lambda: -1.0,
        h: -1.0,
        k: -1.0,
        log_k: f64::INFINITY,
    })?;
    karlin_blk_ungapped_calc(&sfp)
}
