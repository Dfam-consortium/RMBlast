//! Unit tests for the rmstats port.

use rmstats::eff_lengths::{calc_eff_lengths_custom_matrix, SearchSpaceInputs};
use rmstats::evalue::{hsp_bit_score, hsp_evalue, karlin_etos_simple, karlin_stoe_simple};
use rmstats::karlin::{
    karlin_blk_ungapped_calc, karlin_lambda_nr, karlin_lh_to_k, karlin_l_to_h, ScoreFreq,
};
use rmstats::length_adjust::compute_length_adjustment;
use rmstats::nucl_tables::{
    get_nucl_alpha_beta, karlin_blk_nucl_gapped_calc, NuclGappedError,
};
use rmstats::rmblast_tables::{
    karlin_blk_gapped_load_from_tables, kbp_gapped_calc_custom_matrix, MatrixCliOverrides,
    TableLookupError,
};
use rmstats::{KarlinBlk, RmStats};

/// Build a ScoreFreq the way BlastScoreFreqCalc's tail does, from raw
/// (already normalized) per-score probabilities.
fn score_freq(score_min: i32, probs: &[f64]) -> ScoreFreq {
    let score_max = score_min + probs.len() as i32 - 1;
    let mut obs_min = score_min;
    let mut obs_max = score_max;
    for (i, &p) in probs.iter().enumerate() {
        if p > 0.0 {
            obs_min = score_min + i as i32;
            break;
        }
    }
    for (i, &p) in probs.iter().enumerate().rev() {
        if p > 0.0 {
            obs_max = score_min + i as i32;
            break;
        }
    }
    let score_avg = probs
        .iter()
        .enumerate()
        .map(|(i, &p)| (score_min + i as i32) as f64 * p)
        .sum();
    ScoreFreq {
        score_min,
        score_max,
        obs_min,
        obs_max,
        score_avg,
        sprob: probs.to_vec(),
    }
}

/// The worked example in Altschul's comment block in blast_stat.c
/// ("Statistical Significance Parameter Subroutine"): score -2 with
/// probability 0.7, score 0 with probability 0.1, score 3 with probability
/// 0.2 gives lambda = 0.330 and K = 0.154.
#[test]
fn altschul_worked_example() {
    let sfp = score_freq(-2, &[0.7, 0.0, 0.1, 0.0, 0.0, 0.2]);
    let kbp = karlin_blk_ungapped_calc(&sfp).expect("ungapped calc should succeed");
    assert!(
        (kbp.lambda - 0.330).abs() < 5e-4,
        "lambda = {} (expected 0.330)",
        kbp.lambda
    );
    assert!((kbp.k - 0.154).abs() < 5e-4, "K = {} (expected 0.154)", kbp.k);
    assert!(kbp.h > 0.0);
    assert!((kbp.log_k - kbp.k.ln()).abs() < 1e-15);
}

#[test]
fn lambda_h_k_reject_nonnegative_expectation() {
    // Positive average score: KA theory does not apply.
    let sfp = score_freq(-1, &[0.2, 0.0, 0.8]);
    assert!(sfp.score_avg > 0.0);
    assert_eq!(karlin_lambda_nr(&sfp, 0.5), -1.0);
    assert_eq!(karlin_lh_to_k(&sfp, 0.5, 0.5), -1.0);
    assert!(karlin_blk_ungapped_calc(&sfp).is_err());
    // Negative lambda rejected by LtoH.
    assert_eq!(karlin_l_to_h(&sfp, -0.5), -1.0);
}

/// Mode 1: baked table values for 14p35g.matrix at its canonical 29/6 gap
/// costs (rmblast_14p35g_values in blast_stat.c).
#[test]
fn rmblast_table_14p35g() {
    let kbp = karlin_blk_gapped_load_from_tables(29, 6, "14p35g.matrix").unwrap();
    assert_eq!(kbp.lambda, 0.1241860392);
    assert_eq!(kbp.k, 0.2565118357);
    assert_eq!(kbp.h, 0.3349210581);
    assert_eq!(kbp.log_k, 0.2565118357f64.ln());

    // Case-insensitive name match (strcasecmp).
    assert!(karlin_blk_gapped_load_from_tables(29, 6, "14P35G.MATRIX").is_ok());

    // Wrong gap costs: matrix found but combination unsupported.
    assert_eq!(
        karlin_blk_gapped_load_from_tables(20, 5, "14p35g.matrix").unwrap_err(),
        TableLookupError::GapCostsNotSupported
    );
    // Unknown matrix.
    assert_eq!(
        karlin_blk_gapped_load_from_tables(29, 6, "nosuch.matrix").unwrap_err(),
        TableLookupError::MatrixNotFound
    );
}

/// comparison.matrix now has an ALP-fitted entry at 20/5 (new in this port;
/// absent from the C table where it falls through to Mode 3).
#[test]
fn rmblast_table_comparison_matrix() {
    let kbp = karlin_blk_gapped_load_from_tables(20, 5, "comparison.matrix").unwrap();
    assert!((kbp.lambda - 0.0983622300).abs() < 1e-10);
    assert!((kbp.k - 0.0789918084).abs() < 1e-10);
}

#[test]
fn three_mode_hierarchy() {
    let no_cli = MatrixCliOverrides::default();

    // Mode 1: table hit.
    let kbp = kbp_gapped_calc_custom_matrix("25p41g.matrix", 22, 5, &no_cli);
    assert!(kbp.is_valid());

    // Mode 2: table miss, CLI overrides supplied; H = lambda/alpha.
    let cli = MatrixCliOverrides {
        matrix_lambda: 0.1227,
        matrix_k: 0.0919,
        matrix_alpha: 1.5,
        matrix_beta: -4.0,
    };
    let kbp = kbp_gapped_calc_custom_matrix("custom.matrix", 20, 5, &cli);
    assert_eq!(kbp.lambda, 0.1227);
    assert_eq!(kbp.k, 0.0919);
    assert_eq!(kbp.log_k, 0.0919f64.ln());
    assert_eq!(kbp.h, 0.1227 / 1.5);

    // Mode 3: table miss, no CLI params -> sentinel.
    let kbp = kbp_gapped_calc_custom_matrix("custom.matrix", 20, 5, &no_cli);
    assert_eq!(kbp, KarlinBlk::sentinel());
    assert!(!kbp.is_valid());

    // Sentinel semantics downstream (blast_hits.c -RMH- guards).
    assert_eq!(hsp_evalue(1000, Some(&kbp), 1_000_000, true, false), 1.0);
    assert_eq!(hsp_bit_score(1000, Some(&kbp)), 0.0);
    assert_eq!(hsp_evalue(1000, None, 1_000_000, true, false), 1.0);
    assert_eq!(hsp_bit_score(1000, None), 0.0);
}

/// Blast_KarlinBlkNuclGappedCalc against the blastn default scoring
/// (reward 2, penalty -3, gap 5/2) and megablast linear values.
#[test]
fn nucl_gapped_calc() {
    let ungap = KarlinBlk { lambda: 1.28, k: 0.46, log_k: 0.46f64.ln(), h: 0.85 };

    // blastn defaults: 2/-3 with gap costs 5/2 -> {5, 2, 0.625, 0.41, 0.78, ...}
    let (kbp, round_down) = karlin_blk_nucl_gapped_calc(5, 2, 2, -3, &ungap).unwrap();
    assert_eq!(kbp.lambda, 0.625);
    assert_eq!(kbp.k, 0.41);
    assert_eq!(kbp.h, 0.78);
    assert!(round_down, "2/-3 sets round_down");

    // Megablast linear: 1/-3 with gap costs 0/0 -> linear row {1.374, 0.711, 1.31}
    let (kbp, round_down) = karlin_blk_nucl_gapped_calc(0, 0, 1, -3, &ungap).unwrap();
    assert_eq!(kbp.lambda, 1.374);
    assert_eq!(kbp.k, 0.711);
    assert!(!round_down);

    // Infinite gap regime: 1/-3, gap costs >= (2,2) but not tabulated -> ungapped copy.
    let (kbp, _) = karlin_blk_nucl_gapped_calc(20, 5, 1, -3, &ungap).unwrap();
    assert_eq!(kbp.lambda, ungap.lambda);
    assert_eq!(kbp.k, ungap.k);

    // gcd adjustment: 2/-6 reduces to 1/-3 with gap costs doubled; the 1/-3
    // row {2,2,1.37,0.70,1.2} becomes {4,4,0.685,0.70,1.2}.
    let (kbp, _) = karlin_blk_nucl_gapped_calc(4, 4, 2, -6, &ungap).unwrap();
    assert_eq!(kbp.lambda, 1.37 / 2.0);
    assert_eq!(kbp.k, 0.70);
    assert_eq!(kbp.h, 1.2);

    // Unsupported pair.
    assert_eq!(
        karlin_blk_nucl_gapped_calc(5, 2, 7, -3, &ungap).unwrap_err(),
        NuclGappedError::UnsupportedRewardPenalty
    );
    // Unsupported gap costs (1/-3 with 0/1: below the infinite regime,
    // not tabulated).
    assert_eq!(
        karlin_blk_nucl_gapped_calc(0, 1, 1, -3, &ungap).unwrap_err(),
        NuclGappedError::UnsupportedGapCosts
    );
}

/// Blast_GetNuclAlphaBeta: tabulated, linear, and ungapped-fallback modes.
#[test]
fn nucl_alpha_beta() {
    let ungap = KarlinBlk { lambda: 1.33, k: 0.62, log_k: 0.62f64.ln(), h: 1.12 };

    // Tabulated gapped: 1/-2 with 3/1 -> alpha 1.3, beta -1.
    let (alpha, beta) = get_nucl_alpha_beta(1, -2, 3, 1, &ungap, true).unwrap();
    assert_eq!(alpha, 1.3);
    assert_eq!(beta, -1.0);

    // Linear row: 1/-2 with 0/0 -> alpha 1.5, beta -2.
    let (alpha, beta) = get_nucl_alpha_beta(1, -2, 0, 0, &ungap, true).unwrap();
    assert_eq!(alpha, 1.5);
    assert_eq!(beta, -2.0);

    // Fallback (gap costs not in the 1/-3 table): alpha = Lambda/H of the
    // ungapped block, beta = 0 (s_GetUngappedBeta(1,-3) = 0).
    let (alpha, beta) = get_nucl_alpha_beta(1, -3, 20, 5, &ungap, true).unwrap();
    assert_eq!(alpha, ungap.lambda / ungap.h);
    assert_eq!(beta, 0.0);

    // Ungapped search: fallback values, with beta = -2 for 1/-1.
    let (alpha, beta) = get_nucl_alpha_beta(1, -1, 5, 2, &ungap, false).unwrap();
    assert_eq!(alpha, ungap.lambda / ungap.h);
    assert_eq!(beta, -2.0);

    // Unsupported pair errors out.
    assert!(get_nucl_alpha_beta(7, -3, 5, 2, &ungap, true).is_err());
}

/// BLAST_ComputeLengthAdjustment fixed-point behavior: the returned ell
/// satisfies the defining inequality and its successor does not (for a
/// converged gapped case), i.e. ell = floor(ell_fixed).
#[test]
fn length_adjustment_fixed_point() {
    // 25p41g.matrix parameters, a chr20-scale query vs the longlib database.
    let kbp = karlin_blk_gapped_load_from_tables(22, 5, "25p41g.matrix").unwrap();
    let alpha = kbp.lambda / kbp.h;
    let beta = 0.0;
    let (m, n, nseq) = (64_444_167i64, 1_838_298i64, 1_580i32);

    let (ell, converged) = compute_length_adjustment(
        kbp.k,
        kbp.log_k,
        alpha / kbp.lambda,
        beta,
        m as i32,
        n,
        nseq,
    );
    assert!(converged);
    assert!(ell > 0);

    let f = |l: f64| {
        let ss = (m as f64 - l) * (n as f64 - nseq as f64 * l);
        (alpha / kbp.lambda) * (kbp.log_k + ss.ln()) + beta
    };
    // ell is at most the fixed point; ell+1 is beyond it.
    assert!(f(ell as f64) >= ell as f64);
    assert!(f((ell + 1) as f64) < (ell + 1) as f64 + 1.0);

    // Degenerate case: c < 0 (search space too small for K) -> (0, false).
    let (ell, converged) =
        compute_length_adjustment(1e-9, (1e-9f64).ln(), 2.0, 0.0, 10, 10, 1);
    assert_eq!(ell, 0);
    assert!(!converged);
}

/// E-value / bit-score formulas.
#[test]
fn evalue_bitscore_formulas() {
    let kbp = karlin_blk_gapped_load_from_tables(25, 5, "20p41g.matrix").unwrap();
    let searchsp: i64 = 1_000_000_000;

    let s = 300;
    let e = karlin_stoe_simple(s, &kbp, searchsp);
    let expected = searchsp as f64 * (-kbp.lambda * s as f64 + kbp.log_k).exp();
    assert_eq!(e, expected);

    // Round trip through EtoS: score for that E-value is <= s and within one
    // lattice step.
    let s2 = karlin_etos_simple(e, &kbp, searchsp);
    assert!((s2 - s).abs() <= 1, "s2 = {s2}");

    // Bit score uses the raw score.
    let bits = hsp_bit_score(s, Some(&kbp));
    assert!(
        (bits - (s as f64 * kbp.lambda - kbp.log_k) / std::f64::consts::LN_2).abs() < 1e-12
    );

    // round_down applies to the E-value only (SB-2303).
    let e_odd = hsp_evalue(301, Some(&kbp), searchsp, true, true);
    let e_even = hsp_evalue(300, Some(&kbp), searchsp, true, false);
    assert_eq!(e_odd, e_even);
    let b_odd = hsp_bit_score(301, Some(&kbp));
    assert!(b_odd > hsp_bit_score(300, Some(&kbp)));

    // Invalid kbp in the low-level formula returns -1.
    let bad = KarlinBlk { lambda: -1.0, k: -1.0, log_k: -1.0, h: -1.0 };
    assert_eq!(karlin_stoe_simple(300, &bad, searchsp), -1.0);
}

/// End-to-end RmStats: Mode 1 known matrix; effective search space
/// consistency with the parts.
#[test]
fn rmstats_end_to_end() {
    let cli = MatrixCliOverrides::default();
    let stats = RmStats::new_custom_matrix(
        "20p41g.matrix",
        25,
        5,
        &cli,
        1_000_000,   // query length
        1_838_298,   // db length
        1_580,       // db seqs
        0,           // no searchsp override
        None,
    );
    assert!(stats.has_stats());
    assert!(stats.eff.length_adjustment > 0);
    let expected_ss = (1_838_298i64 - 1_580i64 * stats.eff.length_adjustment as i64)
        * (1_000_000i64 - stats.eff.length_adjustment as i64);
    assert_eq!(stats.eff.eff_searchsp, expected_ss);

    let e = stats.evalue(250);
    let b = stats.bit_score(250);
    assert!(e > 0.0 && e.is_finite());
    assert!(b > 0.0 && b.is_finite());

    // Sentinel path: unknown matrix without CLI params.
    let stats = RmStats::new_custom_matrix(
        "comparison.matrix",
        20,
        6, // gap costs not fitted -> Mode 3
        &cli,
        1_000_000,
        1_838_298,
        1_580,
        0,
        None,
    );
    assert!(!stats.has_stats());
    assert_eq!(stats.evalue(250), 1.0);
    assert_eq!(stats.bit_score(250), 0.0);
    assert_eq!(stats.eff.eff_searchsp, 0);
    assert_eq!(stats.eff.length_adjustment, 0);
}

/// The sentinel skip in BLAST_CalcEffLengths keeps a user -searchsp override
/// even when statistics are invalid.
#[test]
fn eff_lengths_sentinel_keeps_override() {
    let inputs = SearchSpaceInputs {
        query_length: 1000,
        db_length: 100_000,
        db_num_seqs: 10,
        eff_searchsp_override: 424242,
    };
    let sentinel = KarlinBlk::sentinel();
    let eff = calc_eff_lengths_custom_matrix(
        &inputs,
        &sentinel,
        None,
        &MatrixCliOverrides::default(),
        20,
        5,
        true,
    );
    assert_eq!(eff.eff_searchsp, 424242);
    assert_eq!(eff.length_adjustment, 0);
}
