//! Tests for the "alp-fit" feature: startup-time Gumbel fitting via the ALP
//! library. Run with `cargo test --features alp-fit`.
#![cfg(feature = "alp-fit")]

use rmstats::alp::{fit_gumbel_gapped, AlpFitOptions};
use rmstats::rmblast_tables::karlin_blk_gapped_load_from_tables;

/// ALP uses process-global RNG state, so the two fits below must not run in
/// parallel test threads; both live in this single test function.
#[test]
fn ffi_fits() {
    ffi_fit_reproduces_14p35g();
    ffi_fit_deterministic_mode();
}

/// Fitting the 4x4 A/C/G/T core of 14p35g.matrix at 29/6 with the same
/// settings used for the baked table (ALP defaults, seed 1) must reproduce
/// rmblast_14p35g_values exactly — the fit is deterministic for a fixed
/// seed.
fn ffi_fit_reproduces_14p35g() {
    // 4x4 core of /usr/local/RepeatMasker/Matrices/ncbi/nt/14p35g.matrix
    // in A,C,G,T order; # FREQS A 0.325 C 0.175 G 0.175 T 0.325.
    let scores: Vec<Vec<i64>> = vec![
        vec![8, -17, -7, -21],
        vec![-18, 12, -16, -10],
        vec![-10, -16, 12, -18],
        vec![-21, -7, -17, 8],
    ];
    let freqs = [0.325, 0.175, 0.175, 0.325];

    let fit = fit_gumbel_gapped(&scores, &freqs, &freqs, 29, 6, 29, 6, &AlpFitOptions::default())
        .expect("ALP fit failed");

    let baked = karlin_blk_gapped_load_from_tables(29, 6, "14p35g.matrix").unwrap();
    let kbp = fit.to_karlin_blk();
    assert!(
        (kbp.lambda - baked.lambda).abs() < 5e-10,
        "lambda {} vs baked {}",
        kbp.lambda,
        baked.lambda
    );
    assert!((kbp.k - baked.k).abs() < 5e-10, "K {} vs baked {}", kbp.k, baked.k);
    assert!((kbp.h - baked.h).abs() < 5e-10, "H {} vs baked {}", kbp.h, baked.h);
    // NCBI-mapped beta = b_J + b_I from the same fit.
    assert!((fit.ncbi_beta() - -4.9556787558).abs() < 5e-10);
    // Error estimates present.
    assert!(fit.raw.lambda_error > 0.0 && fit.raw.k_error > 0.0);
}

/// The LAST-style deterministic mode also converges and lands within the
/// fit's own error bars of the baked values.
fn ffi_fit_deterministic_mode() {
    let scores: Vec<Vec<i64>> = vec![
        vec![8, -17, -7, -21],
        vec![-18, 12, -16, -10],
        vec![-10, -16, 12, -18],
        vec![-21, -7, -17, 8],
    ];
    let freqs = [0.325, 0.175, 0.175, 0.325];

    let opts = AlpFitOptions { deterministic: true, max_time: 60.0, ..Default::default() };
    let fit = fit_gumbel_gapped(&scores, &freqs, &freqs, 29, 6, 29, 6, &opts)
        .expect("deterministic ALP fit failed");
    let baked = karlin_blk_gapped_load_from_tables(29, 6, "14p35g.matrix").unwrap();
    let kbp = fit.to_karlin_blk();
    // Different sampling schedule -> not bit-identical, but must agree to
    // within a few combined error bars.
    assert!(
        (kbp.lambda - baked.lambda).abs() < 6.0 * fit.raw.lambda_error,
        "lambda {} vs baked {} (err {})",
        kbp.lambda,
        baked.lambda,
        fit.raw.lambda_error
    );
    assert!(
        (kbp.k - baked.k).abs() < 6.0 * fit.raw.k_error,
        "K {} vs baked {} (err {})",
        kbp.k,
        baked.k,
        fit.raw.k_error
    );
}
