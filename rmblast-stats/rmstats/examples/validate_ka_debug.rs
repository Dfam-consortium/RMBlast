//! Cross-validate rmstats against NCBI rmblastn BLAST_KA_DEBUG output.
//!
//! Usage:
//!   cargo run --example validate_ka_debug -- \
//!       <matrix> <gapopen> <gapextend> <query_len> <db_len> <db_num_seqs> \
//!       <ka_debug_file>
//!
//! The KA_DEBUG lines (stderr of rmblastn with BLAST_KA_DEBUG=1) look like:
//!   KA_DEBUG: raw_score=1873  lambda=0.1088138918  K=0.1409300313
//!   logK=-1.9594917437  H=0.1778018484  alpha(lambda/H)=0.6119952789
//!   evalue=2.7844e-78  bit_score=296.8603
//! where evalue was printed with %.6g and bit_score with %.4f.

use std::collections::HashMap;

use rmstats::{MatrixCliOverrides, RmStats};

fn field(line: &str, key: &str) -> Option<f64> {
    let start = line.find(key)? + key.len();
    let rest = &line[start..];
    let end = rest.find(char::is_whitespace).unwrap_or(rest.len());
    rest[..end].parse().ok()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 8 && args.len() != 12 {
        eprintln!(
            "usage: validate_ka_debug <matrix> <gapopen> <gapextend> \
             <query_len> <db_len> <db_num_seqs> <ka_debug_file> \
             [<matrix_lambda> <matrix_k> <matrix_alpha> <matrix_beta>]"
        );
        std::process::exit(2);
    }
    let matrix = &args[1];
    let gap_open: i32 = args[2].parse().unwrap();
    let gap_extend: i32 = args[3].parse().unwrap();
    let query_len: i32 = args[4].parse().unwrap();
    let db_len: i64 = args[5].parse().unwrap();
    let db_num_seqs: i32 = args[6].parse().unwrap();
    let path = &args[7];
    let cli = if args.len() == 12 {
        MatrixCliOverrides {
            matrix_lambda: args[8].parse().unwrap(),
            matrix_k: args[9].parse().unwrap(),
            matrix_alpha: args[10].parse().unwrap(),
            matrix_beta: args[11].parse().unwrap(),
        }
    } else {
        MatrixCliOverrides::default()
    };

    let stats = RmStats::new_custom_matrix(
        matrix,
        gap_open,
        gap_extend,
        &cli,
        query_len,
        db_len,
        db_num_seqs,
        0,
        None,
    );
    println!(
        "rmstats: lambda={:.10} K={:.10} logK={:.10} H={:.10} \
         length_adjustment={} eff_searchsp={}",
        stats.kbp_gap.lambda,
        stats.kbp_gap.k,
        stats.kbp_gap.log_k,
        stats.kbp_gap.h,
        stats.eff.length_adjustment,
        stats.eff.eff_searchsp
    );

    let text = std::fs::read_to_string(path).unwrap();
    let mut n = 0usize;
    let mut params_checked = false;
    let mut mismatches = 0usize;
    // score -> (worst evalue rel diff, worst bits abs diff)
    let mut worst: HashMap<i32, (f64, f64)> = HashMap::new();
    let (mut worst_e, mut worst_b) = (0.0f64, 0.0f64);

    for line in text.lines() {
        if !line.starts_with("KA_DEBUG:") {
            continue;
        }
        n += 1;
        let score = field(line, "raw_score=").unwrap() as i32;
        let lambda = field(line, "lambda=").unwrap();
        let k = field(line, "K=").unwrap();
        let log_k = field(line, "logK=").unwrap();
        let h = field(line, "H=").unwrap();
        let evalue = field(line, "evalue=").unwrap();
        let bits = field(line, "bit_score=").unwrap();

        if !params_checked {
            params_checked = true;
            assert!(
                (lambda - stats.kbp_gap.lambda).abs() < 5e-11
                    && (k - stats.kbp_gap.k).abs() < 5e-11
                    && (log_k - stats.kbp_gap.log_k).abs() < 5e-11
                    && (h - stats.kbp_gap.h).abs() < 5e-11,
                "KA parameter mismatch: NCBI lambda={lambda} K={k} logK={log_k} H={h}"
            );
            println!("KA parameters match NCBI (lambda, K, logK, H).");
        }

        let my_e = stats.evalue(score);
        let my_b = stats.bit_score(score);

        // NCBI printed evalue with %.6g (6 significant digits) and
        // bit_score with %.4f.
        let e_rel = if evalue != 0.0 {
            ((my_e - evalue) / evalue).abs()
        } else {
            my_e.abs()
        };
        let b_abs = (my_b - bits).abs();
        if e_rel > 1e-5 || b_abs > 1.5e-4 {
            mismatches += 1;
            let w = worst.entry(score).or_insert((0.0, 0.0));
            w.0 = w.0.max(e_rel);
            w.1 = w.1.max(b_abs);
            if mismatches <= 5 {
                println!(
                    "MISMATCH score={score}: ncbi evalue={evalue:e} mine={my_e:e} \
                     (rel {e_rel:.2e}); ncbi bits={bits} mine={my_b:.4} (abs {b_abs:.2e})"
                );
            }
        }
        worst_e = worst_e.max(e_rel);
        worst_b = worst_b.max(b_abs);
    }

    println!(
        "{n} KA_DEBUG records; {mismatches} mismatches; \
         worst evalue rel diff {worst_e:.3e} (tol 1e-5, %.6g print), \
         worst bits abs diff {worst_b:.3e} (tol 1.5e-4, %.4f print)"
    );
    if mismatches > 0 {
        std::process::exit(1);
    }
    println!("PASS");
}
