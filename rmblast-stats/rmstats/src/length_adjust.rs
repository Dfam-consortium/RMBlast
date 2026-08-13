//! BLAST_ComputeLengthAdjustment, ported exactly from blast_stat.c:5555.

/// BLAST_ComputeLengthAdjustment: compute the length adjustment ("edge-effect
/// correction") for the effective search space.
///
/// `alpha_d_lambda` is alpha/Lambda (for ungapped statistics, 1/H) and `beta`
/// the corresponding intercept. Returns `(length_adjustment, converged)`;
/// the C function returns 0 iff converged.
pub fn compute_length_adjustment(
    k: f64,
    log_k: f64,
    alpha_d_lambda: f64,
    beta: f64,
    query_length: i32,
    db_length: i64,
    db_num_seqs: i32,
) -> (i32, bool) {
    const K_MAX_ITERATIONS: i32 = 20;
    let m = query_length as f64;
    let n = db_length as f64;
    let big_n = db_num_seqs as f64;

    let mut ell_min: f64 = 0.0;
    let mut ell_max: f64;
    let mut converged = false;
    let mut ell_next: f64 = 0.0;

    // Choose ell_max to be the largest nonnegative value that satisfies
    //    K * (m - ell) * (n - N * ell) > MAX(m,n)
    // Use quadratic formula: 2 c /( - b + sqrt( b*b - 4 * a * c ))
    {
        let a = big_n;
        let mb = m * big_n + n;
        let c = n * m - m.max(n) / k;

        if c < 0.0 {
            return (0, false);
        }
        ell_max = 2.0 * c / (mb + (mb * mb - 4.0 * a * c).sqrt());
    }

    let mut i = 1;
    while i <= K_MAX_ITERATIONS {
        let ell = ell_next;
        let ss = (m - ell) * (n - big_n * ell);
        let ell_bar = alpha_d_lambda * (log_k + ss.ln()) + beta;
        if ell_bar >= ell {
            // ell is no bigger than the true fixed point
            ell_min = ell;
            if ell_bar - ell_min <= 1.0 {
                converged = true;
                break;
            }
            if ell_min == ell_max {
                // There are no more points to check
                break;
            }
        } else {
            // ell is greater than the true fixed point
            ell_max = ell;
        }
        if ell_min <= ell_bar && ell_bar <= ell_max {
            // ell_bar is in range. Accept it
            ell_next = ell_bar;
        } else {
            // ell_bar is not in range. Reject it
            ell_next = if i == 1 { ell_max } else { (ell_min + ell_max) / 2.0 };
        }
        i += 1;
    }

    let mut length_adjustment: i32;
    if converged {
        // If ell_fixed is the (unknown) true fixed point, then we wish to
        // set length_adjustment to floor(ell_fixed). We assume that
        // floor(ell_min) = floor(ell_fixed)
        length_adjustment = ell_min as i32;
        // But verify that ceil(ell_min) != floor(ell_fixed)
        let ell = ell_min.ceil();
        if ell <= ell_max {
            let ss = (m - ell) * (n - big_n * ell);
            if alpha_d_lambda * (log_k + ss.ln()) + beta >= ell {
                // ceil(ell_min) == floor(ell_fixed)
                length_adjustment = ell as i32;
            }
        }
    } else {
        // Use the best value seen so far
        length_adjustment = ell_min as i32;
    }

    (length_adjustment, converged)
}
