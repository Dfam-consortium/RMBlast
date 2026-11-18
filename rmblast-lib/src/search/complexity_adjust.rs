//! Complexity-adjusted scoring — port of blast_traceback.c lines 522-574.
//!
//! Formula: adj_score = raw_score + t_sum / lambda   (+ 0.999 for rounding)
//! where t_sum corrects for how far the aligned query composition deviates
//! from the background base frequencies in the scoring matrix.

use crate::hits::{EditOp, EditScript};
use crate::matrix::ScoreMatrix;
use crate::encoding::{BLASTNA_COMPLEMENT, BLASTNA_SIZE};

/// Collect the BLASTNA bases from query positions that appear in aligned
/// (Sub or GapInSubject) positions.  NCBI's complexity adjustment counts
/// Sub positions for the composition; it skips GapInSubject positions
/// (query has a base but subject has a gap) as these aren't "aligned pairs".
///
/// `query` = BLASTNA slice (with sentinels at index 0 and end).
/// `q_start` = 0-based offset of the first aligned query base (into the raw bases,
///             not counting the leading sentinel).
pub fn collect_query_aligned_bases(
    query: &[u8],
    q_start: u32,
    edit_script: &EditScript,
) -> Vec<u8> {
    let mut bases = Vec::new();
    let mut q_pos = q_start as usize;

    for &(op, count) in &edit_script.ops {
        let n = count as usize;
        match op {
            EditOp::Sub => {
                for i in 0..n {
                    let idx = 1 + q_pos + i; // +1 for leading sentinel
                    if idx < query.len() {
                        let b = query[idx];
                        if b != 15 {
                            bases.push(b);
                        }
                    }
                }
                q_pos += n;
            }
            EditOp::GapInSubject => {
                // Query advances without subject — skip these positions for composition.
                q_pos += n;
            }
            EditOp::GapInQuery => {
                // Subject advances, query stays — nothing to count.
            }
        }
    }
    bases
}

/// Compute the complexity-adjusted score.
///
/// `query_aligned_bases` — BLASTNA bases from the aligned (Sub-position) query.
pub fn compute_adjusted_score(
    raw_score: i32,
    query_aligned_bases: &[u8],
    matrix: &ScoreMatrix,
) -> f64 {
    if matrix.lambda <= 0.0 || query_aligned_bases.is_empty() {
        return raw_score as f64;
    }

    let mut counts = [0u32; BLASTNA_SIZE];
    for &b in query_aligned_bases {
        counts[(b & 15) as usize] += 1;
    }

    let total_count: u32 = counts.iter().sum();
    if total_count == 0 {
        return raw_score as f64;
    }

    let mut t_factor = 0.0f64;
    let mut t_sum = 0.0f64;
    let mut t_counts = 0.0f64; // matches NCBI's t_counts: only bases with f>0 && ln(f)!=0

    for i in 0..BLASTNA_SIZE {
        let c = counts[i];
        let f = matrix.freqs[i];
        if c > 0 && f > 0.0 && f.ln() != 0.0 {
            let cf = c as f64;
            t_factor += cf * cf.ln();
            t_sum += cf * f.ln();
            t_counts += cf;
        }
    }

    if t_counts == 0.0 {
        // NCBI computes t_counts * log(t_counts) = 0 * log(0) = NaN here, and
        // the NaN score then fails the cutoff, so the hit is dropped.
        return f64::NAN;
    }

    t_factor -= t_counts * t_counts.ln();
    t_sum -= t_factor;

    let adj = raw_score as f64 + t_sum / matrix.lambda + 0.999;
    if adj < 0.0 { 0.0 } else { adj }
}

/// Apply complexity adjustment; return None if score falls below cutoff.
///
/// `query`, `q_start` and `edit_script` are in forward-query terms.  For a
/// minus-strand hit NCBI counts the composition of its context query, the
/// reverse complement, so the aligned bases are complemented: the same
/// positions read from the other strand.  (Matrices with complement-symmetric
/// FREQS, like comparison.matrix, give the same score either way.)
pub fn apply_complexity_adjust(
    raw_score: i32,
    query: &[u8],
    q_start: u32,
    edit_script: &EditScript,
    matrix: &ScoreMatrix,
    cutoff: i32,
    is_minus: bool,
) -> Option<i32> {
    let mut aligned_bases = collect_query_aligned_bases(query, q_start, edit_script);
    if is_minus {
        for b in &mut aligned_bases {
            *b = BLASTNA_COMPLEMENT[(*b & 15) as usize];
        }
    }
    let adj = compute_adjusted_score(raw_score, &aligned_bases, matrix);
    if adj.is_nan() {
        return None;
    }
    let rounded = adj as i32;
    if rounded < cutoff { None } else { Some(rounded) }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::ScoreMatrix;
    use std::io::Cursor;

    const MAT: &str = r"
# FREQS A 0.255 C 0.245 G 0.245 T 0.255
   A   C   G   T
A  9  -6  -5 -13
C -6   9 -13  -5
G -5 -13   9  -6
T -13  -5  -6   9
";

    fn make_matrix() -> ScoreMatrix {
        ScoreMatrix::from_reader("test", Cursor::new(MAT)).unwrap()
    }

    // Expected values below are NCBI's formula (s_ComputeComplexityAdjustedScore
    // as RMBlast applies it) written out for each composition:
    //   adj = raw + (sum c ln f - (sum c ln c - N ln N)) / lambda + 0.999,
    // clamped at 0.  This code matches NCBI on every -complexity_adjust
    // suite combination and on the asymmetric-FREQS fixture (tests/cli.rs).

    #[test]
    fn test_balanced_composition_unchanged() {
        let m = make_matrix();
        let bases: Vec<u8> = (0..16).flat_map(|_| vec![0u8, 1, 2, 3]).collect(); // ACGT×4
        let adj = compute_adjusted_score(100, &bases, &m);
        // 16 of each base: sum c ln c - N ln N = 64 ln(1/4)
        let t_sum = 32.0 * (0.255f64.ln() + 0.245f64.ln()) - 64.0 * 0.25f64.ln();
        let want = 100.0 + t_sum / m.lambda + 0.999;
        assert!((adj - want).abs() < 1e-9, "adj={adj} want={want}");
        assert!((adj - 100.8808).abs() < 1e-3, "adj={adj}");
    }

    #[test]
    fn test_biased_composition_reduced() {
        let m = make_matrix();
        let bases = vec![0u8; 64]; // all A: no composition credit, 64 ln f(A)
        let want = 1000.0 + 64.0 * 0.255f64.ln() / m.lambda + 0.999;
        let adj = compute_adjusted_score(1000, &bases, &m);
        assert!((adj - want).abs() < 1e-9, "adj={adj} want={want}");
        assert!(adj > 0.0 && adj < 1000.0);
        // the same composition on a low score clamps to 0
        assert_eq!(compute_adjusted_score(200, &bases, &m), 0.0);
    }

    /// Mixed counts, ambiguous bases ignored (FREQS gives them 0).
    #[test]
    fn test_mixed_composition_exact() {
        let m = make_matrix();
        let mut bases = Vec::new();
        bases.extend(std::iter::repeat_n(0u8, 30)); // A
        bases.extend(std::iter::repeat_n(1u8, 10)); // C
        bases.extend(std::iter::repeat_n(3u8, 20)); // T
        bases.extend(std::iter::repeat_n(14u8, 5)); // N: not counted
        let (a, c, t) = (30.0f64, 10.0f64, 20.0f64);
        let n = a + c + t;
        let t_factor = a * a.ln() + c * c.ln() + t * t.ln() - n * n.ln();
        let t_sum = a * 0.255f64.ln() + c * 0.245f64.ln() + t * 0.255f64.ln() - t_factor;
        let want = 500.0 + t_sum / m.lambda + 0.999;
        let adj = compute_adjusted_score(500, &bases, &m);
        assert!((adj - want).abs() < 1e-9, "adj={adj} want={want}");
    }

    /// Only bases without a frequency: NCBI computes 0 ln 0 = NaN, and the
    /// hit is dropped.
    #[test]
    fn test_no_counted_bases_is_nan_and_dropped() {
        let m = make_matrix();
        assert!(compute_adjusted_score(500, &[14u8; 8], &m).is_nan());
    }
}
