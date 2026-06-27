//! Complexity-adjusted scoring — port of blast_traceback.c lines 522-574.
//!
//! Formula: adj_score = raw_score + t_sum / lambda   (+ 0.999 for rounding)
//! where t_sum corrects for how far the aligned query composition deviates
//! from the background base frequencies in the scoring matrix.

use crate::hits::{EditOp, EditScript};
use crate::matrix::ScoreMatrix;
use crate::encoding::BLASTNA_SIZE;

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
        return raw_score as f64;
    }

    t_factor -= t_counts * t_counts.ln();
    t_sum -= t_factor;

    let adj = raw_score as f64 + t_sum / matrix.lambda + 0.999;
    if adj < 0.0 { 0.0 } else { adj }
}

/// Apply complexity adjustment; return None if score falls below cutoff.
pub fn apply_complexity_adjust(
    raw_score: i32,
    query: &[u8],
    q_start: u32,
    edit_script: &EditScript,
    matrix: &ScoreMatrix,
    cutoff: i32,
) -> Option<i32> {
    let aligned_bases = collect_query_aligned_bases(query, q_start, edit_script);
    let adj = compute_adjusted_score(raw_score, &aligned_bases, matrix);
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

    #[test]
    fn test_balanced_composition_unchanged() {
        let m = make_matrix();
        let bases: Vec<u8> = (0..16).flat_map(|_| vec![0u8, 1, 2, 3]).collect(); // ACGT×4
        let adj = compute_adjusted_score(100, &bases, &m);
        assert!((adj - 100.0).abs() < 5.0, "adj={}", adj);
    }

    #[test]
    fn test_biased_composition_reduced() {
        let m = make_matrix();
        let bases = vec![0u8; 64]; // all A
        let adj = compute_adjusted_score(200, &bases, &m);
        assert!(adj < 200.0, "biased composition should reduce score, got {}", adj);
    }
}
