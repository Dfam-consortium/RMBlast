//! Nucleotide substitution matrix in BLASTNA encoding.
//!
//! Reads RMBlast-format custom nucleotide matrices:
//!   - Optional FREQS line: "# FREQS A 0.255 C 0.245 G 0.245 T 0.255"
//!   - Header row of column symbols (IUPAC)
//!   - Rows: row-symbol followed by integer scores
//!
//! After loading, the 16×16 BLASTNA matrix is filled by remapping each IUPAC
//! symbol to its BLASTNA index.  Lambda is computed (or read from file) for
//! use by the complexity adjustment algorithm.

use std::io::{self, BufRead};
use thiserror::Error;

use crate::encoding::{BLASTNA_SIZE, IUPAC_TO_BLASTNA};

#[derive(Debug, Error)]
pub enum MatrixError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    #[error("matrix parse error: {0}")]
    Parse(String),
    #[error("missing required matrix symbol '{0}'")]
    MissingSymbol(char),
}

/// A 16×16 scoring matrix in BLASTNA encoding, plus optional base frequencies
/// and lambda (used by complexity adjustment).
#[derive(Debug, Clone)]
pub struct ScoreMatrix {
    /// scores[query_blastna][subject_blastna]
    pub scores: [[i32; BLASTNA_SIZE]; BLASTNA_SIZE],
    /// Base frequencies in BLASTNA order (A,C,G,T at indices 0-3; rest 0).
    pub freqs: [f64; BLASTNA_SIZE],
    /// Lambda parameter estimated from the matrix + frequencies.
    pub lambda: f64,
    pub name: String,
}

impl ScoreMatrix {
    /// Parse a matrix from any `BufRead` source.
    pub fn from_reader<R: BufRead>(name: &str, reader: R) -> Result<Self, MatrixError> {
        let mut freqs = [0f64; BLASTNA_SIZE];
        let mut col_order: Vec<u8> = Vec::new(); // BLASTNA indices of columns
        let mut raw_rows: Vec<(u8, Vec<i32>)> = Vec::new(); // (blastna_row, scores)

        for line in reader.lines() {
            let line = line?;
            let trimmed = line.trim();

            if trimmed.is_empty() || trimmed.starts_with('#') {
                // Look for FREQS annotation
                if let Some(rest) = trimmed.strip_prefix("# FREQS") {
                    parse_freqs(rest, &mut freqs)?;
                } else if let Some(rest) = trimmed.strip_prefix("#FREQS") {
                    parse_freqs(rest, &mut freqs)?;
                }
                continue;
            }

            let tokens: Vec<&str> = trimmed.split_whitespace().collect();
            if tokens.is_empty() {
                continue;
            }

            // First non-comment non-empty line without a leading symbol that
            // parses as integer is the column header.
            if col_order.is_empty() && tokens[0].len() == 1 {
                // Heuristic: if the first token is a single letter and the
                // second token is also a single letter (or the line has many
                // single letters), it's the header row.
                let all_letters = tokens.iter().all(|t| t.len() == 1 && t.chars().next().map_or(false, |c| c.is_ascii_alphabetic() || *t == "-"));
                if all_letters {
                    for tok in &tokens {
                        let ch = tok.chars().next().unwrap() as u8;
                        col_order.push(matrix_symbol_to_blastna(ch));
                    }
                    continue;
                }
            }

            if !col_order.is_empty() {
                // Two supported formats:
                //   1. With row label:    "A  9  1 -6 ..."  — first token is an IUPAC letter
                //   2. Without row label: "9  1 -6 ..."     — all tokens are integers, row
                //                          position in the file determines the IUPAC symbol
                //                          (same order as the column header).
                //
                // Format 2 is used by NCBI's FREQS matrices (comparison.matrix, 14p35g.matrix,
                // etc.).  Detect by checking whether the first token parses as an integer.
                let ch = tokens[0].chars().next().unwrap() as u8;
                let has_label = (ch.is_ascii_alphabetic() || ch == b'-') && tokens[0].len() == 1;
                let (row_idx, score_tokens): (u8, &[&str]) = if has_label {
                    let idx = matrix_symbol_to_blastna(ch);
                    (idx, &tokens[1..])
                } else if raw_rows.len() < col_order.len() {
                    // No row label: assign row in col_order sequence order.
                    (col_order[raw_rows.len()], &tokens[..])
                } else {
                    continue; // more rows than columns — ignore trailing garbage
                };
                let scores: Result<Vec<i32>, _> = score_tokens.iter().map(|t| t.parse::<i32>()).collect();
                match scores {
                    Ok(s) => raw_rows.push((row_idx, s)),
                    Err(_) => continue, // skip unparsable rows (e.g. trailing whitespace lines)
                }
            }
        }

        if col_order.is_empty() {
            return Err(MatrixError::Parse("no column header found".into()));
        }
        if raw_rows.is_empty() {
            return Err(MatrixError::Parse("no data rows found".into()));
        }

        // NCBI BlastScoreBlkNucleotideMatrixRead (blast_stat.c) initializes every
        // entry to BLAST_SCORE_MIN (= INT2_MIN = -32768), then overwrites the
        // entries present in the file.  Undefined codes (B/D/H/V) keep -32768.
        // Matching this value exactly matters: i32::MIN/2 (~-1.07e9) changes the
        // gapped-DP arithmetic at ambiguous positions and can diverge from NCBI.
        const BLAST_SCORE_MIN: i32 = -32768; // INT2_MIN
        let mut scores = [[BLAST_SCORE_MIN; BLASTNA_SIZE]; BLASTNA_SIZE];
        // Fill in parsed values
        for (row_idx, row_scores) in &raw_rows {
            for (ci, &col_idx) in col_order.iter().enumerate() {
                if ci < row_scores.len() {
                    scores[*row_idx as usize][col_idx as usize] = row_scores[ci];
                }
            }
        }

        // The GAP code (index 15) is a sentinel between strands.  NCBI sets its
        // whole row and column to INT4_MIN/2 (after the BLAST_SCORE_MIN init).
        for i in 0..BLASTNA_SIZE {
            scores[BLASTNA_SIZE - 1][i] = i32::MIN / 2;
            scores[i][BLASTNA_SIZE - 1] = i32::MIN / 2;
        }

        let lambda = estimate_lambda(&scores, &freqs);

        Ok(ScoreMatrix {
            scores,
            freqs,
            lambda,
            name: name.to_owned(),
        })
    }

    /// Parse a matrix from a file path.
    pub fn from_file(path: &str) -> Result<Self, MatrixError> {
        let file = std::fs::File::open(path).map_err(MatrixError::Io)?;
        let reader = io::BufReader::new(file);
        Self::from_reader(path, reader)
    }

    #[inline]
    pub fn score(&self, q: u8, s: u8) -> i32 {
        self.scores[(q & 15) as usize][(s & 15) as usize]
    }
}

/// Map a matrix column/row symbol to its BLASTNA index, matching NCBI's
/// `IUPACNA_TO_BLASTNA` table (blast_encoding.c) as used by the custom matrix
/// reader (`BlastScoreBlkNucleotideMatrixRead`).  This DIFFERS from sequence
/// encoding for two symbols:
///   - '-'  -> 15 (gap slot)
///   - 'X'  -> 15 (gap slot)   [sequence encoding sends 'X' -> 14 (N)]
/// Sending 'X' to 15 here is essential: matrices like 20p43g/18p43g carry BOTH
/// an 'N' column (-1) and an 'X' column (-30).  If 'X' shared N's index 14 the
/// X row/column would clobber N (bug #33).  At 15 the X entries land in the gap
/// row/column, which `from_reader` overwrites with the sentinel -- discarded,
/// exactly as NCBI does.
#[inline]
fn matrix_symbol_to_blastna(ch: u8) -> u8 {
    match ch {
        b'-' | b'X' | b'x' => 15,
        _ => IUPAC_TO_BLASTNA[ch as usize],
    }
}

fn parse_freqs(rest: &str, freqs: &mut [f64; BLASTNA_SIZE]) -> Result<(), MatrixError> {
    use crate::encoding::IUPAC_TO_BLASTNA;
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    let mut i = 0;
    while i + 1 < tokens.len() {
        let sym_bytes = tokens[i].as_bytes();
        if sym_bytes.len() == 1 {
            let idx = IUPAC_TO_BLASTNA[sym_bytes[0] as usize] as usize;
            let f: f64 = tokens[i + 1].parse().map_err(|_| {
                MatrixError::Parse(format!("bad frequency value '{}'", tokens[i + 1]))
            })?;
            freqs[idx] = f;
        }
        i += 2;
    }
    Ok(())
}

/// Estimate lambda matching NCBI's blast_stat.c algorithm exactly.
/// Uses frequency-weighted score sum: Σ_{i,j} f_i * f_j * exp(λ * s_{ij}) = 1
/// Only uses the unambiguous bases (indices 0-3, i.e. A,C,G,T).
/// Falls back to 0.0 if frequencies are not set.
///
/// Algorithm: start at lambda=0.5, double until sum>=1 to find upper bound,
/// then bisect with tolerance 0.00001 (matching NCBI's BLAST_KARLIN_LAMBDA_ACCURACY_DEFAULT).
fn estimate_lambda(scores: &[[i32; BLASTNA_SIZE]; BLASTNA_SIZE], freqs: &[f64; BLASTNA_SIZE]) -> f64 {
    let total_freq: f64 = freqs[0..4].iter().sum();
    if total_freq < 1e-9 {
        return 0.0;
    }

    let eval = |lambda: f64| -> f64 {
        let mut sum = 0f64;
        for i in 0..BLASTNA_SIZE {
            for j in 0..BLASTNA_SIZE {
                if freqs[i] > 0.0 && freqs[j] > 0.0 {
                    sum += freqs[i] * freqs[j] * (lambda * scores[i][j] as f64).exp();
                }
            }
        }
        sum
    };

    // Phase 1: start at 0.5, double until sum >= 1.0 to find upper bound
    let mut lambda_lower = 0.0f64;
    let mut lambda = 0.5f64;
    loop {
        let sum = eval(lambda);
        if sum < 1.0 {
            lambda_lower = lambda;
            lambda *= 2.0;
            if lambda > 1000.0 {
                return 0.0; // no solution
            }
        } else {
            break;
        }
    }
    let mut lambda_upper = lambda;

    // Phase 2: bisect until interval <= 0.00001 (NCBI's tolerance)
    while lambda_upper - lambda_lower > 0.00001 {
        lambda = (lambda_lower + lambda_upper) / 2.0;
        let sum = eval(lambda);
        if sum >= 1.0 {
            lambda_upper = lambda;
        } else {
            lambda_lower = lambda;
        }
    }

    // Return the last midpoint computed (same as NCBI's sbp->matrix->lambda = lambda)
    lambda
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const SAMPLE_MATRIX: &str = r"
# FREQS A 0.255 C 0.245 G 0.245 T 0.255
   A   C   G   T   N
A  9  -6  -5 -13  -1
C -6   9 -13  -5  -1
G -5 -13   9  -6  -1
T -13  -5  -6   9  -1
N -1  -1  -1  -1  -1
";

    #[test]
    fn test_parse_sample() {
        let mat = ScoreMatrix::from_reader("test", Cursor::new(SAMPLE_MATRIX)).unwrap();
        assert_eq!(mat.score(0, 0), 9);  // A vs A
        assert_eq!(mat.score(0, 1), -6); // A vs C
        assert!(mat.lambda > 0.0);
        assert!((mat.freqs[0] - 0.255).abs() < 1e-9); // A freq
    }

    // Matrix carrying BOTH an N column/row (-1) and an X column/row (-30), like
    // the RepeatMasker species matrices (20p43g/18p43g).  Regression for bug #33:
    // 'X' must map to the gap slot (15), NOT to N's index (14), so it cannot
    // clobber the N row/column.
    const NX_MATRIX: &str = r"
# FREQS A 0.255 C 0.245 G 0.245 T 0.255
   A   C   G   T   N   X
A  9  -6  -5 -13  -1 -30
C -6   9 -13  -5  -1 -30
G -5 -13   9  -6  -1 -30
T -13  -5  -6   9  -1 -30
N -1  -1  -1  -1  -1 -30
X -30 -30 -30 -30 -30 -30
";

    #[test]
    fn test_x_column_does_not_clobber_n() {
        let mat = ScoreMatrix::from_reader("test", Cursor::new(NX_MATRIX)).unwrap();
        // Subject N (BLASTNA 14) must keep the N column's -1, not the X column's -30.
        assert_eq!(mat.score(0, 14), -1, "A vs N should be -1 (N col), not -30 (X col)");
        assert_eq!(mat.score(3, 14), -1, "T vs N should be -1");
        // Query N row preserved too.
        assert_eq!(mat.score(14, 0), -1, "N vs A should be -1");
        // Core ACGT scores unaffected by the trailing X column.
        assert_eq!(mat.score(0, 0), 9);
        assert_eq!(mat.score(3, 0), -13);
    }

    #[test]
    fn test_sequence_x_encodes_as_n() {
        // In a *sequence*, 'X' (RepeatMasker hard-mask char) must encode as N (14),
        // not the gap slot (15) — otherwise it renders as '-' and drops from the
        // alignment.  This is the inverse of the matrix-column mapping above.
        use crate::encoding::{IUPAC_TO_BLASTNA, encode_iupac};
        assert_eq!(IUPAC_TO_BLASTNA[b'X' as usize], 14);
        assert_eq!(IUPAC_TO_BLASTNA[b'x' as usize], 14);
        let enc = encode_iupac(b"ACXTN"); // [sentinel, A, C, X->N, T, N, sentinel]
        let inner = &enc[1..enc.len() - 1];
        assert_eq!(inner, &[0u8, 1, 14, 3, 14]);
    }
}
