//! Nucleotide substitution matrix in BLASTNA encoding.
//!
//! Reads RMBlast-format custom nucleotide matrices, as NCBI's
//! `BlastScoreBlkNucleotideMatrixRead` (blast_stat.c) does:
//!   - A FREQS line (required, as in NCBI): a `#` comment containing
//!     "FREQS" followed by A, C, G and T frequencies, e.g.
//!     "# FREQS A 0.255 C 0.245 G 0.245 T 0.255"
//!   - A header row of column symbols
//!   - One row per column, in header order: an optional row label, then one
//!     integer score per column
//!
//! Rows are assigned by position, as in NCBI.  Header symbols map through
//! NCBI's `IUPACNA_TO_BLASTNA`: symbols outside the IUPAC codes (X, Z, U, -)
//! land in the gap slot, which is then overwritten, so their scores are
//! discarded exactly as NCBI discards them.
//!
//! Where NCBI accepts a malformed file, this parser refuses it with the line
//! number: a value that is not an integer (NCBI's sscanf reads `9.5` as 9), a
//! row label that is not the header symbol for that row (NCBI ignores labels),
//! a repeated header symbol, FREQS for anything but exactly A, C, G and T
//! once each, or frequencies that do not sum to 1 (NCBI only asserts that in
//! debug builds).  Like NCBI it refuses a missing FREQS line, short or long
//! rows and a non-square matrix.

use std::io::{self, BufRead};
use thiserror::Error;

use crate::encoding::BLASTNA_SIZE;

#[derive(Debug, Error)]
pub enum MatrixError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    #[error("matrix parse error: {0}")]
    Parse(String),
}

/// Gapped Karlin-Altschul parameters declared in a `# KARLIN` comment line of
/// a matrix file (an rmblastn-rs extension; NCBI ignores `#` comments):
///
/// ```text
/// # KARLIN lambda 0.1276 k 0.0179 h 0.0362 alpha 3.52 beta -26.8
/// ```
///
/// Keys are case-insensitive key/value pairs in any order (same style as the
/// `# FREQS` line).  `lambda` and `k` are required for the entry to be usable;
/// exactly one of `h` or `alpha` suffices (alpha = lambda/H); `beta` defaults
/// to 0.  Unset fields are stored as 0.0.  These values feed E-value/bit-score
/// *reporting only* — never the seeding cutoff (see `ka_stats`).
#[derive(Debug, Clone, Copy, Default)]
pub struct KarlinComment {
    pub lambda: f64,
    pub k: f64,
    pub h: f64,
    pub alpha: f64,
    pub beta: f64,
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
    /// Karlin-Altschul parameters from a `# KARLIN` comment line, if present.
    pub karlin: Option<KarlinComment>,
}

impl ScoreMatrix {
    /// Parse a matrix from any `BufRead` source.
    pub fn from_reader<R: BufRead>(name: &str, reader: R) -> Result<Self, MatrixError> {
        let err = |line: usize, msg: String| MatrixError::Parse(format!("line {line}: {msg}"));
        let mut freqs = [0f64; BLASTNA_SIZE];
        let mut freqs_line: Option<usize> = None;
        let mut karlin: Option<KarlinComment> = None;
        let mut header: Option<Vec<u8>> = None; // column symbols, upper-cased
        let mut rows: Vec<Vec<i32>> = Vec::new();

        for (i, line) in reader.lines().enumerate() {
            let lineno = i + 1;
            let line = line?;
            let trimmed = line.trim();
            if trimmed.is_empty() {
                continue;
            }

            if trimmed.starts_with('#') {
                // NCBI: strstr(line, "FREQS") on any comment line
                if let Some(at) = trimmed.find("FREQS") {
                    if let Some(first) = freqs_line {
                        return Err(err(lineno, format!("a second FREQS line (the first is line {first})")));
                    }
                    parse_freqs(&trimmed[at + 5..], &mut freqs).map_err(|m| err(lineno, m))?;
                    freqs_line = Some(lineno);
                } else if let Some(rest) = trimmed.strip_prefix("# KARLIN")
                    .or_else(|| trimmed.strip_prefix("#KARLIN"))
                {
                    karlin = Some(parse_karlin(rest)?);
                }
                continue;
            }

            let tokens: Vec<&str> = trimmed.split_whitespace().collect();
            let Some(cols) = &header else {
                // The first non-comment line is the column header.
                let mut syms: Vec<u8> = Vec::with_capacity(tokens.len());
                for tok in &tokens {
                    let b = tok.as_bytes();
                    if b.len() != 1 || !(b[0].is_ascii_alphabetic() || b[0] == b'-') {
                        return Err(err(lineno, format!(
                            "expected the column header (single-letter symbols), found '{tok}'")));
                    }
                    let sym = b[0].to_ascii_uppercase();
                    if syms.contains(&sym) {
                        return Err(err(lineno, format!("column symbol '{}' appears twice", sym as char)));
                    }
                    syms.push(sym);
                }
                header = Some(syms);
                continue;
            };

            let row = rows.len();
            if row == cols.len() {
                return Err(err(lineno, format!(
                    "more rows than the {} columns in the header", cols.len())));
            }
            // An optional row label, which must be the header symbol for this row.
            let values = match tokens[0].parse::<i32>() {
                Ok(_) => &tokens[..],
                Err(_) => {
                    let want = cols[row] as char;
                    if !tokens[0].eq_ignore_ascii_case(&want.to_string()) {
                        return Err(err(lineno, format!(
                            "row {} is labelled '{}', but column {} of the header is '{want}'",
                            row + 1, tokens[0], row + 1)));
                    }
                    &tokens[1..]
                }
            };
            let mut scores = Vec::with_capacity(values.len());
            for tok in values {
                scores.push(tok.parse::<i32>()
                    .map_err(|_| err(lineno, format!("score '{tok}' is not an integer")))?);
            }
            if scores.len() != cols.len() {
                return Err(err(lineno, format!(
                    "{} scores for {} columns", scores.len(), cols.len())));
            }
            rows.push(scores);
        }

        let Some(cols) = header else {
            return Err(MatrixError::Parse("no column header found".into()));
        };
        if rows.len() != cols.len() {
            return Err(MatrixError::Parse(format!(
                "{} rows for {} columns: the matrix must be square", rows.len(), cols.len())));
        }
        if freqs_line.is_none() {
            return Err(MatrixError::Parse(
                "no FREQS line (\"# FREQS A <f> C <f> G <f> T <f>\"); NCBI requires one".into()));
        }

        // NCBI BlastScoreBlkNucleotideMatrixRead (blast_stat.c) initializes every
        // entry to BLAST_SCORE_MIN (= INT2_MIN = -32768), then overwrites the
        // entries present in the file.  Undefined codes (B/D/H/V) keep -32768.
        // Matching this value exactly matters: i32::MIN/2 (~-1.07e9) changes the
        // gapped-DP arithmetic at ambiguous positions and can diverge from NCBI.
        const BLAST_SCORE_MIN: i32 = -32768; // INT2_MIN
        let mut scores = [[BLAST_SCORE_MIN; BLASTNA_SIZE]; BLASTNA_SIZE];
        // Rows by position, in file order (a later write to a slot wins, as in NCBI)
        for (ri, row_scores) in rows.iter().enumerate() {
            let r = matrix_symbol_to_blastna(cols[ri]) as usize;
            for (ci, &v) in row_scores.iter().enumerate() {
                scores[r][matrix_symbol_to_blastna(cols[ci]) as usize] = v;
            }
        }

        // The GAP code (index 15) is a sentinel between strands.  NCBI sets its
        // whole row and column to INT4_MIN/2 (after the BLAST_SCORE_MIN init).
        for i in 0..BLASTNA_SIZE {
            scores[BLASTNA_SIZE - 1][i] = i32::MIN / 2;
            scores[i][BLASTNA_SIZE - 1] = i32::MIN / 2;
        }

        // Local alignment needs a negative expected score per aligned pair.
        // With a non-negative one there is no Karlin-Altschul lambda, and NCBI
        // 2.17.1 silently reports no hits (RepeatModeler's at.matrix.4.2).
        let expected: f64 = (0..4).flat_map(|i| (0..4).map(move |j| (i, j)))
            .map(|(i, j)| freqs[i] * freqs[j] * scores[i][j] as f64)
            .sum();
        if expected >= 0.0 {
            return Err(MatrixError::Parse(format!(
                "the expected score of an aligned pair is {expected:.3} (from the A/C/G/T scores \
                 and FREQS); it must be negative for local alignment")));
        }

        let lambda = estimate_lambda(&scores, &freqs);

        Ok(ScoreMatrix {
            scores,
            freqs,
            lambda,
            name: name.to_owned(),
            karlin,
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

/// Map a matrix header symbol (already upper-cased) to its BLASTNA index:
/// NCBI's `IUPACNA_TO_BLASTNA` (blast_encoding.c), which the custom matrix
/// reader uses.  The IUPAC codes get their own slots; every other symbol --
/// X, Z, U, '-' -- goes to 15, the gap slot, which `from_reader` then
/// overwrites with the sentinel, so its scores are discarded as in NCBI.
///
/// This differs from sequence encoding, where X is N.  Matrices like
/// 20p43g/18p43g carry both an N column (-1) and an X column (-30), and
/// ctools.matrix a Z column (+3): sharing N's slot would clobber N (bug #33,
/// and for Z the same bug until 3.0.13).
#[inline]
fn matrix_symbol_to_blastna(sym: u8) -> u8 {
    match sym {
        b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3,
        b'R' => 4, b'Y' => 5, b'M' => 6, b'K' => 7, b'W' => 8, b'S' => 9,
        b'B' => 10, b'D' => 11, b'H' => 12, b'V' => 13, b'N' => 14,
        _ => 15,
    }
}

fn parse_karlin(rest: &str) -> Result<KarlinComment, MatrixError> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    let mut kc = KarlinComment::default();
    let mut i = 0;
    while i + 1 < tokens.len() {
        let val: f64 = tokens[i + 1].parse().map_err(|_| {
            MatrixError::Parse(format!("bad KARLIN value '{}' for key '{}'", tokens[i + 1], tokens[i]))
        })?;
        match tokens[i].to_ascii_lowercase().as_str() {
            "lambda" => kc.lambda = val,
            "k"      => kc.k = val,
            "h"      => kc.h = val,
            "alpha"  => kc.alpha = val,
            "beta"   => kc.beta = val,
            other    => return Err(MatrixError::Parse(format!("unknown KARLIN key '{}'", other))),
        }
        i += 2;
    }
    Ok(kc)
}

/// The FREQS pairs after the word "FREQS": exactly A, C, G and T, once each
/// (any order, any case), with frequencies that sum to 1 within NCBI's
/// debug-build tolerance (0.999..1.001).  NCBI requires 4 pairs.
fn parse_freqs(rest: &str, freqs: &mut [f64; BLASTNA_SIZE]) -> Result<(), String> {
    let tokens: Vec<&str> = rest.split_whitespace().collect();
    if tokens.len() != 8 {
        return Err(format!(
            "FREQS needs 4 base/frequency pairs (A, C, G and T), found {} tokens", tokens.len()));
    }
    let mut seen = [false; 4];
    for pair in tokens.chunks(2) {
        let idx = match pair[0].to_ascii_uppercase().as_str() {
            "A" => 0, "C" => 1, "G" => 2, "T" => 3,
            other => return Err(format!("FREQS base '{other}' is not A, C, G or T")),
        };
        if seen[idx] {
            return Err(format!("FREQS gives base '{}' twice", pair[0]));
        }
        seen[idx] = true;
        let f: f64 = pair[1].parse().ok().filter(|f: &f64| f.is_finite() && *f >= 0.0)
            .ok_or_else(|| format!("FREQS value '{}' for {} is not a frequency", pair[1], pair[0]))?;
        freqs[idx] = f;
    }
    let sum: f64 = freqs[..4].iter().sum();
    if !(0.999..=1.001).contains(&sum) {
        return Err(format!("FREQS frequencies sum to {sum}, not 1"));
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
    fn test_karlin_comment_parsed() {
        let with_karlin = format!(
            "# KARLIN lambda 0.1276 K 0.0179 H 0.3562 beta -26.84\n{}",
            SAMPLE_MATRIX
        );
        let mat = ScoreMatrix::from_reader("test", Cursor::new(with_karlin)).unwrap();
        let kc = mat.karlin.expect("KARLIN line should be parsed");
        assert!((kc.lambda - 0.1276).abs() < 1e-12);
        assert!((kc.k - 0.0179).abs() < 1e-12);
        assert!((kc.h - 0.3562).abs() < 1e-12);
        assert_eq!(kc.alpha, 0.0);
        assert!((kc.beta - -26.84).abs() < 1e-12);
        // Matrix body still parsed normally alongside the comment.
        assert_eq!(mat.score(0, 0), 9);

        // No KARLIN line -> None; unknown key -> parse error.
        assert!(ScoreMatrix::from_reader("test", Cursor::new(SAMPLE_MATRIX)).unwrap().karlin.is_none());
        let bad = format!("# KARLIN lambdo 0.1\n{}", SAMPLE_MATRIX);
        assert!(ScoreMatrix::from_reader("test", Cursor::new(bad)).is_err());
    }

    fn parse(text: &str) -> Result<ScoreMatrix, MatrixError> {
        ScoreMatrix::from_reader("test", Cursor::new(text))
    }

    /// `SAMPLE_MATRIX` with `from` replaced by `to` (once).
    fn sample_with(from: &str, to: &str) -> String {
        assert!(SAMPLE_MATRIX.contains(from), "{from:?}");
        SAMPLE_MATRIX.replacen(from, to, 1)
    }

    /// Each malformed matrix is refused with a message naming the problem.
    /// 3.0.12 accepted all but one of these and mis-scored silently.
    #[test]
    fn malformed_matrices_are_rejected() {
        let cases: Vec<(&str, String, &str)> = vec![
            // NCBI's sscanf reads 9.5 as 9; 3.0.12 dropped the whole row
            ("non-integer score", sample_with("A  9  -6", "A  9.5  -6"), "line 4: score '9.5' is not an integer"),
            ("short row", sample_with("C -6   9 -13  -5  -1", "C -6   9 -13  -5"), "line 5: 4 scores for 5 columns"),
            ("long row", sample_with("C -6   9 -13  -5  -1", "C -6   9 -13  -5  -1 7"), "6 scores for 5 columns"),
            ("missing row", sample_with("N -1  -1  -1  -1  -1\n", ""), "4 rows for 5 columns"),
            ("extra row", format!("{SAMPLE_MATRIX}N -1 -1 -1 -1 -1\n"), "more rows than the 5 columns"),
            // NCBI ignores labels and assigns rows by position
            ("label out of order", sample_with("G -5", "T -5"), "row 3 is labelled 'T', but column 3 of the header is 'G'"),
            ("repeated header symbol", sample_with("   A   C   G   T   N", "   A   C   G   T   A"), "'A' appears twice"),
            ("no FREQS", sample_with("# FREQS A 0.255 C 0.245 G 0.245 T 0.255\n", ""), "no FREQS line"),
            ("FREQS missing T", sample_with(" T 0.255", ""), "4 base/frequency pairs"),
            ("FREQS base twice", sample_with("G 0.245", "A 0.245"), "base 'A' twice"),
            ("FREQS not a number", sample_with("C 0.245", "C x"), "'x' for C is not a frequency"),
            ("FREQS sum", sample_with("A 0.255", "A 0.355"), "sum to"),
            ("two FREQS lines", format!("# FREQS A 0.25 C 0.25 G 0.25 T 0.25\n{SAMPLE_MATRIX}"), "a second FREQS line"),
            ("data before header", "# FREQS A 0.25 C 0.25 G 0.25 T 0.25\n1 2 3\n".to_string(), "expected the column header"),
            ("positive expected score", sample_with("A  9  -6  -5 -13", "A  9  60  50 -13"), "must be negative"),
        ];
        for (label, text, needle) in cases {
            match parse(&text) {
                Err(e) => assert!(e.to_string().contains(needle), "{label}: {e}"),
                Ok(_) => panic!("{label}: accepted"),
            }
        }
    }

    /// NCBI finds FREQS with strstr on any comment line, so `#  FREQS ...`
    /// counts (3.0.12 missed it, leaving lambda 0 and -complexity_adjust off).
    #[test]
    fn freqs_found_anywhere_on_a_comment_line() {
        let m = parse(&sample_with("# FREQS", "#  my matrix FREQS")).unwrap();
        assert!((m.freqs[0] - 0.255).abs() < 1e-12);
        assert!(m.lambda > 0.0);
        assert_eq!(m.lambda, parse(SAMPLE_MATRIX).unwrap().lambda);
    }

    /// Unlabelled rows (assigned by position) and lower-case symbols load the
    /// same matrix as the labelled upper-case form.
    #[test]
    fn unlabelled_rows_and_lower_case_symbols() {
        let unlabelled: String = SAMPLE_MATRIX.lines()
            .map(|l| match l.chars().next() { Some('A'|'C'|'G'|'T'|'N') => &l[1..], _ => l })
            .collect::<Vec<_>>().join("\n");
        let lower = sample_with("   A   C   G   T   N", "   a   c   g   t   n").replacen("\nG -5", "\ng -5", 1);
        let want = parse(SAMPLE_MATRIX).unwrap().scores;
        assert_eq!(parse(&unlabelled).unwrap().scores, want);
        assert_eq!(parse(&lower).unwrap().scores, want);
    }

    /// Non-IUPAC header symbols (Z here, as in ctools.matrix) go to the gap
    /// slot and are discarded, as in NCBI; 3.0.12 mapped Z to N, so the Z
    /// column's +3 replaced N's -1.
    #[test]
    fn z_column_does_not_clobber_n() {
        let text = "# FREQS A 0.25 C 0.25 G 0.25 T 0.25
   A   C   G   T   N   Z
A  9  -6  -5 -13  -1   3
C -6   9 -13  -5  -1   3
G -5 -13   9  -6  -1   3
T -13  -5  -6   9  -1   3
N -1  -1  -1  -1  -1   3
Z  3   3   3   3   3   3
";
        let m = parse(text).unwrap();
        assert_eq!(m.score(0, 14), -1);
        assert_eq!(m.score(14, 0), -1);
        assert_eq!(m.score(14, 14), -1);
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
