//! Ungapped extension from a seed hit.
//!
//! Ported from NCBI's s_NuclUngappedExtendExact (na_ungapped.c).
//!
//! NCBI operates in (RC_query, FWD_subject) for both strands.  Its left extension
//! uses fixed xdrop; its right extension uses adaptive xdrop (X_current tightens
//! once the combined score exceeds x_dropoff, preventing the running total from
//! going negative).
//!
//! Rust uses (FWD_query, RC_subject) internally.  For plus-strand hits the directions
//! align directly with NCBI's (Rust left = NCBI left, Rust right = NCBI right).
//! For minus-strand hits the conventions are anti-parallel: Rust's RIGHT extension
//! physically traverses the same bases as NCBI's LEFT extension, and vice versa.
//! The `plus_strand` flag controls which Rust side gets fixed vs adaptive xdrop so
//! that the algorithm always matches NCBI regardless of the coordinate frame used.

use crate::encoding::BLASTNA_COMPLEMENT;
use crate::matrix::ScoreMatrix;

/// Result of an ungapped extension attempt.
#[derive(Debug, Clone)]
pub struct UngappedResult {
    pub score: i32,
    /// 0-based offset of the first aligned query base.
    pub q_start: u32,
    /// 0-based exclusive end offset in the query.
    pub q_end: u32,
    /// 0-based offset of the first aligned subject base.
    pub s_start: u32,
    /// 0-based exclusive end offset in the subject.
    pub s_end: u32,
}

/// Extend ungapped alignment from seed (q_off, s_off).
///
/// `query` and `subject` are BLASTNA-encoded with a leading sentinel at index 0.
/// `q_off` and `s_off` are 0-based offsets into the real bases.
/// `plus_strand`: true for plus-strand hits (FWD_q/FWD_s), false for minus-strand
///   hits in (FWD_q/RC_s) space — controls which direction gets adaptive xdrop.
///
/// Returns `None` if the total ungapped score is below `min_score`.
pub fn extend_ungapped(
    query: &[u8],
    subject: &[u8],
    q_off: u32,
    s_off: u32,
    matrix: &ScoreMatrix,
    x_dropoff: i32,
    min_score: i32,
    plus_strand: bool,
) -> Option<UngappedResult> {
    // X is negative: threshold for breaking (break when sum < X).
    let x = -x_dropoff;

    let q_real_len = (query.len() - 1) as i64;
    let s_real_len = (subject.len() - 1) as i64;
    let q_off_i = q_off as i64;
    let s_off_i = s_off as i64;

    if plus_strand {
        // ── Plus strand: NCBI LEFT = Rust LEFT (fixed xdrop first) ──────────────
        //   Left extension (fixed xdrop) runs first; score carries into right.
        //   Right extension uses adaptive xdrop.

        let mut score = 0i32;
        let mut sum = 0i32;
        let mut q_beg = q_off_i;

        {
            let mut q = q_off_i - 1;
            let mut s = s_off_i - 1;
            loop {
                if q < 0 || s < 0 { break; }
                let qb = query[1 + q as usize];
                let sb = subject[1 + s as usize];
                if qb == 15 || sb == 15 { break; }
                sum += matrix.score(qb, sb);
                if sum > 0 {
                    q_beg = q;
                    score += sum;
                    sum = 0;
                } else if sum < x {
                    break;
                }
                q -= 1;
                s -= 1;
            }
        }

        let q_start = q_beg as u32;
        let s_start = (s_off_i - (q_off_i - q_beg)) as u32;

        // Right extension (adaptive xdrop).
        let mut q_end = q_off_i;
        let mut x_current = x;
        sum = 0;

        {
            let mut q = q_off_i;
            let mut s = s_off_i;
            loop {
                if q >= q_real_len || s >= s_real_len { break; }
                let qb = query[1 + q as usize];
                let sb = subject[1 + s as usize];
                if qb == 15 || sb == 15 { break; }
                sum += matrix.score(qb, sb);
                if sum > 0 {
                    q_end = q + 1;
                    score += sum;
                    x_current = if -score > x { -score } else { x };
                    sum = 0;
                } else if sum < x_current {
                    break;
                }
                q += 1;
                s += 1;
            }
        }

        let s_end = (s_off_i + (q_end - q_off_i)) as u32;

        if score < min_score {
            None
        } else {
            Some(UngappedResult { score, q_start, q_end: q_end as u32, s_start, s_end })
        }

    } else {
        // ── Minus strand: NCBI LEFT = Rust RIGHT (fixed xdrop first) ────────────
        //   In (FWD_q, RC_s) space, Rust's RIGHT extension traverses the same
        //   physical bases as NCBI's LEFT extension, so it gets fixed xdrop.
        //   Rust's LEFT extension = NCBI's RIGHT extension → adaptive xdrop.

        // Step 1: Right extension with fixed xdrop (= NCBI's left extension).
        let mut score = 0i32;
        let mut sum = 0i32;
        let mut q_end = q_off_i;

        {
            let mut q = q_off_i;
            let mut s = s_off_i;
            loop {
                if q >= q_real_len || s >= s_real_len { break; }
                let qb = query[1 + q as usize];
                let sb = subject[1 + s as usize];
                if qb == 15 || sb == 15 { break; }
                // Score in NCBI's (RC_query, FWD_subject) frame: this loop walks the
                // alignment in (FWD_query, RC_subject) space, so complement both bases
                // before the matrix lookup.  No-op for the complement-symmetric ACGT
                // core; required because the matrix is NOT complement-symmetric for
                // ambiguity codes (e.g. S/A=-10 vs S/T=-11), which otherwise shifts the
                // committed endpoint in low-complexity regions and diverges from NCBI.
                sum += matrix.score(BLASTNA_COMPLEMENT[qb as usize], BLASTNA_COMPLEMENT[sb as usize]);
                if sum > 0 {
                    q_end = q + 1;
                    score += sum;
                    sum = 0;
                } else if sum < x {
                    break;
                }
                q += 1;
                s += 1;
            }
        }

        // Step 2: Left extension with adaptive xdrop (= NCBI's right extension).
        let mut q_beg = q_off_i;
        let mut x_current = x;
        sum = 0;

        {
            let mut q = q_off_i - 1;
            let mut s = s_off_i - 1;
            loop {
                if q < 0 || s < 0 { break; }
                let qb = query[1 + q as usize];
                let sb = subject[1 + s as usize];
                if qb == 15 || sb == 15 { break; }
                // See Step 1: complement both bases to score in NCBI's frame.
                sum += matrix.score(BLASTNA_COMPLEMENT[qb as usize], BLASTNA_COMPLEMENT[sb as usize]);
                if sum > 0 {
                    q_beg = q;
                    score += sum;
                    x_current = if -score > x { -score } else { x };
                    sum = 0;
                } else if sum < x_current {
                    break;
                }
                q -= 1;
                s -= 1;
            }
        }

        let q_start = q_beg as u32;
        let s_start = (s_off_i - (q_off_i - q_beg)) as u32;
        let s_end = (s_off_i + (q_end - q_off_i)) as u32;

        if score < min_score {
            None
        } else {
            Some(UngappedResult { score, q_start, q_end: q_end as u32, s_start, s_end })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encoding::encode_iupac;
    use crate::matrix::ScoreMatrix;
    use std::io::Cursor;

    const MAT: &str = r"
# FREQS A 0.25 C 0.25 G 0.25 T 0.25
   A   C   G   T
A  5  -4  -4  -4
C -4   5  -4  -4
G -4  -4   5  -4
T -4  -4  -4   5
";

    fn make_matrix() -> ScoreMatrix {
        ScoreMatrix::from_reader("test", Cursor::new(MAT)).unwrap()
    }

    #[test]
    fn test_extend_perfect() {
        let q = encode_iupac(b"ACGTACGT");
        let s = encode_iupac(b"ACGTACGT");
        let m = make_matrix();
        let r = extend_ungapped(&q, &s, 0, 0, &m, 20, 0, true).unwrap();
        assert_eq!(r.score, 5 * 8);
        assert_eq!(r.q_start, 0);
        assert_eq!(r.q_end, 8);
    }

    #[test]
    fn test_extend_too_low_score() {
        let q = encode_iupac(b"AAAAAAAA");
        let s = encode_iupac(b"CCCCCCCC");
        let m = make_matrix();
        let r = extend_ungapped(&q, &s, 0, 0, &m, 20, 1, true);
        assert!(r.is_none(), "all-mismatch should fail min_score=1 cutoff");
    }

    #[test]
    fn test_left_extension() {
        let q = encode_iupac(b"ACGTACGT");
        let s = encode_iupac(b"ACGTACGT");
        let m = make_matrix();
        let r = extend_ungapped(&q, &s, 4, 4, &m, 20, 0, true).unwrap();
        assert_eq!(r.score, 5 * 8);
        assert_eq!(r.q_start, 0);
        assert_eq!(r.q_end, 8);
    }

    #[test]
    fn test_score_carries_from_left_to_right() {
        let q = encode_iupac(b"AAAAAAAA");
        let s = encode_iupac(b"AAAAAAAA");
        let m = make_matrix();
        let r = extend_ungapped(&q, &s, 4, 4, &m, 20, 0, true).unwrap();
        assert_eq!(r.score, 5 * 8);
        assert_eq!(r.q_start, 0);
        assert_eq!(r.q_end, 8);
    }

    // ── Bug #36 regression: minus-strand must extend in NCBI's (RC_q, FWD_s) frame ──
    //
    // NCBI runs minus-strand ungapped extension in (RC_query, FWD_subject) space with the
    // SAME left-then-right plus-strand code: the adaptive RIGHT extension begins AT the
    // seed base, so a matching seed/word commits positive score and tightens X_current
    // before the extension reaches any low-scoring region.
    //
    // The legacy mirror frame (FWD_query, RC_subject) with plus_strand=false runs the
    // adaptive extension starting just PAST the word, so X_current never tightens from the
    // word.  Across a long, mildly-negative stretch (which never trips the loose -xdrop
    // threshold) the legacy path OVER-extends into an adjacent good region, while NCBI
    // stops short.  This is bug #36: the over-extension set a diagonal-dedup endpoint that
    // suppressed a real seed, dropping a whole HSP.  The production code (collect_ungapped*)
    // therefore extends minus strand via the NCBI frame; this test pins that they differ
    // and that the NCBI frame is the one that stops early.

    /// Production minus-strand extension: map the (FWD_q, RC_s) anchor to (RC_q, FWD_s)
    /// and run the plus-strand code, exactly as collect_ungapped_combined does.
    fn minus_via_ncbi_frame(
        q_fwd: &[u8], s_rc: &[u8], q_off: u32, s_off: u32,
        m: &ScoreMatrix, xdrop: i32,
    ) -> Option<UngappedResult> {
        let q_len = (q_fwd.len() - 2) as u32;
        let s_len = (s_rc.len() - 2) as u32;
        let q_rc  = crate::encoding::revcomp_blastna(q_fwd);   // sentinel-wrapped RC query
        let s_fwd = crate::encoding::revcomp_blastna(s_rc);    // sentinel-wrapped FWD subject
        // Same anchor mapping as production: qn = q_len-1-q_off, sn = s_len-1-s_off.
        let qn = q_len - 1 - q_off;
        let sn = s_len - 1 - s_off;
        let r = extend_ungapped(&q_rc, &s_fwd, qn, sn, m, xdrop, i32::MIN, true)?;
        Some(UngappedResult {
            score: r.score,
            q_start: q_len - r.q_end,
            q_end:   q_len - r.q_start,
            s_start: s_len - r.s_end,
            s_end:   s_len - r.s_start,
        })
    }

    #[test]
    fn test_minus_ncbi_frame_stops_early_bug36() {
        // Build the bug-#36 shape in the FWD_q/RC_s frame on the diagonal (q_off==s_off).
        // Reading toward decreasing offset from the anchor (the adaptive physical
        // direction): [seed/word: 6 matches][mild stretch: short negatives][long good run].
        // make_matrix scores match=+5, mismatch=-4.  A large xdrop keeps the legacy loose
        // threshold from ever tripping across the mild stretch.
        let m = make_matrix();
        // Layout (low→high index): a long good run, a deep mismatch stretch that exceeds the
        // NCBI-frame's tightened X_current but stays within the loose -xdrop, then a
        // match-rich word/context block.  match=+5, mismatch=-4, xdrop large (=200).
        //   good run G×24 | mismatch ×20 | word/context ACGTACGTACGT
        let qs = b"GGGGGGGGGGGGGGGGGGGGGGGGACACACACACACACACACACACGTACGTACGT";
        let ss = b"GGGGGGGGGGGGGGGGGGGGGGGGCACACACACACACACACACAACGTACGTACGT";
        let q = encode_iupac(qs);
        let s = encode_iupac(ss);
        let xdrop = 200;

        // Scan interior anchors; the divergence appears for anchors sitting in/just above
        // the word block, where the NCBI frame commits the word early (tightening X_current)
        // and then stops in the mismatch stretch while the legacy frame sails across it.
        let n = qs.len() as u32;
        let mut found_divergence = false;
        for anchor in 1..n - 1 {
            let legacy = extend_ungapped(&q, &s, anchor, anchor, &m, xdrop, i32::MIN, false).unwrap();
            let ncbi   = minus_via_ncbi_frame(&q, &s, anchor, anchor, &m, xdrop).unwrap();
            if legacy.score != ncbi.score {
                // When they differ, the legacy mirror frame is the one that over-extends
                // (higher score because it crossed the mismatch stretch into more matches).
                assert!(legacy.score > ncbi.score,
                    "anchor={anchor}: legacy should over-extend: ncbi={ncbi:?} legacy={legacy:?}");
                found_divergence = true;
            }
        }
        assert!(found_divergence,
            "expected at least one anchor where the legacy (FWD_q,RC_s) frame over-extends \
             relative to NCBI's (RC_q,FWD_s) frame — this is the bug #36 mechanism");
    }

    #[test]
    fn test_minus_ncbi_frame_clean_match_roundtrips() {
        // On a clean all-match region the two frames agree (no adaptive divergence),
        // and the NCBI-frame mapping round-trips coordinates correctly.
        let m = make_matrix();
        let q = encode_iupac(b"ACGTACGTACGTACGT");
        let s = encode_iupac(b"ACGTACGTACGTACGT");
        let legacy = extend_ungapped(&q, &s, 8, 8, &m, 30, i32::MIN, false).unwrap();
        let ncbi   = minus_via_ncbi_frame(&q, &s, 8, 8, &m, 30).unwrap();
        assert_eq!(legacy.score, ncbi.score);
        assert_eq!((legacy.q_start, legacy.q_end), (ncbi.q_start, ncbi.q_end));
        assert_eq!((legacy.s_start, legacy.s_end), (ncbi.s_start, ncbi.s_end));
        assert_eq!(ncbi.score, 5 * 16);
    }
}
