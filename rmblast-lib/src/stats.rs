//! Post-alignment statistics: Kimura divergence, CpG-adjusted divergence,
//! transition/transversion counts, and percent gap computation.
//!
//! Ported from objtools/align_format/showalign.cpp lines 2143-2361.

use crate::encoding::BLASTNA_TO_IUPAC;

/// Classification of an aligned base pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PairType {
    Identity,
    Transition,
    Transversion,
    /// IUPAC ambiguity pair where the two bases share at least one base in common — not a substitution.
    IubMatch,
    /// IUPAC ambiguity pair where the two bases have NO bases in common — counts as a substitution
    /// for perc_sub but not as a transition/transversion (following NCBI's showalign behavior).
    IubSub,
    Gap,
    Unknown,
}


/// Statistics computed from one aligned pair of sequences.
#[derive(Debug, Clone, Default)]
pub struct AlignStats {
    pub matches: u32,
    pub mismatches: u32,
    pub query_gaps: u32,    // inserted in subject
    pub subject_gaps: u32,  // deleted from subject
    pub transitions: u32,
    pub transversions: u32,
    pub cpg_sites: u32,
    /// Kimura divergence (× 100, i.e. percentage).
    pub kdiv: f64,
    /// CpG-adjusted Kimura divergence (× 100).
    pub cpg_kdiv: f64,
    /// % substitutions relative to query length.
    pub perc_sub: f32,
    /// % query gaps relative to query length.
    pub perc_query_gap: f32,
    /// % subject gaps relative to subject length.
    pub perc_db_gap: f32,
}

/// Compute alignment statistics from IUPAC-decoded aligned strings.
///
/// `query_seq` and `subj_seq` are the aligned sequences as ASCII IUPAC
/// (A/C/G/T/- etc.), equal length (gaps are '-').
/// `query_is_ancestral` matches rmblastn's convention: if `true`, the query
/// is the ancestral (repeat consensus) sequence.
pub fn compute_align_stats(query_seq: &[u8], subj_seq: &[u8], query_is_ancestral: bool) -> AlignStats {
    assert_eq!(query_seq.len(), subj_seq.len());

    let (anc, der) = if query_is_ancestral {
        (query_seq, subj_seq)
    } else {
        (subj_seq, query_seq)
    };

    let n = anc.len();
    let mut matches = 0u32;
    let mut mismatches = 0u32;
    let mut q_gaps = 0u32;
    let mut s_gaps = 0u32;
    let mut transi = 0u32;
    let mut transv = 0u32;
    let mut cpg_sites = 0u32;
    let mut cpg_transi = 0.0f64;
    let mut unambig_pairs = 0u32;
    let mut q_len = 0u32;
    let mut d_len = 0u32;

    let mut prev_anc = b' ';
    let mut prev_pair = PairType::Unknown;

    for i in 0..n {
        let ac = anc[i].to_ascii_uppercase();
        let dc = der[i].to_ascii_uppercase();

        // Track query/subject lengths (not counting gaps)
        if query_seq[i] != b'-' { q_len += 1; }
        if subj_seq[i] != b'-' { d_len += 1; }

        // CpG detection: ancestral C followed by G
        if prev_anc == b'C' && ac == b'G' {
            cpg_sites += 1;
            if prev_pair == PairType::Transition && classify_pair_ascii(ac, dc) == PairType::Transition {
                cpg_transi -= 1.0;
            } else if prev_pair == PairType::Transition || classify_pair_ascii(ac, dc) == PairType::Transition {
                cpg_transi -= 0.9;
            }
        }

        if ac == b'-' || dc == b'-' {
            if query_seq[i] == b'-' {
                q_gaps += 1;
            } else {
                s_gaps += 1;
            }
            // NCBI updates prev_anc whenever ancestral is non-gap, even when derived has a gap.
            // Matching that behaviour here prevents false CpG detection in the next iteration.
            if ac != b'-' {
                prev_anc = ac;
                prev_pair = PairType::Gap; // gap in derived → treat as unknown pair type
            }
        } else {
            let pair = classify_pair_ascii(ac, dc);
            match pair {
                PairType::Identity => {
                    matches += 1;
                    if is_unambig(ac) && is_unambig(dc) {
                        unambig_pairs += 1;
                    }
                }
                PairType::Transition => {
                    mismatches += 1;
                    if is_unambig(ac) && is_unambig(dc) {
                        unambig_pairs += 1;
                        transi += 1;
                        cpg_transi += 1.0;
                    }
                }
                PairType::Transversion => {
                    mismatches += 1;
                    if is_unambig(ac) && is_unambig(dc) {
                        unambig_pairs += 1;
                        transv += 1;
                    }
                }
                PairType::IubMatch => { /* treated as match for perc_sub purposes */ }
                PairType::IubSub => { mismatches += 1; /* IUPAC pair with no base overlap: counts as substitution, not trans/transv */ }
                _ => {}
            }
            prev_anc = ac;
            prev_pair = pair;
        }
    }

    let (kdiv, cpg_kdiv) = kimura_divergence(transi, cpg_transi, transv, unambig_pairs);

    let perc_sub = if q_len > 0 { (mismatches as f32 / q_len as f32) * 100.0 } else { 0.0 };
    let perc_query_gap = if q_len > 0 { (q_gaps as f32 / q_len as f32) * 100.0 } else { 0.0 };
    let perc_db_gap = if d_len > 0 { (s_gaps as f32 / d_len as f32) * 100.0 } else { 0.0 };

    AlignStats {
        matches,
        mismatches,
        query_gaps: q_gaps,
        subject_gaps: s_gaps,
        transitions: transi,
        transversions: transv,
        cpg_sites,
        kdiv,
        cpg_kdiv,
        perc_sub,
        perc_query_gap,
        perc_db_gap,
    }
}

fn classify_pair_ascii(a: u8, b: u8) -> PairType {
    if a == b'-' || b == b'-' {
        return PairType::Gap;
    }
    // NCBI checks literal equality FIRST (tabular.cpp x_fillAlignStatsRMBlast: `if(anc==der)
    // match++`), so any identical pair — including ambiguity-vs-same-ambiguity (R/R) or N/N —
    // is an exact match, never a substitution.
    if a.to_ascii_uppercase() == b.to_ascii_uppercase() {
        return PairType::Identity;
    }
    if !is_unambig(a) || !is_unambig(b) {
        // At least one side is an ambiguity code (and the two differ literally).
        // NCBI's `mutType` IUBMatch table lists ONLY pure-base↔ambiguity pairs (e.g. AW, GK)
        // and pure-base↔N (AN, CN, GN, TN) — it has NO ambiguity↔ambiguity or ambiguity↔N
        // entries.  So an IUBMatch (not counted as substitution) requires EXACTLY ONE side to
        // be a pure A/C/G/T base whose value lies in the other side's IUB set.  Ambiguity↔N
        // (e.g. K vs N) and ambiguity↔ambiguity pairs are NOT in the table → substitutions.
        // (The previous symmetric set-overlap rule wrongly treated K-vs-N etc. as matches,
        // undercounting perc_sub on subject-ambiguity × query-N columns.)
        let one_pure = is_unambig(a) ^ is_unambig(b);
        if one_pure && (iupac_bitmask(a) & iupac_bitmask(b)) != 0 {
            return PairType::IubMatch;
        } else {
            return PairType::IubSub;
        }
    }
    let a_pur = a == b'A' || a == b'G';
    let b_pur = b == b'A' || b == b'G';
    if a_pur == b_pur { PairType::Transition } else { PairType::Transversion }
}

/// 4-bit IUPAC bitmask: bit0=A, bit1=C, bit2=G, bit3=T.
#[inline]
fn iupac_bitmask(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => 0b0001,
        b'C' => 0b0010,
        b'G' => 0b0100,
        b'T' | b'U' => 0b1000,
        b'R' => 0b0101, // A,G
        b'Y' => 0b1010, // C,T
        b'M' => 0b0011, // A,C
        b'K' => 0b1100, // G,T
        b'W' => 0b1001, // A,T
        b'S' => 0b0110, // C,G
        b'B' => 0b1110, // C,G,T
        b'D' => 0b1101, // A,G,T
        b'H' => 0b1011, // A,C,T
        b'V' => 0b0111, // A,C,G
        _    => 0b1111, // N, X, '-', unknown → treat as compatible with all
    }
}

#[inline]
fn is_unambig(c: u8) -> bool {
    matches!(c, b'A' | b'C' | b'G' | b'T')
}

/// Kimura two-parameter divergence.
/// Returns (kimura × 100, cpg_adjusted_kimura × 100).
fn kimura_divergence(transi: u32, cpg_transi: f64, transv: u32, unambig: u32) -> (f64, f64) {
    if unambig < 1 {
        return (100.0, 100.0);
    }
    let n = unambig as f64;
    let p = transi as f64 / n;
    let cpg_p = cpg_transi / n;
    let q = transv as f64 / n;

    let k = kimura_formula(p, q);
    let cpg_k = kimura_formula(cpg_p, q);
    (k, cpg_k)
}

fn kimura_formula(p: f64, q: f64) -> f64 {
    let operand = (1.0 - 2.0 * p - q) * (1.0 - 2.0 * q).sqrt();
    if operand > 0.0 {
        (-0.5 * operand.ln()).abs() * 100.0
    } else {
        100.0
    }
}

/// Build IUPAC ASCII strings from BLASTNA encoded aligned sequences.
/// `q_bases` and `s_bases` should be the same length (aligned, with gap=15).
pub fn blastna_to_iupac_aligned(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| {
        if b == 15 { b'-' } else { BLASTNA_TO_IUPAC[(b & 15) as usize] }
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_perfect_match() {
        let q = b"ACGTACGT";
        let s = b"ACGTACGT";
        let st = compute_align_stats(q, s, false);
        assert_eq!(st.mismatches, 0);
        assert_eq!(st.matches, 8);
        assert!((st.kdiv - 0.0).abs() < 1.0);
    }

    #[test]
    fn test_transition() {
        // A→G is a transition
        let q = b"A";
        let s = b"G";
        let st = compute_align_stats(q, s, false);
        assert_eq!(st.transitions, 1);
        assert_eq!(st.transversions, 0);
    }

    #[test]
    fn test_transversion() {
        // A→C is a transversion
        let q = b"A";
        let s = b"C";
        let st = compute_align_stats(q, s, false);
        assert_eq!(st.transitions, 0);
        assert_eq!(st.transversions, 1);
    }
}
