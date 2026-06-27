//! Karlin-Altschul–based ungapped pre-filter cutoff for RMBlast custom matrices.
//!
//! NCBI's `BlastInitialWordParametersNew` (blast_parameters.c) sets the ungapped
//! keep threshold (`cutoff_score`) from KA statistics whenever valid gapped KA
//! params (lambda, K) exist for the matrix + gap-penalty combination, falling
//! back to `min_raw_gapped_score / 2` (= `cutoff_score_max`) only when they do
//! not.  rmblastn enters this path because `BlastInitialWordParametersNew` reads
//! the *gapped* KA block (`kbp_gap`), which IS populated for the RepeatMasker
//! `p##g` matrices from the hardcoded ALP table in blast_stat.c.
//!
//! For the matrix_only_scoring / blastn path the computation reduces to a closed
//! form (no length adjustment, no gap-decay since `gap_decay_rate == 0` for the
//! ungapped pass, and K taken straight from the table — no Karlin K iteration):
//!
//! ```text
//!   E        = CUTOFF_E_BLASTN = 0.05
//!   searchsp = MIN(avg_subj_len, 2*query_len) * avg_subj_len        (Int8)
//!   es       = ceil( ln(K * searchsp / E) / lambda )                (BlastKarlinEtoS_simple)
//!   cutoff   = MIN( MAX(es, 1), min_raw_gapped_score / 2 )          (capped at cutoff_score_max)
//! ```
//!
//! `comparison.matrix` (and any matrix/gap combo not in the table) has no KA
//! entry → `lookup_ka` returns `None` → callers keep the `min_raw_gapped_score/2`
//! fallback, exactly as NCBI does when KA params are unavailable.
//!
//! Stage 1 (this module) hardcodes NCBI's ALP table verbatim.  A later stage may
//! compute the parameters with ALP directly; only the `(lambda, K)` *source*
//! would change — the lookup/formula/logging stay put.

/// `CUTOFF_E_BLASTN` from blast_parameters.h.
pub const CUTOFF_E_BLASTN: f64 = 0.05;

/// Hardcoded RMBlast ALP parameters: `(matrix_name, gap_open, gap_extend, lambda, K)`.
///
/// Transcribed verbatim from `rmblast_*_values` in NCBI blast_stat.c (2.17.0).
/// Each RepeatMasker `p##g` matrix supports exactly one gap-penalty combination.
/// Only `lambda` and `K` are needed for the ungapped cutoff (`BlastKarlinEtoS_simple`
/// uses neither H, alpha, beta, ... ).
#[rustfmt::skip]
pub static RMBLAST_KA_PARAMS: &[(&str, i32, i32, f64, f64)] = &[
    // 14p family (gap 29/6)
    ("14p35g.matrix", 29, 6, 0.1241860392, 0.2565118357),
    ("14p37g.matrix", 29, 6, 0.1243509013, 0.2538327145),
    ("14p39g.matrix", 29, 6, 0.1181081970, 0.2381799563),
    ("14p41g.matrix", 29, 6, 0.1209164962, 0.2516961980),
    ("14p43g.matrix", 29, 6, 0.1204338701, 0.2484920137),
    ("14p45g.matrix", 29, 6, 0.1281012415, 0.2651519027),
    ("14p47g.matrix", 29, 6, 0.1169571490, 0.2334091021),
    ("14p49g.matrix", 29, 6, 0.1167697593, 0.2289233061),
    ("14p51g.matrix", 29, 6, 0.1238779017, 0.2467119379),
    ("14p53g.matrix", 29, 6, 0.1211022467, 0.2307304663),
    // 18p family (gap 28/5)
    ("18p35g.matrix", 28, 5, 0.1204776630, 0.2111410357),
    ("18p37g.matrix", 28, 5, 0.1213800379, 0.2057492194),
    ("18p39g.matrix", 28, 5, 0.1228853168, 0.2111281814),
    ("18p41g.matrix", 28, 5, 0.1096830298, 0.1715966796),
    ("18p43g.matrix", 28, 5, 0.1179924113, 0.1957371364),
    ("18p45g.matrix", 28, 5, 0.1184247414, 0.2068045787),
    ("18p47g.matrix", 28, 5, 0.1158783320, 0.1842046578),
    ("18p49g.matrix", 28, 5, 0.1248964218, 0.2141853759),
    ("18p51g.matrix", 28, 5, 0.1106793836, 0.1621592148),
    ("18p53g.matrix", 28, 5, 0.1082347513, 0.1408677318),
    // 20p family (gap 25/5)
    ("20p35g.matrix", 25, 5, 0.1112832877, 0.1466667286),
    ("20p37g.matrix", 25, 5, 0.1134511409, 0.1501236494),
    ("20p39g.matrix", 25, 5, 0.1233401001, 0.1853070864),
    ("20p41g.matrix", 25, 5, 0.1088138918, 0.1409300313),
    ("20p43g.matrix", 25, 5, 0.1097495426, 0.1408421343),
    ("20p45g.matrix", 25, 5, 0.1101110565, 0.1402262486),
    ("20p47g.matrix", 25, 5, 0.1198508511, 0.1752301830),
    ("20p49g.matrix", 25, 5, 0.1150022018, 0.1500562415),
    ("20p51g.matrix", 25, 5, 0.1101356126, 0.1310353805),
    ("20p53g.matrix", 25, 5, 0.1239625969, 0.1752022586),
    // 25p family (gap 22/5)
    ("25p35g.matrix", 22, 5, 0.0965155465, 0.0677355511),
    ("25p37g.matrix", 22, 5, 0.0989264764, 0.0763719670),
    ("25p39g.matrix", 22, 5, 0.1005843197, 0.0805691511),
    ("25p41g.matrix", 22, 5, 0.1030661339, 0.0840485952),
    ("25p43g.matrix", 22, 5, 0.1155434360, 0.1226554129),
    ("25p45g.matrix", 22, 5, 0.0973402625, 0.0742020190),
    ("25p47g.matrix", 22, 5, 0.0967346256, 0.0743885008),
    ("25p49g.matrix", 22, 5, 0.1070577141, 0.0941947942),
    ("25p51g.matrix", 22, 5, 0.1044972455, 0.0789036639),
    ("25p53g.matrix", 22, 5, 0.0977828435, 0.0607500307),
    // 30p family (gap 22/5) — NCBI ships placeholder ALP values for this one.
    ("30p53g.matrix", 22, 5, 0.1,          0.1),
];

/// Details of a KA-derived cutoff, for run logging.
#[derive(Debug, Clone, Copy)]
pub struct KaCutoffInfo {
    pub lambda: f64,
    pub k: f64,
    pub searchsp: u64,
    /// The raw KA score `es` before the `cutoff_score_max` cap.
    pub es: i32,
    pub avg_subj_length: u64,
}

/// Look up `(lambda, K)` for a matrix path + gap-penalty combination.
///
/// Matching mirrors NCBI's `strcasecmp(matrix_info->name, matrix)`: the matrix
/// *base name* (e.g. `18p43g.matrix`) is compared case-insensitively, and the
/// gap-open / gap-extend penalties must match the table entry exactly.
pub fn lookup_ka(matrix_name: &str, gap_open: i32, gap_extend: i32) -> Option<(f64, f64)> {
    let base = std::path::Path::new(matrix_name)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(matrix_name);
    RMBLAST_KA_PARAMS
        .iter()
        .find(|&&(n, o, e, _, _)| o == gap_open && e == gap_extend && n.eq_ignore_ascii_case(base))
        .map(|&(_, _, _, lambda, k)| (lambda, k))
}

/// Compute the ungapped pre-filter cutoff exactly as NCBI's
/// `BlastInitialWordParametersNew` does for the blastn / matrix_only path.
///
/// Returns `(cutoff, Some(info))` when KA params matched, else
/// `(min_raw_gapped_score / 2, None)` — the historical fixed fallback.
///
/// `query_length` is the single-strand query length (it is doubled here for the
/// reverse-complement strand, matching NCBI's `query_length *= 2`).
/// `avg_subj_length` is the average DB sequence length (`BlastSeqSrcGetAvgSeqLen`,
/// i.e. total residues / number of sequences, integer-truncated).
pub fn ungapped_cutoff(
    matrix_name: &str,
    gap_open: i32,
    gap_extend: i32,
    min_raw_gapped_score: i32,
    query_length: u64,
    avg_subj_length: u64,
) -> (i32, Option<KaCutoffInfo>) {
    let fallback = min_raw_gapped_score / 2;
    match lookup_ka(matrix_name, gap_open, gap_extend) {
        Some((lambda, k)) if avg_subj_length > 0 && lambda > 0.0 && k > 0.0 => {
            // searchsp = MIN(subj_len, 2*query_len) * subj_len   (Int8 in NCBI)
            let qlen2 = query_length.saturating_mul(2);
            let searchsp = avg_subj_length.min(qlen2).saturating_mul(avg_subj_length);
            // BlastKarlinEtoS_simple: S = ceil( ln(K * searchsp / E) / lambda )
            let es = ((k * searchsp as f64 / CUTOFF_E_BLASTN).ln() / lambda).ceil() as i32;
            // BLAST_Cutoffs keeps MAX(es, 1); BlastInitialWordParametersNew then
            // caps at cutoff_score_max (= min_raw_gapped_score / 2).
            let cutoff = es.max(1).min(fallback);
            (
                cutoff,
                Some(KaCutoffInfo { lambda, k, searchsp, es, avg_subj_length }),
            )
        }
        _ => (fallback, None),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lookup_matches_basename_and_gaps_case_insensitively() {
        assert_eq!(lookup_ka("18p43g.matrix", 28, 5), Some((0.1179924113, 0.1957371364)));
        assert_eq!(lookup_ka("/some/path/18p43g.matrix", 28, 5), Some((0.1179924113, 0.1957371364)));
        assert_eq!(lookup_ka("18P43G.MATRIX", 28, 5), Some((0.1179924113, 0.1957371364)));
        // wrong gap combo -> no match (NCBI Mode 3 -> invalid -> fallback)
        assert_eq!(lookup_ka("18p43g.matrix", 20, 5), None);
        // not in the table (e.g. comparison.matrix) -> no match
        assert_eq!(lookup_ka("comparison.matrix", 28, 5), None);
    }

    #[test]
    fn cutoff_falls_back_without_ka() {
        // comparison.matrix is not in the KA table -> fixed min_raw_gapped_score/2.
        let (c, info) = ungapped_cutoff("comparison.matrix", 20, 5, 200, 1_000_000, 500);
        assert_eq!(c, 100);
        assert!(info.is_none());
    }

    #[test]
    fn cutoff_ka_based_18p43g_below_fixed() {
        // The bug #35 case: 18p43g {28,5}, min_raw_gapped_score 300 -> fixed would be 150.
        // KA cutoff must be < 150 (so the ungapped-140 seed survives).
        let (c, info) = ungapped_cutoff("18p43g.matrix", 28, 5, 300, 49802, 919);
        let info = info.expect("KA should apply for 18p43g {28,5}");
        // closed form check
        let searchsp = 919u64.min(2 * 49802).saturating_mul(919);
        let es = ((0.1957371364 * searchsp as f64 / 0.05).ln() / 0.1179924113).ceil() as i32;
        assert_eq!(info.es, es);
        assert_eq!(c, es.max(1).min(150));
        assert!(c < 150, "KA cutoff {} must be below the fixed 150", c);
        assert!(c <= 140, "KA cutoff {} must let the ungapped-140 anchor seed through", c);
    }
}
