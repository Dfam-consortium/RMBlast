//! Output formatting for rmblastn: tabular (-outfmt 6) and pairwise (-outfmt 0, default).
//!
//! Tabular output (outfmt 6):
//!   Plus strand:  sstrand="plus"  sstart=1-based-start  send=1-based-end
//!   Minus strand: sstrand="minus" sstart=1-based-right-end  send=1-based-left-start
//!   (sstart > send for minus — standard NCBI blastn tabular convention)
//!
//!   Note: NCBI's tabular.cpp swaps the transi/transv/cpg_sites member names, so
//!   the outfmt field "transi" maps to st.transversions and "cpg_sites" maps to
//!   st.transitions.  This port preserves that quirk for output compatibility.
//!
//! Pairwise output (outfmt 0):
//!   Matches NCBI rmblastn default pairwise alignment format exactly.

use std::io::Write;

use crate::hits::{Hsp, Strand};
use crate::stats::{blastna_to_iupac_aligned, AlignStats};

// ── Program identity string ───────────────────────────────────────────────────
pub const RMBLASTN_VERSION: &str = "RMBLASTN 2.17.0+";

// ── AlignResult ───────────────────────────────────────────────────────────────

/// A fully decorated alignment result ready for output.
#[derive(Clone)]
pub struct AlignResult {
    pub hsp: Hsp,
    pub query_id: String,
    pub subject_id: String,
    pub stats: AlignStats,
}

// ── Tabular output ────────────────────────────────────────────────────────────

/// A recognized output column in -outfmt 6 field list.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OutField {
    Score,
    PercSub,
    PercQueryGap,
    PercDbGap,
    QSeqId,
    QStart,
    QEnd,
    QLen,
    Sstrand,
    SSeqId,
    SStart,
    SSend,
    SLen,
    Kdiv,
    CpgKdiv,
    /// NCBI quirk: outfmt field "transi" outputs st.transversions.
    Transi,
    /// NCBI quirk: outfmt field "transv" outputs st.cpg_sites.
    Transv,
    /// NCBI quirk: outfmt field "cpg_sites" outputs st.transitions.
    CpgSites,
    QSeq,
    SSeq,
    /// Expect value (requires Karlin-Altschul statistics; prints "1.0" when
    /// unavailable, matching the NCBI Mode-3 sentinel).
    Evalue,
    /// Bit score (prints " 0.0" when statistics are unavailable).
    Bitscore,
}

impl OutField {
    pub fn from_str(s: &str) -> Option<Self> {
        match s {
            "score"          => Some(Self::Score),
            "perc_sub"       => Some(Self::PercSub),
            "perc_query_gap" => Some(Self::PercQueryGap),
            "perc_db_gap"    => Some(Self::PercDbGap),
            "qseqid"         => Some(Self::QSeqId),
            "qstart"         => Some(Self::QStart),
            "qend"           => Some(Self::QEnd),
            "qlen"           => Some(Self::QLen),
            "sstrand"        => Some(Self::Sstrand),
            "sseqid"         => Some(Self::SSeqId),
            "sstart"         => Some(Self::SStart),
            "send"           => Some(Self::SSend),
            "slen"           => Some(Self::SLen),
            "kdiv"           => Some(Self::Kdiv),
            "cpg_kdiv"       => Some(Self::CpgKdiv),
            "transi"         => Some(Self::Transi),
            "transv"         => Some(Self::Transv),
            "cpg_sites"      => Some(Self::CpgSites),
            "qseq"           => Some(Self::QSeq),
            "sseq"           => Some(Self::SSeq),
            "evalue"         => Some(Self::Evalue),
            "bitscore"       => Some(Self::Bitscore),
            _                => None,
        }
    }
}

/// True if any requested field needs Karlin-Altschul statistics — the gate
/// that keeps `ka_stats::resolve` (including a potential ALP fit) and all
/// per-HSP E-value math from running when the output doesn't ask for them.
pub fn outfmt_needs_stats(fields: &[OutField]) -> bool {
    fields.iter().any(|f| matches!(f, OutField::Evalue | OutField::Bitscore))
}

/// Parse an outfmt string like "6 score perc_sub ... qseq sseq" into a field list.
/// The leading format number (e.g. "6") is skipped.
pub fn parse_outfmt(s: &str) -> Vec<OutField> {
    s.split_whitespace()
        .filter(|t| t.parse::<u32>().is_err())
        .filter_map(|t| OutField::from_str(t))
        .collect()
}

/// C `%e`-style scientific formatting: signed exponent, minimum two exponent
/// digits ("1.23e-05"), unlike Rust's `{:e}` ("1.23e-5").
fn c_sci(value: f64, precision: usize) -> String {
    let s = format!("{:.*e}", precision, value);
    match s.split_once('e') {
        Some((mant, exp)) => {
            let exp: i32 = exp.parse().unwrap_or(0);
            let sign = if exp < 0 { '-' } else { '+' };
            format!("{}e{}{:02}", mant, sign, exp.abs())
        }
        None => s,
    }
}

/// Format an E-value exactly as NCBI's tabular output does
/// (CBlastTabularInfo::SetScores: CAlignFormatUtil::GetScoreString, with the
/// whole [1e-180, 0.0009) band overridden by NStr::DoubleToString(evalue, 2,
/// fDoubleScientific) = C "%.2e").  The Mode-3 sentinel 1.0 prints as "1.0".
pub fn format_evalue(evalue: f64) -> String {
    if evalue < 1.0e-180 {
        "0.0".to_string()
    } else if evalue < 0.0009 {
        c_sci(evalue, 2)                  // tabular override of %2.0le/%3.0le
    } else if evalue < 0.1 {
        format!("{:.3}", evalue)          // %4.3lf (width always exceeded)
    } else if evalue < 1.0 {
        format!("{:.2}", evalue)          // %3.2lf
    } else if evalue < 10.0 {
        format!("{:.1}", evalue)          // %2.1lf
    } else {
        format!("{:.0}", evalue)          // %2.0lf
    }
}

/// Format a bit score exactly as CAlignFormatUtil::GetScoreString does.  Note
/// the truncating (long) cast for the >99.9 arm and the width-4 space padding
/// of %4.1lf (the sentinel 0.0 prints as " 0.0").
pub fn format_bit_score(bits: f64) -> String {
    if bits > 99999.0 {
        c_sci(bits, 3)                    // %5.3le (width always exceeded)
    } else if bits > 99.9 {
        format!("{:>3}", bits as i64)     // %3.0ld — truncation toward zero
    } else {
        format!("{:>4.1}", bits)          // %4.1lf
    }
}

/// Write a single HSP in NCBI rmblastn outfmt-6 style (tab-separated).
///
/// `ka` supplies Karlin-Altschul statistics for the `evalue`/`bitscore`
/// fields; pass `None` when those fields are absent (or to force the NCBI
/// Mode-3 sentinel rendering: evalue "1.0", bitscore " 0.0").
pub fn write_tabular<W: Write>(
    w: &mut W,
    r: &AlignResult,
    delimiter: char,
    fields: &[OutField],
    ka: Option<&rmstats::RmStats>,
) -> std::io::Result<()> {
    let h  = &r.hsp;
    let st = &r.stats;

    let sstart_1b = h.s_start + 1;
    let send_1b   = h.s_end;

    let (ss, se) = match h.strand {
        Strand::Plus  => (sstart_1b, send_1b),
        Strand::Minus => (send_1b, sstart_1b),
    };

    let sstrand_str = match h.strand {
        Strand::Plus  => "plus",
        Strand::Minus => "minus",
    };

    let mut first = true;
    for field in fields {
        if !first { write!(w, "{}", delimiter)?; }
        first = false;

        match field {
            OutField::Score        => write!(w, "{}", h.score)?,
            OutField::PercSub      => write!(w, "{:.2}", st.perc_sub)?,
            OutField::PercQueryGap => write!(w, "{:.2}", st.perc_query_gap)?,
            OutField::PercDbGap    => write!(w, "{:.2}", st.perc_db_gap)?,
            OutField::QSeqId       => write!(w, "{}", r.query_id)?,
            OutField::QStart       => write!(w, "{}", h.q_start + 1)?,
            OutField::QEnd         => write!(w, "{}", h.q_end)?,
            OutField::QLen         => write!(w, "{}", h.q_len)?,
            OutField::Sstrand      => write!(w, "{}", sstrand_str)?,
            OutField::SSeqId       => write!(w, "{}", r.subject_id)?,
            OutField::SStart       => write!(w, "{}", ss)?,
            OutField::SSend        => write!(w, "{}", se)?,
            OutField::SLen         => write!(w, "{}", h.s_len)?,
            OutField::Kdiv         => write!(w, "{:.2}", st.kdiv)?,
            OutField::CpgKdiv      => write!(w, "{:.2}", st.cpg_kdiv)?,
            OutField::Transi       => write!(w, "{}", st.transversions)?,
            OutField::Transv       => write!(w, "{}", st.cpg_sites)?,
            OutField::CpgSites     => write!(w, "{}", st.transitions)?,
            OutField::QSeq => {
                let iupac = blastna_to_iupac_aligned(&h.q_seq);
                w.write_all(&iupac)?;
            }
            OutField::SSeq => {
                let iupac = blastna_to_iupac_aligned(&h.s_seq);
                w.write_all(&iupac)?;
            }
            // The final (complexity-adjusted) score is what NCBI computes
            // E-values from at traceback time — h.score is exactly that.
            OutField::Evalue => {
                let e = ka.map_or(1.0, |s| s.evalue(h.score));
                write!(w, "{}", format_evalue(e))?;
            }
            OutField::Bitscore => {
                let b = ka.map_or(0.0, |s| s.bit_score(h.score));
                write!(w, "{}", format_bit_score(b))?;
            }
        }
    }
    write!(w, "\n")?;
    Ok(())
}

// ── Pairwise output ───────────────────────────────────────────────────────────

/// NCBI's GetPercentMatch formula.
/// Returns 100 only for a perfect match; otherwise rounds 100*n/d, capped at 99.
fn get_percent_match(n: u32, d: u32) -> u32 {
    if d == 0 { return 0; }
    if n == d { return 100; }
    let r = (0.5 + 100.0 * n as f64 / d as f64) as u32;
    r.min(99)
}

/// Number of decimal digits needed to represent `n`.
fn num_digits(n: u32) -> usize {
    if n == 0 { 1 } else { n.ilog10() as usize + 1 }
}

/// Convert a Unix timestamp (seconds since 1970-01-01 UTC) to an NCBI date string.
/// Format: "Mon D, YYYY  H:MM AM/PM" (no leading zeros on day or hour).
fn format_ncbi_date(secs: u64) -> String {
    const MONTH_NAMES: [&str; 12] = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];

    let time_of_day = secs % 86400;
    let minute = (time_of_day / 60) % 60;
    let hour   = time_of_day / 3600;

    // Howard Hinnant's civil_from_days (public domain).
    let z = secs as i64 / 86400 + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = z - era * 146_097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146_096) / 365;
    let y_raw = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp  = (5 * doy + 2) / 153;
    let d   = doy - (153 * mp + 2) / 5 + 1;
    let m   = if mp < 10 { mp + 3 } else { mp - 9 };
    let y   = if m <= 2 { y_raw + 1 } else { y_raw };

    let (h12, ampm) = if hour == 0       { (12u64, "AM") }
                      else if hour < 12  { (hour, "AM") }
                      else if hour == 12 { (12, "PM") }
                      else               { (hour - 12, "PM") };

    format!("{} {}, {}  {}:{:02} {}",
        MONTH_NAMES[(m - 1) as usize], d, y, h12, minute, ampm)
}

/// Write the program/database header (printed once, before any queries).
///
/// ```text
/// RMBLASTN VERSION
///
///
/// Reference: ...
/// ...
///
///
///
/// Database: db_name
///            N sequences; M total letters
///
///
///
/// ```
pub fn write_pairwise_program_header<W: Write>(
    w: &mut W,
    db_name: &str,
    n_seqs: usize,
    total_letters: u64,
) -> std::io::Result<()> {
    write!(w, "{}\n\n\n", RMBLASTN_VERSION)?;
    write!(w, "Reference: Robert M. Hubley, Arian Smit\n")?;
    write!(w, "RMBlast - RepeatMasker Search Engine\n")?;
    write!(w, "2010 <http://www.repeatmasker.org>\n")?;
    write!(w, "\n")?;
    write!(w, "Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A.\n")?;
    write!(w, "Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J.\n")?;
    write!(w, "Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of\n")?;
    write!(w, "protein database search programs\", Nucleic Acids Res. 25:3389-3402.\n")?;
    write!(w, "\n\n\n")?;
    write!(w, "Database: {}\n", db_name)?;
    write!(w, "           {} sequences; {} total letters\n", n_seqs, total_letters)?;
    write!(w, "\n\n\n")?;
    Ok(())
}

/// Write the per-query header block.
///
/// ```text
/// Query= defline
///
/// Length=N
///
/// ```
pub fn write_pairwise_query_header<W: Write>(
    w: &mut W,
    query_defline: &str,
    query_len: u32,
) -> std::io::Result<()> {
    write!(w, "Query= {}\n\nLength={}\n\n", query_defline, query_len)?;
    Ok(())
}

/// Write the per-subject header block (before the first HSP for a subject).
///
/// ```text
/// >subject_id
/// Length=N
///
/// ```
pub fn write_pairwise_subject_header<W: Write>(
    w: &mut W,
    subject_id: &str,
    subject_len: u32,
) -> std::io::Result<()> {
    write!(w, ">{} \nLength={}\n\n", subject_id, subject_len)?;
    Ok(())
}

/// Write one HSP in pairwise alignment format.
///
/// Outputs the stats block, then the wrapped alignment rows (60 bases/row),
/// then two blank lines (as required by NCBI's between-HSP and end-of-subject
/// blank-line conventions).
pub fn write_pairwise_hsp<W: Write>(w: &mut W, r: &AlignResult) -> std::io::Result<()> {
    let h  = &r.hsp;
    let st = &r.stats;

    let aln_len = h.q_seq.len() as u32;
    let n_ident = st.matches;
    let n_gaps  = st.query_gaps + st.subject_gaps;

    let pct_ident = get_percent_match(n_ident, aln_len);
    let pct_gaps  = get_percent_match(n_gaps,  aln_len);

    let strand_str = match h.strand {
        Strand::Plus  => "Plus/Plus",
        Strand::Minus => "Plus/Minus",
    };

    write!(w, " Score = {}\n", h.score)?;
    write!(w, " Substitutions = {:.2}%, Query Gaps = {:.2}%, DB Gaps = {:.2}%\n",
        st.perc_sub, st.perc_query_gap, st.perc_db_gap)?;
    write!(w, " TransI/TransV = {}/{},  CpG_sites = {}\n",
        st.transitions, st.transversions, st.cpg_sites)?;
    write!(w, " Kimura Div = {:.2}%, Kimura Div CpG Adjusted = {:.2}%\n",
        st.kdiv, st.cpg_kdiv)?;
    write!(w, " Identities = {}/{} ({}%), Gaps = {}/{} ({}%)\n",
        n_ident, aln_len, pct_ident, n_gaps, aln_len, pct_gaps)?;
    write!(w, " Strand={}\n\n", strand_str)?;

    // ── Alignment rows ────────────────────────────────────────────────────────
    let q_iupac = blastna_to_iupac_aligned(&h.q_seq);
    let s_iupac = blastna_to_iupac_aligned(&h.s_seq);

    // Column width for position numbers: digits of the maximum coordinate.
    let max_coord = h.q_end.max(h.s_end);
    let pos_w = num_digits(max_coord);

    // Starting display positions (1-based).
    let mut q_pos = h.q_start + 1;
    let mut s_pos: u32 = match h.strand {
        Strand::Plus  => h.s_start + 1,
        Strand::Minus => h.s_end,        // rightmost 1-based position (counts down)
    };

    let total_len = q_iupac.len();
    let mut i = 0;
    while i < total_len {
        let end = (i + 60).min(total_len);
        let q_row = &q_iupac[i..end];
        let s_row = &s_iupac[i..end];

        let q_nongap = q_row.iter().filter(|&&c| c != b'-').count() as u32;
        let s_nongap = s_row.iter().filter(|&&c| c != b'-').count() as u32;

        let q_row_start = q_pos;
        let q_row_end   = if q_nongap > 0 { q_pos + q_nongap - 1 } else { q_pos };

        let (s_row_start, s_row_end) = match h.strand {
            Strand::Plus  => {
                let ss = s_pos;
                let se = if s_nongap > 0 { s_pos + s_nongap - 1 } else { s_pos };
                (ss, se)
            }
            Strand::Minus => {
                let ss = s_pos;
                let se = if s_nongap > 0 { s_pos - (s_nongap - 1) } else { s_pos };
                (ss, se)
            }
        };

        // Build the match indicator line.
        let match_line: Vec<u8> = q_row.iter().zip(s_row.iter()).map(|(&q, &s)| {
            if q == b'-' || s == b'-' { b' ' }
            else if q.to_ascii_uppercase() == s.to_ascii_uppercase() { b'|' }
            else { b' ' }
        }).collect();

        // Prefix width: "Query" (5) + 2 spaces + pos_w digits + 2 spaces = 9 + pos_w.
        let prefix_spaces = 9 + pos_w;

        write!(w, "Query  {:<width$}  ", q_row_start, width = pos_w)?;
        w.write_all(q_row)?;
        write!(w, "  {}\n", q_row_end)?;

        for _ in 0..prefix_spaces { w.write_all(b" ")?; }
        w.write_all(&match_line)?;
        write!(w, "\n")?;

        write!(w, "Sbjct  {:<width$}  ", s_row_start, width = pos_w)?;
        w.write_all(s_row)?;
        write!(w, "  {}\n", s_row_end)?;

        write!(w, "\n")?;   // blank line after each chunk

        q_pos += q_nongap;
        match h.strand {
            Strand::Plus  => s_pos += s_nongap,
            Strand::Minus => s_pos = s_pos.saturating_sub(s_nongap),
        }

        i = end;
    }

    // One extra blank line → total 2 blank lines after last chunk.
    write!(w, "\n")?;

    Ok(())
}

/// Write the program/database footer (printed once, after all queries).
///
/// ```text
///
///
///   Database: db_name
///     Posted date:  Mon D, YYYY  H:MM AM/PM
///   Number of letters in database: M
///   Number of sequences in database:  N
///
///
///
/// Matrix: blastn matrix 0 0
/// Gap Penalties: Existence: GAP_OPEN, Extension: GAP_EXTEND.00
/// ```
pub fn write_pairwise_footer<W: Write>(
    w: &mut W,
    db_name: &str,
    n_seqs: usize,
    total_letters: u64,
    db_mtime_unix: Option<u64>,
    gap_open: i32,
    gap_extend: i32,
) -> std::io::Result<()> {
    // Two extra blank lines (the last HSP already contributed 2, giving 4 total).
    write!(w, "\n\n")?;
    write!(w, "  Database: {}\n", db_name)?;
    if let Some(secs) = db_mtime_unix {
        write!(w, "    Posted date:  {}\n", format_ncbi_date(secs))?;
    }
    write!(w, "  Number of letters in database: {}\n", total_letters)?;
    write!(w, "  Number of sequences in database:  {}\n", n_seqs)?;
    write!(w, "\n\n\n")?;
    write!(w, "Matrix: blastn matrix 0 0\n")?;
    write!(w, "Gap Penalties: Existence: {}, Extension: {:.2}\n",
        gap_open, gap_extend as f64)?;
    Ok(())
}

/// Write all pairwise results for one query.  Results must be pre-sorted by score
/// descending.  Subjects are output in order of first appearance (best-score-first).
pub fn write_pairwise_results<W: Write>(
    w: &mut W,
    results: &[AlignResult],
    query_defline: &str,
    query_len: u32,
) -> std::io::Result<()> {
    if results.is_empty() {
        return Ok(());
    }

    write_pairwise_query_header(w, query_defline, query_len)?;

    // Group by subject, preserving order of first appearance.
    let mut seen_subjects: Vec<&str> = Vec::new();
    for r in results {
        if !seen_subjects.contains(&r.subject_id.as_str()) {
            seen_subjects.push(&r.subject_id);
        }
    }

    for subj_id in seen_subjects {
        let subj_results: Vec<&AlignResult> = results.iter()
            .filter(|r| r.subject_id == subj_id)
            .collect();
        if subj_results.is_empty() { continue; }

        let s_len = subj_results[0].hsp.s_len;
        write_pairwise_subject_header(w, subj_id, s_len)?;

        for r in subj_results {
            write_pairwise_hsp(w, r)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    // Expected strings mirror C snprintf with the NCBI tabular format rules
    // (SetScores: GetScoreString + the %.2e override for [1e-180, 0.0009)).
    #[test]
    fn evalue_formatting_matches_ncbi_tabular() {
        assert_eq!(format_evalue(1.0e-200), "0.0");
        assert_eq!(format_evalue(2.5e-105), "2.50e-105");
        assert_eq!(format_evalue(1.234e-30), "1.23e-30");
        assert_eq!(format_evalue(5.0e-5), "5.00e-05");
        assert_eq!(format_evalue(8.9e-4), "8.90e-04");
        assert_eq!(format_evalue(0.0009), "0.001");
        assert_eq!(format_evalue(0.05), "0.050");
        assert_eq!(format_evalue(0.55), "0.55");
        assert_eq!(format_evalue(1.0), "1.0");    // Mode-3 sentinel rendering
        assert_eq!(format_evalue(3.7), "3.7");
        assert_eq!(format_evalue(15.0), "15");
        assert_eq!(format_evalue(2000.0), "2000");
    }

    #[test]
    fn bit_score_formatting_matches_ncbi() {
        assert_eq!(format_bit_score(0.0), " 0.0");    // sentinel: %4.1lf pads
        assert_eq!(format_bit_score(7.25), " 7.2");
        assert_eq!(format_bit_score(45.67), "45.7");
        assert_eq!(format_bit_score(99.9), "99.9");
        assert_eq!(format_bit_score(99.95), " 99");   // (long) truncation arm
        assert_eq!(format_bit_score(234.7), "234");
        assert_eq!(format_bit_score(12345.9), "12345");
        assert_eq!(format_bit_score(123456.0), "1.235e+05");
    }
}
