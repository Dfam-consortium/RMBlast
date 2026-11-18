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
//!   NCBI rmblastn 2.17.1's default pairwise format, for a database built with
//!   -parse_seqids (as RepeatMasker and RepeatModeler build them).  Checked
//!   line for line against NCBI on the CLI fixtures; the one known difference
//!   is "Posted date", which NCBI takes from the database build time and this
//!   port from the database file's modification time, in UTC.

use std::io::Write;

use crate::hits::{Hsp, Strand};
use crate::stats::{blastna_to_iupac_aligned, AlignStats};

// ── Program identity string ───────────────────────────────────────────────────
pub const RMBLASTN_VERSION: &str = "RMBLASTN 2.17.1+";

/// `NStr::Int8ToString(n, fWithCommas)`: 11207 -> "11,207".
pub fn with_commas(n: u64) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// `NStr::Wrap(str, width, fWrap_FlatFile)` (corelib/ncbistr.cpp) without a
/// prefix, for text with no backspaces: breaks after a space, a comma or a
/// hyphen, preferring spaces, and forces a break at `width` when there is
/// none.  Used for the database title (width 68, kFormatLineLength).
pub fn ncbi_wrap_flat(s: &str, width: usize) -> Vec<String> {
    #[derive(PartialEq, PartialOrd, Clone, Copy)]
    enum Score { Forced, Punct, Comma, Space, Newline }
    let b = s.as_bytes();
    let len = b.len();
    let mut lines = Vec::new();
    let mut pos = 0;
    while pos < len {
        let mut column = 0usize;
        let mut best_pos = usize::MAX;
        let mut best = Score::Forced;
        let nl = s[pos..].find('\n').map_or(len, |i| pos + i);
        let pos0 = if nl - pos <= width { nl } else { pos };
        let mut pos2 = pos0;
        while pos2 < len && column <= width {
            let c = b[pos2];
            let mut score = Score::Forced;
            let mut score_pos = pos2;
            if c == b'\n' {
                best_pos = pos2;
                best = Score::Newline;
                break;
            } else if c.is_ascii_whitespace() {
                score = Score::Space;
            } else if c == b',' && column < width && score_pos < len - 1 {
                score = Score::Comma;
                score_pos += 1;
            } else if c == b'-' && score_pos < len - 1 && column < width {
                score = Score::Punct;
                score_pos += 1;
            }
            if score >= best && score_pos > pos0 {
                best_pos = score_pos;
                best = score;
            }
            pos2 += 1;
            column += 1;
        }
        if best != Score::Newline && column <= width {
            best_pos = len; // the rest fits
        }
        let best_pos = best_pos.min(len).max(pos + 1);
        lines.push(s[pos..best_pos].to_string());
        pos = best_pos;
        if best == Score::Space {
            while pos < len && b[pos] == b' ' {
                pos += 1;
            }
            if pos < len && b[pos] == b'\n' {
                pos += 1;
            }
        } else if best == Score::Newline {
            pos += 1;
        }
    }
    if lines.is_empty() {
        lines.push(String::new());
    }
    lines
}

/// `s_WrapOutputLine` (objtools/align_format/showalign.cpp), used for the
/// subject title: past every 60th character, break after the next whitespace
/// (which stays on the line).  A long word is never broken.
fn ncbi_wrap_title(s: &str) -> String {
    const LINE_LEN: usize = 60;
    if s.len() <= LINE_LEN {
        return s.to_string();
    }
    let mut out = String::with_capacity(s.len() + s.len() / LINE_LEN);
    let mut do_wrap = false;
    for (i, c) in s.char_indices() {
        if i > 0 && i % LINE_LEN == 0 {
            do_wrap = true;
        }
        out.push(c);
        if do_wrap && c.is_ascii_whitespace() {
            out.push('\n');
            do_wrap = false;
        }
    }
    out
}

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

/// The field names `-outfmt 6` accepts, in the order the help lists them.
pub const OUT_FIELD_NAMES: &[&str] = &[
    "score", "perc_sub", "perc_query_gap", "perc_db_gap", "qseqid", "qstart", "qend",
    "qlen", "sstrand", "sseqid", "sstart", "send", "slen", "kdiv", "cpg_kdiv", "transi",
    "transv", "cpg_sites", "qseq", "sseq", "evalue", "bitscore",
];

impl OutField {
    pub fn parse(s: &str) -> Option<Self> {
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

/// Parse a tabular outfmt string, "6 score perc_sub ... qseq sseq", into its
/// field list.  Anything this port does not implement is an error: another
/// format number, a bare "6" (NCBI's default 12 columns are not implemented),
/// or an unknown field name.  3.0.12 dropped unknown fields silently, shifting
/// the columns, and treated "7 ..." and "10 ..." as 6.
pub fn parse_outfmt(s: &str) -> Result<Vec<OutField>, String> {
    let mut tokens = s.split_whitespace();
    match tokens.next() {
        Some("6") => {}
        Some(fmt) => return Err(format!(
            "output format '{fmt}' is not supported; use 0 (pairwise) or 6 <fields> (tabular)")),
        None => return Err("empty output format".into()),
    }
    let fields: Vec<&str> = tokens.collect();
    if fields.is_empty() {
        return Err(format!(
            "-outfmt 6 needs a field list; NCBI's default columns are not supported. Fields: {}",
            OUT_FIELD_NAMES.join(" ")));
    }
    fields.iter().map(|f| OutField::parse(f).ok_or_else(|| format!(
        "unknown -outfmt field '{f}'. Supported fields: {}", OUT_FIELD_NAMES.join(" "))))
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
    writeln!(w)?;
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
    writeln!(w, "Reference: Robert M. Hubley, Arian Smit")?;
    writeln!(w, "RMBlast - RepeatMasker Search Engine")?;
    writeln!(w, "2010 <http://www.repeatmasker.org>")?;
    writeln!(w)?;
    writeln!(w, "Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A.")?;
    writeln!(w, "Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J.")?;
    writeln!(w, "Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new generation of")?;
    writeln!(w, "protein database search programs\", Nucleic Acids Res. 25:3389-3402.")?;
    write!(w, "\n\n\n")?;
    write!(w, "Database: ")?;
    for line in ncbi_wrap_flat(db_name, 68) {
        writeln!(w, "{}", line)?;
    }
    writeln!(w, "           {} sequences; {} total letters",
        with_commas(n_seqs as u64), with_commas(total_letters))?;
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
    // NCBI trims a query's defline (CFastaReader), unlike a subject title
    write!(w, "Query= {}\n\nLength={}\n\n", query_defline.trim_end(), query_len)?;
    Ok(())
}

/// Write the per-subject header block (before the first HSP for a subject).
///
/// ```text
/// >subject_id title
/// Length=N
///
/// ```
///
/// As NCBI prints it for a -parse_seqids database: the id, one space, then
/// the rest of the defline (after the id and one whitespace character)
/// wrapped by `ncbi_wrap_title`.  An id-only defline gives ">id " with the
/// trailing space, as in NCBI.
pub fn write_pairwise_subject_header<W: Write>(
    w: &mut W,
    subject_id: &str,
    subject_defline: Option<&str>,
    subject_len: u32,
) -> std::io::Result<()> {
    let title = subject_defline
        .and_then(|d| d.strip_prefix(subject_id))
        .map(|rest| { let mut c = rest.chars(); c.next(); c.as_str() })
        .unwrap_or("");
    write!(w, ">{} {}\nLength={}\n\n", subject_id, ncbi_wrap_title(title), subject_len)?;
    Ok(())
}

/// Write one HSP in pairwise alignment format.
///
/// Outputs the stats block, then the wrapped alignment rows (60 bases/row),
/// then two blank lines (as required by NCBI's between-HSP and end-of-subject
/// blank-line conventions).
///
/// `query_mask`: the query's DUST intervals (0-based, inclusive, sorted).
/// Query bases inside them print in lower case, as NCBI prints filtered
/// query regions.
pub fn write_pairwise_hsp<W: Write>(
    w: &mut W,
    r: &AlignResult,
    query_mask: &[(u32, u32)],
) -> std::io::Result<()> {
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

    writeln!(w, " Score = {}", h.score)?;
    writeln!(w, " Substitutions = {:.2}%, Query Gaps = {:.2}%, DB Gaps = {:.2}%",
        st.perc_sub, st.perc_query_gap, st.perc_db_gap)?;
    writeln!(w, " TransI/TransV = {}/{},  CpG_sites = {}",
        st.transitions, st.transversions, st.cpg_sites)?;
    writeln!(w, " Kimura Div = {:.2}%, Kimura Div CpG Adjusted = {:.2}%",
        st.kdiv, st.cpg_kdiv)?;
    writeln!(w, " Identities = {}/{} ({}%), Gaps = {}/{} ({}%)",
        n_ident, aln_len, pct_ident, n_gaps, aln_len, pct_gaps)?;
    write!(w, " Strand={}\n\n", strand_str)?;

    // ── Alignment rows ────────────────────────────────────────────────────────
    let mut q_iupac = blastna_to_iupac_aligned(&h.q_seq);
    if !query_mask.is_empty() {
        let mut q = h.q_start; // 0-based query position of the next base
        let mut m = query_mask.partition_point(|&(_, e)| e < q);
        for c in q_iupac.iter_mut() {
            if *c == b'-' {
                continue;
            }
            while m < query_mask.len() && query_mask[m].1 < q {
                m += 1;
            }
            if m < query_mask.len() && query_mask[m].0 <= q {
                c.make_ascii_lowercase();
            }
            q += 1;
        }
    }
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
            else if q.eq_ignore_ascii_case(&s) { b'|' }
            else { b' ' }
        }).collect();

        // Prefix width: "Query" (5) + 2 spaces + pos_w digits + 2 spaces = 9 + pos_w.
        let prefix_spaces = 9 + pos_w;

        write!(w, "Query  {:<width$}  ", q_row_start, width = pos_w)?;
        w.write_all(q_row)?;
        writeln!(w, "  {}", q_row_end)?;

        for _ in 0..prefix_spaces { w.write_all(b" ")?; }
        w.write_all(&match_line)?;
        writeln!(w)?;

        write!(w, "Sbjct  {:<width$}  ", s_row_start, width = pos_w)?;
        w.write_all(s_row)?;
        writeln!(w, "  {}", s_row_end)?;

        writeln!(w)?;   // blank line after each chunk

        q_pos += q_nongap;
        match h.strand {
            Strand::Plus  => s_pos += s_nongap,
            Strand::Minus => s_pos = s_pos.saturating_sub(s_nongap),
        }

        i = end;
    }

    // One extra blank line → total 2 blank lines after last chunk.
    writeln!(w)?;

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
    // (each query block already ended with the blank lines NCBI puts here)
    write!(w, "  Database: ")?;
    for line in ncbi_wrap_flat(db_name, 68) {
        writeln!(w, "{}", line)?;
    }
    if let Some(secs) = db_mtime_unix {
        writeln!(w, "    Posted date:  {}", format_ncbi_date(secs))?;
    }
    writeln!(w, "  Number of letters in database: {}", with_commas(total_letters))?;
    writeln!(w, "  Number of sequences in database:  {}", with_commas(n_seqs as u64))?;
    write!(w, "\n\n\n")?;
    writeln!(w, "Matrix: blastn matrix 0 0")?;
    writeln!(w, "Gap Penalties: Existence: {}, Extension: {:.2}",
        gap_open, gap_extend as f64)?;
    Ok(())
}

/// Write all pairwise results for one query.  Results must be pre-sorted by score
/// descending.  Subjects are output in order of first appearance (best-score-first).
///
/// `subject_defline` gives a subject's full defline (None for a 2bit
/// database, which stores only names); `query_mask` is the query's DUST
/// intervals (see `write_pairwise_hsp`).  A query without hits gets NCBI's
/// "***** No hits found *****" block.
pub fn write_pairwise_results<'a, W: Write>(
    w: &mut W,
    results: &[AlignResult],
    query_defline: &str,
    query_len: u32,
    subject_defline: &dyn Fn(&str) -> Option<&'a str>,
    query_mask: &[(u32, u32)],
) -> std::io::Result<()> {
    write_pairwise_query_header(w, query_defline, query_len)?;
    if results.is_empty() {
        write!(w, "\n***** No hits found *****\n\n\n\n\n")?;
        return Ok(());
    }

    // Group by subject, preserving order of first appearance, in one pass.
    let mut order: Vec<&str> = Vec::new();
    let mut groups: std::collections::HashMap<&str, Vec<&AlignResult>> = std::collections::HashMap::new();
    for r in results {
        groups.entry(r.subject_id.as_str())
            .or_insert_with(|| { order.push(r.subject_id.as_str()); Vec::new() })
            .push(r);
    }

    for subj_id in order {
        let subj_results = &groups[subj_id];
        write_pairwise_subject_header(w, subj_id, subject_defline(subj_id), subj_results[0].hsp.s_len)?;
        for r in subj_results {
            write_pairwise_hsp(w, r, query_mask)?;
        }
    }
    // NCBI ends each query's block with two more blank lines.
    write!(w, "\n\n")?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn with_commas_like_ncbi() {
        for (n, want) in [(0, "0"), (17, "17"), (999, "999"), (1000, "1,000"),
                          (11207, "11,207"), (3099734149, "3,099,734,149")] {
            assert_eq!(with_commas(n), want);
        }
    }

    /// NCBI 2.17.1's "Database:" lines for these titles (kFormatLineLength 68).
    #[test]
    fn database_title_wraps_like_ncbi() {
        let path = "databases/repeat-libraries/2026-11-release/curated-families-with-long-names/db.fa";
        assert_eq!(ncbi_wrap_flat(path, 68), [
            "databases/repeat-libraries/2026-11-release/curated-families-with-",
            "long-names/db.fa",
        ]);
        assert_eq!(ncbi_wrap_flat("db.fa", 68), ["db.fa"]);
    }

    /// NCBI 2.17.1's subject header for a -parse_seqids database: the title
    /// breaks after the first whitespace past each 60th character, and a long
    /// word (hyphens included) is never broken.
    #[test]
    fn subject_title_wraps_like_ncbi() {
        let mut out = Vec::new();
        let d = format!("id3 {}", "a very long description with many words ".repeat(4));
        write_pairwise_subject_header(&mut out, "id3", Some(&d), 262).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(),
            ">id3 a very long description with many words a very long description \n\
             with many words a very long description with many words a \n\
             very long description with many words \nLength=262\n\n");
        let mut out = Vec::new();
        let long = "id5 x averyveryveryveryveryveryveryveryveryveryveryveryverylongwordwithnospaces-and-hyphens-here-too";
        write_pairwise_subject_header(&mut out, "id5", Some(long), 554).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), format!(">{long}\nLength=554\n\n"));
        let mut out = Vec::new();
        write_pairwise_subject_header(&mut out, "onlyid", Some("onlyid"), 302).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), ">onlyid \nLength=302\n\n");
    }

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
