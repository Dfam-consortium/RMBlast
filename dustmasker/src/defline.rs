//! The defline NCBI dustmasker prints for a FASTA record.
//!
//! Without `-parse_seqids`, NCBI's dustmasker prints
//! `CDeflineGenerator().GenerateDefline()` (objtools/seqmasks_io/mask_writer.cpp),
//! and the FASTA reader stores the whole defline as the title, so the id goes
//! through the same clean-up as the description.  This module reproduces that
//! clean-up for the plain-FASTA case:
//!
//! 1. leading whitespace is dropped and the line is cut at its first tab
//!    (observed from dustmasker 2.17.0; done by the FASTA reader)
//! 2. trailing `.,;~ ` are stripped (create_defline.cpp, after x_SetFlags)
//! 3. the TPA/TSA/MAG/... prefixes are removed (`s_tpaPrefixList`)
//! 4. HTML entities and character references are decoded (`NStr::HtmlDecode`)
//! 5. trailing `,;~ ` are stripped again
//! 6. `x_CleanAndCompress`: runs of spaces collapse, `" ,"`/`" ;"` swap,
//!    spaces next to brackets go, and so on
//!
//! Checked against NCBI dustmasker 2.17.0 on randomly generated deflines
//! (see `run_dustmasker_compare.sh`).

use crate::html_entities::HTML_ENTITIES;

/// The defline (without '>') that NCBI dustmasker prints for `raw`, as
/// bytes: NCBI can cut a final non-ASCII character short, and so do we.
pub fn ncbi_defline(raw: &str) -> Vec<u8> {
    // 1. FASTA reader: leading whitespace, cut at the first tab
    let s = raw.trim_start_matches(|c: char| c.is_ascii_whitespace());
    let s = s.split('\t').next().unwrap_or("");

    // 2. strip trailing periods, commas, semicolons, tildes and spaces
    // (only ASCII is removed, so the rest is still UTF-8)
    let mut title = String::from_utf8(strip_trailing(s.as_bytes(), b".,;~ ")).unwrap_or_default();

    // 3. TPA and similar prefixes (case-insensitive), then leading spaces
    for prefix in TPA_PREFIXES {
        if title.len() >= prefix.len()
            && title.as_bytes()[..prefix.len()].eq_ignore_ascii_case(prefix.as_bytes())
        {
            title.drain(..prefix.len());
            title.drain(..title.len() - title.trim_start_matches(' ').len());
        }
    }
    let title = title.trim_start_matches(' ');

    // 4. HTML entities, 5. trailing commas, semicolons, tildes, spaces
    let decoded = html_decode(title);
    let decoded = strip_trailing(&decoded, b",;~ ");

    // 6.
    clean_and_compress(&decoded)
}

/// `s.erase(s.find_last_not_of(set) + 1)` guarded by `pos != NPOS`, as NCBI
/// writes it: a string made only of `set` characters is left whole (so a
/// defline of "..." prints as "...").
fn strip_trailing(s: &[u8], set: &[u8]) -> Vec<u8> {
    match s.iter().rposition(|c| !set.contains(c)) {
        Some(p) => s[..=p].to_vec(),
        None => s.to_vec(),
    }
}

/// `s_tpaPrefixList` (create_defline.cpp), in NCBI's order.
const TPA_PREFIXES: &[&str] = &[
    "MAG ", "MAG:", "MULTISPECIES:", "TLS:", "TPA:", "TPA_exp:", "TPA_inf:",
    "TPA_reasm:", "TPA_asm:", "TSA:", "UNVERIFIED_ORG:", "UNVERIFIED_ASMBLY:",
    "UNVERIFIED_CONTAM:", "UNVERIFIED:",
];

/// Port of `NStr::HtmlDecode` (corelib/ncbistr.cpp) for UTF-8 input, quirks
/// included: an entity name is at most 16 characters, a '&' or '#' inside it
/// ends the attempt, and a name starting with x + hex digit is looked up
/// without its first letter (so `&xi;` is not decoded, as in NCBI).
fn html_decode(s: &str) -> Vec<u8> {
    let b = s.as_bytes();
    let n = b.len();
    let mut out: Vec<u8> = Vec::with_capacity(n);
    let mut i = 0;
    while i < n {
        let ch = b[i];
        i += 1;
        if i != n && ch == b'&' {
            let mut itmp = i;
            let mut end_of_entity = i;
            let mut ent = b[itmp].is_ascii_alphabetic();
            let mut dec = false;
            if !ent && b[itmp] == b'#' {
                itmp += 1;
                dec = itmp != n && b[itmp].is_ascii_digit();
            }
            let mut hex = false;
            if !dec && itmp != n && (b[itmp] == b'x' || b[itmp] == b'X') {
                itmp += 1;
                hex = itmp != n && b[itmp].is_ascii_hexdigit();
            }
            let start_of_entity = itmp;
            if itmp != n && (ent || dec || hex) {
                let mut len = 0;
                while len < 16 && itmp != n {
                    let c = b[itmp];
                    if c == b'&' || c == b'#' {
                        break;
                    }
                    if c == b';' {
                        end_of_entity = itmp;
                        break;
                    }
                    ent = ent && c.is_ascii_alphanumeric();
                    dec = dec && c.is_ascii_digit();
                    hex = hex && c.is_ascii_hexdigit();
                    len += 1;
                    itmp += 1;
                }
                if end_of_entity != i && (ent || dec || hex) {
                    let body = &s[start_of_entity..end_of_entity];
                    let uch: Option<u32> = if ent {
                        HTML_ENTITIES.iter().find(|(_, name)| *name == body).map(|(u, _)| *u)
                    } else {
                        let radix = if dec { 10 } else { 16 };
                        let mut v: u32 = 0;
                        for c in body.chars() {
                            v = v.wrapping_mul(radix).wrapping_add(c.to_digit(radix).unwrap_or(0));
                        }
                        Some(v)
                    };
                    if let Some(u) = uch {
                        append_utf8(&mut out, u);
                        i = end_of_entity + 1;
                        continue;
                    }
                }
            }
        }
        out.push(ch);
    }
    out
}

/// `CUtf8::x_AppendChar` (corelib/ncbistr.cpp): the UTF-8 pattern applied to
/// any 32-bit value with no range check, so surrogates and values above
/// U+10FFFF come out as NCBI writes them (`&#x110000;` -> F4 90 80 80).
fn append_utf8(out: &mut Vec<u8>, ch: u32) {
    if ch < 0x80 {
        out.push(ch as u8);
    } else if ch < 0x800 {
        out.extend([((ch >> 6) | 0xC0) as u8, ((ch & 0x3F) | 0x80) as u8]);
    } else if ch < 0x10000 {
        out.extend([((ch >> 12) | 0xE0) as u8, (((ch >> 6) & 0x3F) | 0x80) as u8,
                    ((ch & 0x3F) | 0x80) as u8]);
    } else {
        out.extend([((ch >> 18) | 0xF0) as u8, (((ch >> 12) & 0x3F) | 0x80) as u8,
                    (((ch >> 6) & 0x3F) | 0x80) as u8, ((ch & 0x3F) | 0x80) as u8]);
    }
}

/// Port of `x_CleanAndCompress` (create_defline.cpp, itself copied from
/// `CleanAndCompress` in objtools/format/utils.cpp) for nucleotides.  Bytes,
/// not a String: like NCBI's, the result can end in a truncated UTF-8
/// sequence.
fn clean_and_compress(input: &[u8]) -> Vec<u8> {
    let b = input;
    // skip spaces at both ends
    let mut start = 0;
    let mut left = b.len();
    while left > 0 && b[start] == b' ' {
        start += 1;
        left -= 1;
    }
    while left > 0 && b[start + left - 1] == b' ' {
        left -= 1;
    }
    if left < 1 {
        return Vec::new();
    }
    let b = &b[start..start + left];
    let mut out: Vec<u8> = Vec::with_capacity(left);
    let pair = |a: u8, c: u8| ((a as u16) << 8) | c as u16;

    let mut pos = 0;
    let mut curr = b[pos];
    pos += 1;
    left -= 1;
    let mut two_chars: u16 = curr as u16;
    while left > 0 {
        let mut next = b[pos];
        pos += 1;
        two_chars = (two_chars << 8) | next as u16;
        if two_chars == pair(b',', b',') {
            // replace double commas with comma + space
            out.push(curr);
            next = b' ';
        } else if two_chars == pair(b' ', b' ') {
            // skip runs of spaces (keep the last)
        } else if two_chars == pair(b'(', b' ') {
            // skip a space after an opening bracket
            next = curr;
            two_chars = curr as u16;
        } else if two_chars == pair(b' ', b')') {
            // skip a space before a closing bracket
        } else if two_chars == pair(b' ', b',') || two_chars == pair(b' ', b';') {
            // swap the space and the comma/semicolon
            out.push(next);
            next = curr;
            two_chars = curr as u16;
        } else if two_chars == pair(b',', b' ') || two_chars == pair(b';', b' ') {
            out.push(curr);
            out.push(b' ');
            // the pair's first byte, not `curr`: after ",," the C has already
            // replaced `curr` with a space
            let skip = (two_chars >> 8) as u8;
            // the C reads one past the end here and relies on the terminating
            // NUL; reading 0 at the end does the same
            while next == b' ' || next == skip {
                next = b.get(pos).copied().unwrap_or(0);
                pos += 1;
                left = left.wrapping_sub(1);
            }
            two_chars = next as u16;
        } else {
            out.push(curr);
        }
        curr = next;
        left = left.wrapping_sub(1);
        if left > b.len() {
            break; // ran past the end in the loop above (C: left underflows to 0 checks)
        }
    }
    // C: `if (curr > 0 && curr != ' ')` on a signed char, so a final byte of
    // 0x80 or above is dropped: NCBI cuts the last UTF-8 sequence short.
    if curr > 0 && curr < 0x80 && curr != b' ' {
        out.push(curr);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// NCBI dustmasker 2.17.0 output for these deflines.
    #[test]
    fn matches_ncbi_examples() {
        for (raw, want) in [
            ("s1  two  spaces\ttab trailing.  ", "s1 two spaces"),
            ("s2 ends with period.", "s2 ends with period"),
            (" lead", "lead"),
            ("s4\tTabbed desc", "s4"),
            ("s6 two periods..", "s6 two periods"),
            ("s7 a. b.", "s7 a. b"),
            ("s12 ;semi", "s12; semi"),
            ("s13 x\ty\tz", "s13 x"),
            ("\ts17 tab lead", "s17 tab lead"),
            ("s18 ...", "s18"),
            ("s20 ends,comma,", "s20 ends,comma"),
            ("s21 a  b   c", "s21 a b c"),
            ("gi|123|ref|NC_1.1| desc", "gi|123|ref|NC_1.1| desc"),
            ("s27 (paren)", "s27 (paren)"),
            ("id292  &#x42;&#65;,, ,&#9;\t\t", "id292 BA,  &#9"),
            ("...", "..."),
            ("~", "~"),
            ("TPA: ,", ""),
            ("MAG ", "MAG"),
        ] {
            assert_eq!(ncbi_defline(raw), want.as_bytes(), "{raw:?}");
        }
    }

    /// NCBI drops a final byte >= 0x80 (signed char), so "café" prints as
    /// "caf\xC3".
    #[test]
    fn final_non_ascii_byte_is_dropped_like_ncbi() {
        assert_eq!(ncbi_defline("s caf\u{e9}"), b"s caf\xC3");
        assert_eq!(ncbi_defline("s caf\u{e9} x"), "s caf\u{e9} x".as_bytes());
    }
}
