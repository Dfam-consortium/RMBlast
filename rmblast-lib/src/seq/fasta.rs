//! FASTA reader that produces BLASTNA-encoded query sequences.
//!
//! Each record carries the full BLASTNA encoding with leading/trailing
//! sentinels (15), ready for the lookup table and alignment engine.
//!
//! Supports multi-line FASTA; reads from any BufRead.  Accepted as NCBI
//! accepts them: CRLF line ends, blank lines, whitespace at either end of a
//! line, lower case, IUPAC codes (plus U as T, X and '-' as N), empty records
//! and repeated ids.  Refused with the record and line number, where NCBI
//! guesses or is inconsistent: sequence before the first '>', a '>' with no
//! id, ';' comment lines, a defline that is not UTF-8, whitespace inside a
//! sequence line (NCBI counts it as residues in a query but drops it in a
//! database), and any other character in a sequence line.

use std::io::{self, BufRead};
use thiserror::Error;

use crate::encoding::encode_iupac;

#[derive(Debug, Error)]
pub enum FastaError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    /// A malformed line; the message says what and where.
    #[error("FASTA line {line}: {msg}")]
    Format { line: usize, msg: String },
}

/// A single FASTA record in BLASTNA encoding.
#[derive(Debug, Clone)]
pub struct FastaRecord {
    /// Sequence identifier (everything after '>' up to first whitespace).
    pub id: String,
    /// Full defline (everything after '>').
    pub defline: String,
    /// BLASTNA-encoded sequence with leading/trailing sentinel (15).
    /// Sequence bases start at index 1.
    pub seq: Vec<u8>,
}

impl FastaRecord {
    /// Length of the sequence (excluding sentinels).
    #[inline]
    pub fn len(&self) -> usize {
        self.seq.len().saturating_sub(2)
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Slice of BLASTNA bases (no sentinels).
    #[inline]
    pub fn bases(&self) -> &[u8] {
        if self.seq.len() < 2 {
            return &[];
        }
        &self.seq[1..self.seq.len() - 1]
    }
}

/// Streaming FASTA reader.  Call `next_record()` to retrieve records one by one.
pub struct FastaReader<R: BufRead> {
    reader: R,
    peeked_header: Option<String>,
    done: bool,
    /// Lines read so far (1-based number of the last line read).
    line_no: usize,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(reader: R) -> Self {
        FastaReader { reader, peeked_header: None, done: false, line_no: 0 }
    }

    fn format_err(&self, msg: String) -> FastaError {
        FastaError::Format { line: self.line_no, msg }
    }

    /// Read the next FASTA record, or `None` at EOF.
    pub fn next_record(&mut self) -> Result<Option<FastaRecord>, FastaError> {
        if self.done {
            return Ok(None);
        }

        let mut header: Option<String> = self.peeked_header.take();
        let mut raw_seq: Vec<u8> = Vec::new();

        let mut buf: Vec<u8> = Vec::new();
        loop {
            buf.clear();
            let n = self.reader.read_until(b'\n', &mut buf)?;
            if n == 0 {
                self.done = true;
                break;
            }
            self.line_no += 1;
            if let Some(rest) = buf.strip_prefix(b">") {
                let defline = std::str::from_utf8(rest)
                    .map_err(|_| self.format_err("the defline is not UTF-8 text".into()))?
                    .trim_end()
                    .to_owned();
                if defline.split_whitespace().next().is_none() {
                    return Err(self.format_err("'>' with no sequence id".into()));
                }
                if header.is_none() {
                    header = Some(defline);
                } else {
                    self.peeked_header = Some(defline);
                    break;
                }
                continue;
            }
            let line = buf.trim_ascii();
            if line.is_empty() {
                continue;
            }
            if line[0] == b';' {
                return Err(self.format_err("';' comment lines are not supported".into()));
            }
            if line[0] == b'>' {
                return Err(self.format_err("a header's '>' must be the first character of the line".into()));
            }
            let Some(defline) = &header else {
                return Err(self.format_err("sequence data before the first '>' header".into()));
            };
            if let Some(&bad) = line.iter().find(|&&c| !is_sequence_char(c)) {
                let id = defline.split_whitespace().next().unwrap_or("");
                let what = if bad.is_ascii_whitespace() {
                    "whitespace inside a sequence line".to_string()
                } else if bad.is_ascii_graphic() {
                    format!("'{}' is not a nucleotide code", bad as char)
                } else {
                    format!("byte 0x{bad:02X} is not a nucleotide code")
                };
                return Err(self.format_err(format!("{what} (in record '{id}')")));
            }
            raw_seq.extend_from_slice(line);
        }

        match header {
            None => Ok(None),
            Some(defline) => {
                let id = defline.split_whitespace().next().unwrap_or("").to_owned();
                let seq = encode_iupac(&raw_seq);
                Ok(Some(FastaRecord { id, defline, seq }))
            }
        }
    }
}

/// Letters that `encode_iupac` gives a meaning: the IUPAC codes, U (as T),
/// X (as N), and '-' (as N).
fn is_sequence_char(c: u8) -> bool {
    matches!(c.to_ascii_uppercase(),
        b'A' | b'C' | b'G' | b'T' | b'U' | b'R' | b'Y' | b'M' | b'K' | b'W' | b'S'
        | b'B' | b'D' | b'H' | b'V' | b'N' | b'X' | b'-')
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = Result<FastaRecord, FastaError>;
    fn next(&mut self) -> Option<Self::Item> {
        match self.next_record() {
            Ok(Some(r)) => Some(Ok(r)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_simple_fasta() {
        let data = b">seq1 test\nACGT\nACGT\n>seq2\nTTTT\n";
        let mut reader = FastaReader::new(Cursor::new(data));
        let r1 = reader.next_record().unwrap().unwrap();
        assert_eq!(r1.id, "seq1");
        assert_eq!(r1.len(), 8);
        assert_eq!(r1.bases(), &[0, 1, 2, 3, 0, 1, 2, 3]); // ACGTACGT

        let r2 = reader.next_record().unwrap().unwrap();
        assert_eq!(r2.id, "seq2");
        assert_eq!(r2.bases(), &[3, 3, 3, 3]); // TTTT

        assert!(reader.next_record().unwrap().is_none());
    }

    fn read_all(data: &[u8]) -> Result<Vec<FastaRecord>, FastaError> {
        FastaReader::new(Cursor::new(data.to_vec())).collect()
    }

    /// The variations NCBI accepts, with the coordinates NCBI gives them.
    #[test]
    fn accepted_variations() {
        let want = read_all(b">a desc\nACGT\nacgt\n>b\n>a\nNNRY\n").unwrap();
        for data in [
            &b">a desc\r\nACGT\r\nacgt\r\n>b\r\n>a\r\nNNRY\r\n"[..], // CRLF
            &b">a desc\n\nACGT  \n   acgt\n\n>b\n>a\n\tNNRY"[..],     // blank lines, end whitespace, no final newline
        ] {
            let got = read_all(data).unwrap();
            assert_eq!(got.len(), 3);
            for (g, w) in got.iter().zip(&want) {
                assert_eq!((&g.id, &g.defline, &g.seq), (&w.id, &w.defline, &w.seq));
            }
        }
        assert!(want[1].is_empty()); // an empty record is a record
        assert_eq!(want[0].bases(), &[0, 1, 2, 3, 0, 1, 2, 3]);
        assert_eq!(read_all(b">x\nAU-X\n").unwrap()[0].bases(), &[0, 3, 14, 14]);
    }

    /// Each malformed input is an error naming the line; 3.0.12 skipped,
    /// guessed or aborted with "stream did not contain valid UTF-8".
    #[test]
    fn malformed_input_is_rejected() {
        let cases: &[(&[u8], &str)] = &[
            (b"ACGT\n>a\nACGT\n", "line 1: sequence data before the first '>' header"),
            (b">a\nACGT\n>\nACGT\n", "line 3: '>' with no sequence id"),
            (b">a\nACGT\n> \nACGT\n", "line 3: '>' with no sequence id"),
            (b">a\n;comment\nACGT\n", "line 2: ';' comment lines"),
            (b">a caf\xe9\nACGT\n", "line 1: the defline is not UTF-8"),
            (b">a\nAC GT\n", "line 2: whitespace inside a sequence line (in record 'a')"),
            (b">a\nACGT*\n", "line 2: '*' is not a nucleotide code"),
            (b">a\nACGT\n>b\nAC1GT\n", "line 4: '1' is not a nucleotide code (in record 'b')"),
            (b">a\nAC\x00GT\n", "byte 0x00"),
            (b"  >a\nACGT\n", "line 1: a header's '>' must be the first character"),
        ];
        for (data, needle) in cases {
            match read_all(data) {
                Err(e) => assert!(e.to_string().contains(needle), "{:?}: {e}", String::from_utf8_lossy(data)),
                Ok(_) => panic!("{:?} accepted", String::from_utf8_lossy(data)),
            }
        }
    }

    #[test]
    fn test_empty_fasta() {
        let data = b"";
        let mut reader = FastaReader::new(Cursor::new(data));
        assert!(reader.next_record().unwrap().is_none());
    }
}
