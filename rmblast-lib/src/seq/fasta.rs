//! FASTA reader that produces BLASTNA-encoded query sequences.
//!
//! Each record carries the full BLASTNA encoding with leading/trailing
//! sentinels (15), ready for the lookup table and alignment engine.
//!
//! Supports multi-line FASTA; reads from any BufRead.

use std::io::{self, BufRead};
use thiserror::Error;

use crate::encoding::encode_iupac;

#[derive(Debug, Error)]
pub enum FastaError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    #[error("unexpected data before first '>' header")]
    NoHeader,
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
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(reader: R) -> Self {
        FastaReader { reader, peeked_header: None, done: false }
    }

    /// Read the next FASTA record, or `None` at EOF.
    pub fn next_record(&mut self) -> Result<Option<FastaRecord>, FastaError> {
        if self.done {
            return Ok(None);
        }

        let mut header: Option<String> = self.peeked_header.take();
        let mut raw_seq: Vec<u8> = Vec::new();

        let mut line = String::new();
        loop {
            line.clear();
            let n = self.reader.read_line(&mut line)?;
            if n == 0 {
                self.done = true;
                break;
            }
            let trimmed = line.trim_end();
            if trimmed.is_empty() {
                continue;
            }
            if trimmed.starts_with('>') {
                let defline = trimmed[1..].to_owned();
                if header.is_none() {
                    header = Some(defline);
                } else {
                    self.peeked_header = Some(defline);
                    break;
                }
            } else {
                if header.is_none() {
                    // Data before any header — skip or error
                    continue;
                }
                raw_seq.extend_from_slice(trimmed.as_bytes());
            }
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

    #[test]
    fn test_empty_fasta() {
        let data = b"";
        let mut reader = FastaReader::new(Cursor::new(data));
        assert!(reader.next_record().unwrap().is_none());
    }
}
