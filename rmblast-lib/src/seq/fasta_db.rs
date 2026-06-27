//! FASTA subject-database reader with full IUB/IUPAC support.
//!
//! Loads all sequences into memory preserving their original IUPAC codes
//! (R/Y/K/M/W/S/B/D/H/V/N).  Ambiguous positions are resolved to specific
//! bases at decode time using the same NCBI-compatible constrained-random
//! algorithm as makeblastdb: NcbiRandom seeded with the sequence length,
//! one GetRand() call per ambiguous position in 5'→3' order.
//!
//! This matches the seeding behaviour produced by makeblastdb exactly,
//! enabling full IUB-aware comparison with NCBI rmblastn.

use std::collections::HashMap;
use std::io::BufReader;
use std::path::Path;

use thiserror::Error;

use crate::seq::fasta::FastaReader;
use crate::seq::twobit::{NcbiRandom, SeqInfo};

#[derive(Debug, Error)]
pub enum FastaDbError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("FASTA parse error: {0}")]
    Fasta(String),
    #[error("sequence '{0}' not found in FASTA database")]
    NotFound(String),
}

/// An in-memory FASTA subject database with full IUB/IUPAC support.
///
/// Sequences are stored with their original BLASTNA codes (including 4–14 for
/// IUPAC ambiguities).  `get_full_sequence_blastna` resolves ambiguous positions
/// to specific bases using NcbiRandom, matching makeblastdb behaviour.
pub struct FastaDb {
    pub sequences: Vec<SeqInfo>,
    name_to_index: HashMap<String, usize>,
    /// Raw BLASTNA sequences (sentinel | bases | sentinel).
    /// IUPAC codes 4–14 are preserved; they are resolved at decode time.
    raw_seqs: Vec<Vec<u8>>,
}

impl FastaDb {
    /// Load all sequences from a FASTA file.
    pub fn open<P: AsRef<Path>>(path: P) -> Result<Self, FastaDbError> {
        let f = std::fs::File::open(path.as_ref())?;
        let mut reader = FastaReader::new(BufReader::new(f));
        let mut sequences: Vec<SeqInfo> = Vec::new();
        let mut name_to_index: HashMap<String, usize> = HashMap::new();
        let mut raw_seqs: Vec<Vec<u8>> = Vec::new();

        while let Some(rec) = reader
            .next_record()
            .map_err(|e| FastaDbError::Fasta(e.to_string()))?
        {
            let dna_size = rec.len() as u32;
            let idx = sequences.len();
            name_to_index.insert(rec.id.clone(), idx);
            sequences.push(SeqInfo { name: rec.id, file_offset: 0, dna_size });
            raw_seqs.push(rec.seq);
        }

        Ok(FastaDb { sequences, name_to_index, raw_seqs })
    }

    /// Decode a sequence to BLASTNA, resolving every IUB ambiguity (codes 4–14)
    /// to a specific base using NCBI's constrained-random algorithm.
    ///
    /// The RNG is seeded with `dna_size` and advanced once per ambiguous position
    /// in 5'→3' order, matching makeblastdb's `CAmbigDataBuilder::x_Random`.
    pub fn get_full_sequence_blastna(&self, name: &str) -> Result<Vec<u8>, FastaDbError> {
        let &idx = self
            .name_to_index
            .get(name)
            .ok_or_else(|| FastaDbError::NotFound(name.to_owned()))?;
        let dna_size = self.sequences[idx].dna_size;
        let mut out = self.raw_seqs[idx].clone();
        let last = out.len().saturating_sub(1);
        let mut rng = NcbiRandom::new(dna_size);
        // Iterate only the bases (skip sentinels at index 0 and last).
        for slot in &mut out[1..last] {
            let code = *slot;
            if code >= 4 && code <= 14 {
                *slot = rng.next_base_iupac(code);
            }
        }
        Ok(out)
    }

    /// Return the ambiguity mask for a named sequence.
    ///
    /// Returns a `Vec<u8>` of length `dna_size` where 0 = unambiguous (A/C/G/T)
    /// and 4–14 = original BLASTNA code at that position.  Returns an empty Vec
    /// if the sequence contains no ambiguous positions.
    pub fn get_n_mask(&self, name: &str) -> Vec<u8> {
        let Some(&idx) = self.name_to_index.get(name) else {
            return Vec::new();
        };
        let raw = &self.raw_seqs[idx];
        let bases = &raw[1..raw.len() - 1];
        if !bases.iter().any(|&c| c >= 4 && c <= 14) {
            return Vec::new();
        }
        let mut mask = vec![0u8; bases.len()];
        for (i, &c) in bases.iter().enumerate() {
            if c >= 4 && c <= 14 {
                mask[i] = c;
            }
        }
        mask
    }
}
