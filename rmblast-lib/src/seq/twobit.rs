//! UCSC 2bit file reader.
//!
//! 2bit format summary:
//!   - 16-byte file header: signature(4) version(4) seqCount(4) reserved(4)
//!   - seqCount × (name_size(1) + name(name_size) + offset(4)) index entries
//!   - Per-sequence records at the offsets:
//!       dnaSize(4)  nBlockCount(4)  nBlockStarts[](4)  nBlockSizes[](4)
//!       maskBlockCount(4)  maskBlockStarts[](4)  maskBlockSizes[](4)
//!       reserved(4)  packed_dna(ceil(dnaSize/4) bytes)
//!
//! Bases are packed 4 per byte, MSB-first: T=00 C=01 A=10 G=11.
//! N-blocks overlay the packed sequence with ambiguous bases.
//!
//! We decode to BLASTNA encoding: T→3, C→1, A→0, G→2.
//! N-block positions use NCBI-compatible randomisation: NcbiRandom seeded with
//! the sequence length, one call per N-block position in 5'→3' order.  This
//! matches makeblastdb's CAmbigDataBuilder and lets the scanner form seeds that
//! span N positions exactly as NCBI does.
//!
//! IUB ambiguity limitation: faToTwoBit converts non-ACGT characters
//! (R/Y/K/M/W/S/B/D/H/V) to N, losing the original code.  The .2bit format
//! therefore only supports N-level ambiguity.  All N-block positions are
//! randomised as pure N (any of A/C/G/T).  For full IUB fidelity, supply the
//! database in FASTA format instead (see SubjectDb / DESIGN_NOTES.md).
//!
//! Lowercase (soft-masked) regions in mask blocks → treated as uppercase here;
//! RepeatMasker provides its own masking via -dust / the query mask.

use std::collections::HashMap;
use std::io::{self, Read};
use thiserror::Error;

use crate::encoding::{NCBI_AMBIG_ALLOWED, NUCL_SENTINEL, UCSC2BIT_TO_BLASTNA};

const TWOBIT_MAGIC: u32 = 0x1A412743;
const TWOBIT_MAGIC_SWAPPED: u32 = 0x4327411A;

// ─── NCBI LFG random number generator ────────────────────────────────────────
//
// Mirrors NCBI's CRandom class (util/random_gen.*) used by makeblastdb's
// CAmbigDataBuilder to assign pseudo-random bases at ambiguous positions.
// Seed = sequence length (base_length).  Key detail: GetRand() returns
// x_GetRand32Bits() >> 1, so the base is (GetRand() & 3) = bits [2:1] of
// the raw LFG output — NOT bits [1:0].

const CRANDOM_STATE_SIZE: usize = 33;
const CRANDOM_STATE_OFFSET: usize = 12;

pub(crate) struct NcbiRandom {
    state: [u32; CRANDOM_STATE_SIZE],
    rj: i32,
    rk: i32,
}

impl NcbiRandom {
    pub(crate) fn new(seed: u32) -> Self {
        let mut state = [0u32; CRANDOM_STATE_SIZE];
        state[0] = seed;
        for i in 1..CRANDOM_STATE_SIZE {
            state[i] = 1103515245u32.wrapping_mul(state[i - 1]).wrapping_add(12345);
        }
        let rj = CRANDOM_STATE_OFFSET as i32;
        let rk = (CRANDOM_STATE_SIZE - 1) as i32;
        let mut rng = NcbiRandom { state, rj, rk };
        for _ in 0..(10 * CRANDOM_STATE_SIZE) {
            rng.x_get_rand32();
        }
        rng
    }

    #[inline(always)]
    fn x_get_rand32(&mut self) -> u32 {
        let r = self.state[self.rk as usize]
            .wrapping_add(self.state[self.rj as usize]);
        self.state[self.rk as usize] = r;
        self.rk -= 1;
        self.rj -= 1;
        if self.rk < 0 {
            self.rk = (CRANDOM_STATE_SIZE - 1) as i32;
        } else if self.rj < 0 {
            self.rj = (CRANDOM_STATE_SIZE - 1) as i32;
        }
        r
    }

    #[inline(always)]
    fn get_rand(&mut self) -> u32 {
        self.x_get_rand32() >> 1
    }

    /// Returns a random BLASTNA base (0=A, 1=C, 2=G, 3=T).
    /// ncbi2na and BLASTNA share the same A/C/G/T encoding so no remapping needed.
    #[inline(always)]
    fn next_base(&mut self) -> u8 {
        (self.get_rand() & 3) as u8
    }

    /// Returns a random BLASTNA base constrained to those allowed by `blastna_code`.
    /// Mirrors NCBI's `x_Random` (writedb_convert.cpp): for N (code 14) uses
    /// `GetRand() & 3`; otherwise `pick = GetRand() % bitcount` selects the `pick`-th
    /// allowed base (bits scanned in A,C,G,T order).  Consumes exactly one `GetRand()`
    /// per call, keeping the random stream in sync with makeblastdb.
    #[inline(always)]
    pub(crate) fn next_base_iupac(&mut self, blastna_code: u8) -> u8 {
        let r = self.get_rand();
        if blastna_code == 14 {
            return (r & 3) as u8;
        }
        let entry = &NCBI_AMBIG_ALLOWED[blastna_code as usize];
        let count = entry[0] as u32;
        entry[1 + (r % count) as usize]
    }
}

// ─────────────────────────────────────────────────────────────────────────────

#[derive(Debug, Error)]
pub enum TwoBitError {
    #[error("I/O error: {0}")]
    Io(#[from] io::Error),
    #[error("invalid 2bit magic number {0:#010x}")]
    BadMagic(u32),
    #[error("unsupported 2bit version {0}")]
    BadVersion(u32),
    #[error("sequence '{0}' not found in 2bit file")]
    NotFound(String),
}

/// Metadata for one sequence stored in the 2bit file.
#[derive(Debug, Clone)]
pub struct SeqInfo {
    pub name: String,
    pub file_offset: u64,
    pub dna_size: u32,
}

/// Memory-mapped or buffered view of a 2bit file.
/// Keeps all data in memory for random access; suitable for genome-scale use
/// with memmap2 (see `from_mmap`).
pub struct TwoBitFile {
    data: Vec<u8>,
    swap: bool, // byte-swap all multi-byte integers
    pub sequences: Vec<SeqInfo>,
    name_to_index: HashMap<String, usize>,
}

impl TwoBitFile {
    /// Load a 2bit file fully into memory.
    pub fn open<P: AsRef<std::path::Path>>(path: P) -> Result<Self, TwoBitError> {
        let mut f = std::fs::File::open(path)?;
        let mut data = Vec::new();
        f.read_to_end(&mut data)?;
        Self::from_bytes(data)
    }

    fn from_bytes(data: Vec<u8>) -> Result<Self, TwoBitError> {
        if data.len() < 16 {
            return Err(TwoBitError::Io(io::Error::new(io::ErrorKind::UnexpectedEof, "file too small")));
        }

        let magic = u32::from_le_bytes(data[0..4].try_into().unwrap());
        let swap = match magic {
            TWOBIT_MAGIC => false,
            TWOBIT_MAGIC_SWAPPED => true,
            m => return Err(TwoBitError::BadMagic(m)),
        };

        let r32 = |pos: usize| -> u32 {
            let b = &data[pos..pos + 4];
            let v = u32::from_le_bytes(b.try_into().unwrap());
            if swap { v.swap_bytes() } else { v }
        };

        let version = r32(4);
        if version != 0 {
            return Err(TwoBitError::BadVersion(version));
        }
        let seq_count = r32(8) as usize;
        // reserved = r32(12), ignored

        let mut pos = 16usize;
        let mut sequences = Vec::with_capacity(seq_count);
        let mut name_to_index = HashMap::with_capacity(seq_count);

        for i in 0..seq_count {
            if pos >= data.len() {
                break;
            }
            let name_len = data[pos] as usize;
            pos += 1;
            let name = String::from_utf8_lossy(&data[pos..pos + name_len]).into_owned();
            pos += name_len;
            let offset = r32(pos) as u64;
            pos += 4;
            name_to_index.insert(name.clone(), i);
            sequences.push(SeqInfo { name, file_offset: offset, dna_size: 0 });
        }

        // Fill in dna_size from each sequence record header
        for seq in &mut sequences {
            if (seq.file_offset as usize) + 4 <= data.len() {
                seq.dna_size = r32(seq.file_offset as usize);
            }
        }

        Ok(TwoBitFile { data, swap, sequences, name_to_index })
    }

    fn r32(&self, pos: usize) -> u32 {
        let b = &self.data[pos..pos + 4];
        let v = u32::from_le_bytes(b.try_into().unwrap());
        if self.swap { v.swap_bytes() } else { v }
    }

    /// Decode a subsequence [start, start+len) (0-based, half-open) to BLASTNA.
    /// Returns a Vec with a leading and trailing sentinel (15).
    /// N-block positions use NCBI-compatible random bases (matching makeblastdb).
    pub fn get_sequence_blastna(
        &self,
        name: &str,
        start: u32,
        len: u32,
    ) -> Result<Vec<u8>, TwoBitError> {
        let idx = *self.name_to_index.get(name)
            .ok_or_else(|| TwoBitError::NotFound(name.to_owned()))?;
        let seq = &self.sequences[idx];
        self.decode_sequence(seq, start, len)
    }

    /// Decode the full sequence to BLASTNA (with sentinels).
    /// N-block positions use NCBI-compatible random bases (matching makeblastdb).
    pub fn get_full_sequence_blastna(&self, name: &str) -> Result<Vec<u8>, TwoBitError> {
        let idx = *self.name_to_index.get(name)
            .ok_or_else(|| TwoBitError::NotFound(name.to_owned()))?;
        let seq = &self.sequences[idx];
        self.decode_sequence(seq, 0, seq.dna_size)
    }

    /// Return an ambiguity mask for the named sequence.
    ///
    /// Returns a `Vec<u8>` of length `dna_size` where 0 = unambiguous (A/C/G/T)
    /// and 14 = N (from an N-block).  Returns an empty Vec if the sequence has
    /// no N-blocks.
    ///
    /// Note: the .2bit format does not preserve non-N IUPAC codes (R/Y/K/M/W/S/
    /// B/D/H/V); faToTwoBit converts them all to N.  For full IUB fidelity
    /// supply the database in FASTA format (see SubjectDb).
    pub fn get_n_mask(&self, name: &str) -> Result<Vec<u8>, TwoBitError> {
        let idx = *self.name_to_index.get(name)
            .ok_or_else(|| TwoBitError::NotFound(name.to_owned()))?;
        let seq = &self.sequences[idx];
        let base = seq.file_offset as usize;
        let n_block_count = self.r32(base + 4) as usize;
        if n_block_count == 0 {
            return Ok(Vec::new());
        }
        let n_starts_off = base + 8;
        let n_sizes_off = n_starts_off + n_block_count * 4;
        let dna_size = seq.dna_size as usize;
        let mut mask = vec![0u8; dna_size];
        for bi in 0..n_block_count {
            let n_start = self.r32(n_starts_off + bi * 4) as usize;
            let n_size  = self.r32(n_sizes_off  + bi * 4) as usize;
            let n_end   = (n_start + n_size).min(dna_size);
            for p in n_start..n_end {
                mask[p] = 14;
            }
        }
        Ok(mask)
    }

    fn decode_sequence(&self, seq: &SeqInfo, start: u32, len: u32) -> Result<Vec<u8>, TwoBitError> {
        let base = seq.file_offset as usize;
        // Layout: dnaSize(4) nBlockCount(4) nBlockStarts[nBlockCount](4) nBlockSizes[nBlockCount](4)
        //         maskBlockCount(4) maskBlockStarts[](4) maskBlockSizes[](4) reserved(4) dna_data
        let n_block_count = self.r32(base + 4) as usize;
        let n_starts_off = base + 8;
        let n_sizes_off = n_starts_off + n_block_count * 4;
        let mask_count_off = n_sizes_off + n_block_count * 4;
        let mask_block_count = self.r32(mask_count_off) as usize;
        let mask_starts_off = mask_count_off + 4;
        let mask_sizes_off = mask_starts_off + mask_block_count * 4;
        let dna_off = mask_sizes_off + mask_block_count * 4 + 4; // +4 for reserved

        let end = start + len;
        let dna_size = seq.dna_size;
        let actual_end = end.min(dna_size);
        let actual_len = if actual_end > start { actual_end - start } else { 0 };

        let mut out = vec![NUCL_SENTINEL; actual_len as usize + 2];

        // Decode packed DNA into BLASTNA
        for i in 0..actual_len {
            let pos = start + i;
            let byte_idx = (pos / 4) as usize;
            let bit_shift = 6 - ((pos % 4) * 2);
            let two_bits = ((self.data[dna_off + byte_idx] >> bit_shift) & 0x3) as usize;
            out[(i + 1) as usize] = UCSC2BIT_TO_BLASTNA[two_bits];
        }

        // Replace N-block positions with NCBI-compatible pseudo-random bases.
        // .2bit files only carry N-level ambiguity (faToTwoBit converts all
        // non-ACGT characters to N), so we randomise uniformly over all four
        // bases, matching makeblastdb's GetRand()&3 for N positions.
        if n_block_count > 0 {
            let mut rng = NcbiRandom::new(dna_size);
            for bi in 0..n_block_count {
                let n_start = self.r32(n_starts_off + bi * 4);
                let n_size  = self.r32(n_sizes_off  + bi * 4);
                let n_end   = n_start + n_size;
                for pos in n_start..n_end.min(dna_size) {
                    let rand_base = rng.next_base();
                    if pos >= start && pos < actual_end {
                        out[(pos - start + 1) as usize] = rand_base;
                    }
                }
            }
        }

        Ok(out)
    }

    /// Iterate over all sequence names.
    pub fn seq_names(&self) -> impl Iterator<Item = &str> {
        self.sequences.iter().map(|s| s.name.as_str())
    }
}

/// Iterator over all sequences in a 2bit file decoded to BLASTNA.
pub struct TwoBitSeqIter<'a> {
    file: &'a TwoBitFile,
    idx: usize,
}

impl<'a> TwoBitSeqIter<'a> {
    pub fn new(file: &'a TwoBitFile) -> Self {
        TwoBitSeqIter { file, idx: 0 }
    }
}

impl<'a> Iterator for TwoBitSeqIter<'a> {
    type Item = Result<(String, Vec<u8>), TwoBitError>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.idx >= self.file.sequences.len() {
            return None;
        }
        let seq = &self.file.sequences[self.idx];
        self.idx += 1;
        let name = seq.name.clone();
        let dna_size = seq.dna_size;
        Some(self.file.decode_sequence(seq, 0, dna_size).map(|v| (name, v)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a minimal valid 2bit byte stream for "ACGT" (4 bases, no N/mask blocks).
    fn make_test_twobit() -> Vec<u8> {
        let mut v: Vec<u8> = Vec::new();

        // File header
        v.extend_from_slice(&TWOBIT_MAGIC.to_le_bytes()); // magic
        v.extend_from_slice(&0u32.to_le_bytes());           // version
        v.extend_from_slice(&1u32.to_le_bytes());           // seqCount = 1
        v.extend_from_slice(&0u32.to_le_bytes());           // reserved

        // Index entry
        v.push(4u8); // name length
        v.extend_from_slice(b"seq1"); // name
        // offset to sequence record = 16 (header) + 1 + 4 + 4 = 25
        let seq_offset: u32 = 25;
        v.extend_from_slice(&seq_offset.to_le_bytes());

        // Sequence record for "ACGT"
        // Packed: A=10, C=01, G=11, T=00 → 0b10011100 = 0x9C
        v.extend_from_slice(&4u32.to_le_bytes()); // dnaSize = 4
        v.extend_from_slice(&0u32.to_le_bytes()); // nBlockCount = 0
        v.extend_from_slice(&0u32.to_le_bytes()); // maskBlockCount = 0
        v.extend_from_slice(&0u32.to_le_bytes()); // reserved
        v.push(0x9C); // packed ACGT

        v
    }

    #[test]
    fn test_decode_acgt() {
        let data = make_test_twobit();
        let tb = TwoBitFile::from_bytes(data).unwrap();
        let seq = tb.get_full_sequence_blastna("seq1").unwrap();
        // sentinels + ACGT in BLASTNA: A=0,C=1,G=2,T=3
        assert_eq!(seq, vec![15, 0, 1, 2, 3, 15]);
    }

    #[test]
    fn test_ncbi_random() {
        // Verify CRandom(23) matches actual makeblastdb output for ACGTACGTNACGTNACGTNACGT
        // (Ns at positions 8, 13, 18). NSQ confirmed values: T(3), A(0), G(2).
        let mut rng = NcbiRandom::new(23);
        assert_eq!(rng.next_base(), 3); // T at pos 8
        assert_eq!(rng.next_base(), 0); // A at pos 13
        assert_eq!(rng.next_base(), 2); // G at pos 18
    }

    #[test]
    fn test_ncbi_random_428() {
        // MSTC#LTR/ERVL-MaLR has length=428, single N at position 343
        // NSQ confirmed: position 343 = C(1)
        let mut rng = NcbiRandom::new(428);
        assert_eq!(rng.next_base(), 1); // C at pos 343
    }
}
