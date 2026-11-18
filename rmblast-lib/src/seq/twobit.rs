//! UCSC 2bit file reader.
//!
//! 2bit format summary:
//!   - 16-byte file header: signature(4) version(4) seqCount(4) reserved(4)
//!   - seqCount × (name_size(1) + name(name_size) + offset(4)) index entries
//!   - Per-sequence records at the offsets:
//!     dnaSize(4)  nBlockCount(4)  nBlockStarts[](4)  nBlockSizes[](4)
//!     maskBlockCount(4)  maskBlockStarts[](4)  maskBlockSizes[](4)
//!     reserved(4)  packed_dna(ceil(dnaSize/4) bytes)
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
//! database in FASTA format instead (see SubjectDb, and "Known differences
//! from NCBI" in README.md).
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

/// An [`NcbiRandom`] seeded on first use.  Seeding costs about 6,000
/// instructions, and most sequences have no ambiguous base to draw for; the
/// draws are the same as with an eagerly seeded generator.
pub(crate) struct LazyRandom {
    seed: u32,
    rng: Option<NcbiRandom>,
}

impl LazyRandom {
    pub(crate) fn new(seed: u32) -> Self {
        LazyRandom { seed, rng: None }
    }

    #[inline]
    pub(crate) fn get(&mut self) -> &mut NcbiRandom {
        let seed = self.seed;
        self.rng.get_or_insert_with(|| NcbiRandom::new(seed))
    }
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
    pub(crate) fn get_rand(&mut self) -> u32 {
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
    /// A truncated or inconsistent file.  The message says what is wrong and
    /// where (byte offsets).
    #[error("corrupt 2bit file: {0}")]
    Corrupt(String),
    /// A sequence longer than the engine supports (`seq::MAX_SEQ_LEN`).
    #[error("{0}")]
    TooLong(String),
    /// Two index entries share a name; see `FastaDbError::DuplicateId`.
    #[error("duplicate sequence name '{name}' (index entries {first} and {second})")]
    DuplicateName { name: String, first: usize, second: usize },
}

/// Metadata for one sequence stored in the 2bit file.
#[derive(Debug, Clone)]
pub struct SeqInfo {
    pub name: String,
    pub file_offset: u64,
    pub dna_size: u32,
}

/// Where one record's parts are in the file, checked against its length at
/// load so decoding never reads out of bounds.
#[derive(Debug, Clone, Copy)]
struct RecordLayout {
    n_block_count: usize,
    n_starts_off: usize,
    n_sizes_off: usize,
    dna_off: usize,
}

/// A 2bit file held in memory for random access.
pub struct TwoBitFile {
    data: Vec<u8>,
    swap: bool, // byte-swap all multi-byte integers
    pub sequences: Vec<SeqInfo>,
    layouts: Vec<RecordLayout>,
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

    /// Parse and check the whole file: every index entry and every record
    /// header, block table and packed-base run must lie inside the file.  A
    /// corrupt file is an error here rather than a panic (or a silently short
    /// database) later.
    fn from_bytes(data: Vec<u8>) -> Result<Self, TwoBitError> {
        let corrupt = |msg: String| TwoBitError::Corrupt(msg);
        if data.len() < 16 {
            return Err(corrupt(format!("{} bytes, shorter than the 16-byte header", data.len())));
        }

        let magic = u32::from_le_bytes(data[0..4].try_into().unwrap());
        let swap = match magic {
            TWOBIT_MAGIC => false,
            TWOBIT_MAGIC_SWAPPED => true,
            m => return Err(TwoBitError::BadMagic(m)),
        };

        // A u32 at `pos`, or an error naming `what` if it runs past the end.
        let r32 = |pos: usize, what: &str| -> Result<u32, TwoBitError> {
            let b = pos.checked_add(4).and_then(|end| data.get(pos..end)).ok_or_else(|| {
                corrupt(format!("{what} at byte {pos} is past the end of the file ({} bytes)", data.len()))
            })?;
            let v = u32::from_le_bytes(b.try_into().unwrap());
            Ok(if swap { v.swap_bytes() } else { v })
        };
        // `count` u32s starting at `pos`: returns the offset just past them.
        let skip_u32s = |pos: usize, count: usize, what: &str| -> Result<usize, TwoBitError> {
            count.checked_mul(4).and_then(|n| pos.checked_add(n))
                .filter(|&end| end <= data.len())
                .ok_or_else(|| corrupt(format!(
                    "{what}: {count} entries at byte {pos} run past the end of the file ({} bytes)",
                    data.len())))
        };

        let version = r32(4, "version")?;
        if version != 0 {
            return Err(TwoBitError::BadVersion(version));
        }
        let seq_count = r32(8, "sequence count")? as usize;
        // reserved = r32(12), ignored

        // Each index entry is at least 5 bytes, so this also bounds the
        // allocations below by the file size.
        if seq_count > (data.len() - 16) / 5 {
            return Err(corrupt(format!(
                "header says {seq_count} sequences, more than a {}-byte file can index", data.len())));
        }

        let mut pos = 16usize;
        let mut sequences = Vec::with_capacity(seq_count);
        let mut name_to_index = HashMap::with_capacity(seq_count);

        for i in 0..seq_count {
            let entry = i + 1;
            let name_len = *data.get(pos).ok_or_else(|| corrupt(format!(
                "index entry {entry} of {seq_count} at byte {pos} is past the end of the file")))? as usize;
            pos += 1;
            let name_bytes = data.get(pos..pos + name_len).ok_or_else(|| corrupt(format!(
                "name of index entry {entry} at byte {pos} runs past the end of the file")))?;
            let name = std::str::from_utf8(name_bytes)
                .map_err(|_| corrupt(format!("name of index entry {entry} at byte {pos} is not UTF-8")))?
                .to_owned();
            pos += name_len;
            let offset = r32(pos, "record offset")? as u64;
            pos += 4;
            if let Some(&first) = name_to_index.get(&name) {
                return Err(TwoBitError::DuplicateName { name, first: first + 1, second: entry });
            }
            name_to_index.insert(name.clone(), i);
            sequences.push(SeqInfo { name, file_offset: offset, dna_size: 0 });
        }

        // Record: dnaSize nBlockCount nBlockStarts[] nBlockSizes[] maskBlockCount
        //         maskBlockStarts[] maskBlockSizes[] reserved packedDna[ceil(dnaSize/4)]
        let mut layouts = Vec::with_capacity(sequences.len());
        for seq in &mut sequences {
            let what = |part: &str| format!("{part} of '{}'", seq.name);
            let base = seq.file_offset as usize;
            seq.dna_size = r32(base, &what("dnaSize"))?;
            crate::seq::check_seq_len(&seq.name, seq.dna_size as usize).map_err(TwoBitError::TooLong)?;
            let n_block_count = r32(base + 4, &what("nBlockCount"))? as usize;
            let n_starts_off = base + 8;
            let n_sizes_off = skip_u32s(n_starts_off, n_block_count, &what("nBlockStarts"))?;
            let mask_count_off = skip_u32s(n_sizes_off, n_block_count, &what("nBlockSizes"))?;
            let mask_block_count = r32(mask_count_off, &what("maskBlockCount"))? as usize;
            let mask_sizes_off = skip_u32s(mask_count_off + 4, mask_block_count, &what("maskBlockStarts"))?;
            let reserved_off = skip_u32s(mask_sizes_off, mask_block_count, &what("maskBlockSizes"))?;
            r32(reserved_off, &what("reserved word"))?;
            let dna_off = reserved_off + 4;
            let packed_len = (seq.dna_size as usize).div_ceil(4);
            if dna_off + packed_len > data.len() {
                return Err(corrupt(format!(
                    "'{}' has {} bases ({packed_len} packed bytes at byte {dna_off}) but the file ends at byte {}",
                    seq.name, seq.dna_size, data.len())));
            }
            layouts.push(RecordLayout { n_block_count, n_starts_off, n_sizes_off, dna_off });
        }

        Ok(TwoBitFile { data, swap, sequences, layouts, name_to_index })
    }

    /// A u32 inside a record that `from_bytes` has already bounds-checked.
    fn r32(&self, pos: usize) -> u32 {
        let b = &self.data[pos..pos + 4];
        let v = u32::from_le_bytes(b.try_into().unwrap());
        if self.swap { v.swap_bytes() } else { v }
    }

    /// N block `bi` of a record as [start, end), clamped to the sequence
    /// (blocks at or past the end are empty, as before; the end is computed
    /// in u64 so a corrupt start + size cannot wrap).
    fn n_block(&self, lay: &RecordLayout, bi: usize, dna_size: u32) -> (u32, u32) {
        let n_start = self.r32(lay.n_starts_off + bi * 4);
        let n_size = self.r32(lay.n_sizes_off + bi * 4);
        let n_end = (n_start as u64 + n_size as u64).min(dna_size as u64) as u32;
        (n_start.min(n_end), n_end)
    }

    /// Decode the full sequence to BLASTNA (with sentinels).
    /// N-block positions use NCBI-compatible random bases (matching makeblastdb).
    pub fn get_full_sequence_blastna(&self, name: &str) -> Result<Vec<u8>, TwoBitError> {
        let idx = *self.name_to_index.get(name)
            .ok_or_else(|| TwoBitError::NotFound(name.to_owned()))?;
        self.decode_sequence(idx, 0, self.sequences[idx].dna_size)
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
        let (lay, dna_size) = (self.layouts[idx], self.sequences[idx].dna_size);
        if lay.n_block_count == 0 {
            return Ok(Vec::new());
        }
        let mut mask = vec![0u8; dna_size as usize];
        for bi in 0..lay.n_block_count {
            let (n_start, n_end) = self.n_block(&lay, bi, dna_size);
            mask[n_start as usize..n_end as usize].fill(14);
        }
        Ok(mask)
    }

    fn decode_sequence(&self, idx: usize, start: u32, len: u32) -> Result<Vec<u8>, TwoBitError> {
        let (lay, dna_size) = (self.layouts[idx], self.sequences[idx].dna_size);

        let actual_end = (start as u64 + len as u64).min(dna_size as u64) as u32;
        let actual_len = actual_end.saturating_sub(start);

        let mut out = vec![NUCL_SENTINEL; actual_len as usize + 2];

        // Decode packed DNA into BLASTNA (from_bytes checked the packed run fits)
        for i in 0..actual_len {
            let pos = start + i;
            let byte_idx = (pos / 4) as usize;
            let bit_shift = 6 - ((pos % 4) * 2);
            let two_bits = ((self.data[lay.dna_off + byte_idx] >> bit_shift) & 0x3) as usize;
            out[(i + 1) as usize] = UCSC2BIT_TO_BLASTNA[two_bits];
        }

        // Replace N-block positions with NCBI-compatible pseudo-random bases.
        // .2bit files only carry N-level ambiguity (faToTwoBit converts all
        // non-ACGT characters to N), so we randomise uniformly over all four
        // bases, matching makeblastdb's GetRand()&3 for N positions.
        if lay.n_block_count > 0 {
            let mut rng = NcbiRandom::new(dna_size);
            for bi in 0..lay.n_block_count {
                let (n_start, n_end) = self.n_block(&lay, bi, dna_size);
                for pos in n_start..n_end {
                    let rand_base = rng.next_base();
                    if pos >= start && pos < actual_end {
                        out[(pos - start + 1) as usize] = rand_base;
                    }
                }
            }
        }

        Ok(out)
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

    /// A valid 2bit byte stream holding `records` (name, sequence of ACGT and
    /// N).  Runs of N become N blocks (packed as T, like faToTwoBit); there
    /// are no mask blocks.
    pub(crate) fn twobit_bytes(records: &[(&str, &str)]) -> Vec<u8> {
        let code = |b: u8| match b { b'T' | b'N' => 0u8, b'C' => 1, b'A' => 2, b'G' => 3, _ => panic!("ACGTN only") };
        let n_runs = |seq: &str| -> Vec<(u32, u32)> {
            let mut runs: Vec<(u32, u32)> = Vec::new();
            for (i, b) in seq.bytes().enumerate() {
                if b == b'N' {
                    match runs.last_mut() {
                        Some((st, sz)) if *st + *sz == i as u32 => *sz += 1,
                        _ => runs.push((i as u32, 1)),
                    }
                }
            }
            runs
        };
        let mut v: Vec<u8> = Vec::new();
        v.extend_from_slice(&TWOBIT_MAGIC.to_le_bytes());
        v.extend_from_slice(&0u32.to_le_bytes());
        v.extend_from_slice(&(records.len() as u32).to_le_bytes());
        v.extend_from_slice(&0u32.to_le_bytes());
        let index_len: usize = records.iter().map(|(n, _)| 1 + n.len() + 4).sum();
        let mut offset = 16 + index_len;
        let mut bodies = Vec::new();
        for (name, seq) in records {
            v.push(name.len() as u8);
            v.extend_from_slice(name.as_bytes());
            v.extend_from_slice(&(offset as u32).to_le_bytes());
            let mut body = Vec::new();
            body.extend_from_slice(&(seq.len() as u32).to_le_bytes());
            let runs = n_runs(seq);
            body.extend_from_slice(&(runs.len() as u32).to_le_bytes()); // nBlockCount
            for (st, _) in &runs { body.extend_from_slice(&st.to_le_bytes()); }
            for (_, sz) in &runs { body.extend_from_slice(&sz.to_le_bytes()); }
            body.extend_from_slice(&0u32.to_le_bytes()); // maskBlockCount
            body.extend_from_slice(&0u32.to_le_bytes()); // reserved
            for chunk in seq.as_bytes().chunks(4) {
                let mut byte = 0u8;
                for (k, &b) in chunk.iter().enumerate() {
                    byte |= code(b) << (6 - 2 * k);
                }
                body.push(byte);
            }
            offset += body.len();
            bodies.push(body);
        }
        for b in bodies {
            v.extend_from_slice(&b);
        }
        v
    }

    #[test]
    fn twobit_bytes_round_trips() {
        let tb = TwoBitFile::from_bytes(twobit_bytes(&[("a", "ACGTA"), ("b", "TTGCA")])).unwrap();
        assert_eq!(tb.get_full_sequence_blastna("a").unwrap(), vec![15, 0, 1, 2, 3, 0, 15]);
        assert_eq!(tb.get_full_sequence_blastna("b").unwrap(), vec![15, 3, 3, 2, 1, 0, 15]);
    }

    /// Subjects are looked up by name, so a duplicate would leave one record
    /// unsearched and the other searched twice: reject it at load.
    #[test]
    fn duplicate_names_are_rejected() {
        let e = TwoBitFile::from_bytes(twobit_bytes(&[("x", "ACGT"), ("y", "ACGT"), ("x", "TTTT")]))
            .err().expect("duplicate accepted");
        assert!(matches!(e, TwoBitError::DuplicateName { ref name, first: 1, second: 3 } if name == "x"), "{e}");
    }

    fn two_records() -> Vec<u8> {
        twobit_bytes(&[("a", "ACGTNNNACGTAC"), ("b", "TTGCANNNNNNNNNNNGGA")])
    }

    #[test]
    fn n_blocks_round_trip() {
        let tb = TwoBitFile::from_bytes(two_records()).unwrap();
        let mask = tb.get_n_mask("a").unwrap();
        assert_eq!(mask, [0, 0, 0, 0, 14, 14, 14, 0, 0, 0, 0, 0, 0]);
        let seq = tb.get_full_sequence_blastna("b").unwrap();
        assert_eq!(seq.len(), 19 + 2);
        assert_eq!(&seq[..6], &[15, 3, 3, 2, 1, 0]);
    }

    /// Every truncation of a valid file is an error, never a panic, and never
    /// a database with fewer sequences (3.0.12 panicked at several offsets and
    /// returned zero sequences for a header-only file).
    #[test]
    fn every_truncation_is_an_error() {
        let full = two_records();
        assert!(TwoBitFile::from_bytes(full.clone()).is_ok());
        for len in 0..full.len() {
            let r = std::panic::catch_unwind(|| TwoBitFile::from_bytes(full[..len].to_vec()));
            match r {
                Ok(Err(_)) => {}
                Ok(Ok(_)) => panic!("{len}-byte prefix of a {}-byte file loaded", full.len()),
                Err(_) => panic!("{len}-byte prefix panicked"),
            }
        }
    }

    fn put_u32(v: &mut [u8], pos: usize, x: u32) {
        v[pos..pos + 4].copy_from_slice(&x.to_le_bytes());
    }

    /// Corrupt counts and offsets are errors that say what is wrong, and none
    /// allocates by the corrupt value (dnaSize 0xFFFFFFF0 used to allocate 4 GB).
    #[test]
    fn corrupt_fields_are_errors() {
        let full = two_records();
        let rec_a = u32::from_le_bytes(full[16 + 2..16 + 6].try_into().unwrap()) as usize; // after "a"'s 1-byte len + name
        let cases: Vec<(&str, Vec<u8>, &str)> = vec![
            ("huge seqCount", { let mut v = full.clone(); put_u32(&mut v, 8, u32::MAX); v }, "sequences"),
            // the bogus 3rd entry parses out of record data; a later check catches it
            ("seqCount past index", { let mut v = full.clone(); put_u32(&mut v, 8, 3); v }, ""),
            ("offset past end", { let mut v = full.clone(); put_u32(&mut v, 18, 1 << 30); v }, "dnaSize of 'a'"),
            ("huge dnaSize", { let mut v = full.clone(); put_u32(&mut v, rec_a, 0x7FFF_FFF0); v }, "bases"),
            ("huge nBlockCount", { let mut v = full.clone(); put_u32(&mut v, rec_a + 4, 0x4000_0000); v }, "nBlockStarts of 'a'"),
            ("header only", { let mut v = full[..16].to_vec(); put_u32(&mut v, 8, 1); v }, "sequences"),
        ];
        for (label, bytes, needle) in cases {
            match TwoBitFile::from_bytes(bytes) {
                Err(e @ TwoBitError::Corrupt(_)) => {
                    assert!(e.to_string().contains(needle), "{label}: {e}");
                }
                Err(e) => panic!("{label}: wrong error {e}"),
                Ok(_) => panic!("{label}: loaded"),
            }
        }
    }

    /// A 2bit can hold up to 2^32 - 1 bases per sequence; the engine's
    /// coordinates stop at 2^31 - 1 (3.0.12 wrapped and found nothing).
    #[test]
    fn sequences_over_2_31_are_refused() {
        let mut v = twobit_bytes(&[("big", "ACGT")]);
        let rec = u32::from_le_bytes(v[16 + 4..16 + 8].try_into().unwrap()) as usize;
        put_u32(&mut v, rec, 1 << 31);
        match TwoBitFile::from_bytes(v) {
            Err(e @ TwoBitError::TooLong(_)) => assert!(e.to_string().contains("'big' has 2147483648 bases"), "{e}"),
            Err(e) => panic!("wrong error: {e}"),
            Ok(_) => panic!("loaded"),
        }
        assert!(crate::seq::check_seq_len("q", i32::MAX as usize).is_ok());
        assert!(crate::seq::check_seq_len("q", i32::MAX as usize + 1).is_err());
    }

    /// N blocks at or past the sequence end are clamped, as before, and a
    /// start + size that overflows u32 cannot wrap into the sequence.
    #[test]
    fn n_blocks_past_the_end_are_clamped() {
        let mut v = twobit_bytes(&[("a", "ACGTNNAC")]);
        let rec = u32::from_le_bytes(v[16 + 2..16 + 6].try_into().unwrap()) as usize;
        put_u32(&mut v, rec + 8, 6);            // block start 6
        put_u32(&mut v, rec + 12, u32::MAX);    // block size: start + size wraps in u32
        let tb = TwoBitFile::from_bytes(v).unwrap();
        assert_eq!(tb.get_n_mask("a").unwrap(), [0, 0, 0, 0, 0, 0, 14, 14]);
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
