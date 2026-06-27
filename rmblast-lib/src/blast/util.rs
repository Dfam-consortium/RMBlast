// blast_util.rs — sequence utility functions.
//
// blast_compress_blastna_sequence replicates makeblastdb's s_Ncbi4naToNcbi2na
// (blast_objmgr_tools.cpp): ambiguous bases get pseudo-random 2-bit values from
// an LFG RNG seeded by the sequence length, exactly as stored in NCBI BLAST databases.

use super::types::SeqBlk;

// ── NCBI2NA encoding / unpacking ─────────────────────────────────────────────

/// NCBI2NA_MASK (blast_util.h line 52).
pub const NCBI2NA_MASK: u8 = 0x03;

/// NCBI2NA_UNPACK_BASE(x, N) — extract 2-bit base N from packed byte x.
///
/// N=3 → bits 7:6 (MSB, first base in the byte)
/// N=0 → bits 1:0 (LSB, last base in the byte)
///
/// Mirrors: `#define NCBI2NA_UNPACK_BASE(x, N) (((x)>>(2*(N))) & NCBI2NA_MASK)`
#[inline(always)]
pub fn ncbi2na_unpack_base(x: u8, n: u32) -> u8 {
    (x >> (2 * n)) & NCBI2NA_MASK
}

// ── NCBI LFG random number generator ─────────────────────────────────────────
//
// Mirrors CRandom (util/random_gen.hpp + util/random_gen.cpp) from NCBI C++ Toolkit.
// Lagged Fibonacci generator with lags 33 and 13 (STATE_SIZE - STATE_OFFSET - 1 = 20
// from the initial rk/rj gap of 32-12=20).

struct NcbiRandom {
    state: [u32; 33],
    rk:    i32,
    rj:    i32,
}

impl NcbiRandom {
    const STATE_SIZE:   usize = 33;
    const STATE_OFFSET: i32   = 12;

    fn new(seed: u32) -> Self {
        let mut rng = NcbiRandom { state: [0u32; 33], rk: 32, rj: Self::STATE_OFFSET };
        rng.state[0] = seed;
        for i in 1..Self::STATE_SIZE {
            rng.state[i] = 1103515245u32.wrapping_mul(rng.state[i - 1]).wrapping_add(12345);
        }
        rng.rk = (Self::STATE_SIZE - 1) as i32;
        rng.rj = Self::STATE_OFFSET;
        for _ in 0..10 * Self::STATE_SIZE {
            rng.x_get_rand32bits();
        }
        rng
    }

    #[inline]
    fn x_get_rand32bits(&mut self) -> u32 {
        let r = self.state[self.rk as usize].wrapping_add(self.state[self.rj as usize]);
        self.state[self.rk as usize] = r;
        self.rk -= 1;
        self.rj -= 1;
        if self.rk < 0 {
            self.rk = (Self::STATE_SIZE - 1) as i32;
        } else if self.rj < 0 {
            self.rj = (Self::STATE_SIZE - 1) as i32;
        }
        r
    }

    #[inline]
    fn get_rand(&mut self) -> u32 {
        self.x_get_rand32bits() >> 1
    }
}

// ── BLASTNA → NCBI4NA table ───────────────────────────────────────────────────
//
// BLASTNA: A=0 C=1 G=2 T=3 R=4 Y=5 M=6 K=7 S=8 W=9 H=10 B=11 V=12 D=13 N=14 gap=15
// NCBI4NA: bit0=A bit1=C bit2=G bit3=T; N=15 gap=0
const BLASTNA_TO_NCBI4NA: [u8; 16] = [1, 2, 4, 8, 5, 10, 3, 12, 6, 9, 11, 14, 7, 13, 15, 0];

/// Convert a single BLASTNA base to NCBI2NA, advancing the RNG for ambiguous bases.
///
/// Mirrors s_Ncbi4naToNcbi2na (blast_objmgr_tools.cpp):
///   - Unambiguous (A/C/G/T): returns the NCBI2NA code directly (no RNG call).
///   - N (NCBI4NA=15) or gap (NCBI4NA=0): returns `rng.GetRand() & 0x3`.
///   - Other ambiguous (2–3 bits set): returns `rng.GetRand() % bitcount`-th set-bit index.
#[inline]
fn blastna_to_ncbi2na(b: u8, rng: &mut NcbiRandom) -> u8 {
    if b < 4 { return b; }  // A=0, C=1, G=2, T=3: BLASTNA == NCBI2NA
    let ncbi4na = BLASTNA_TO_NCBI4NA[b as usize];
    if ncbi4na == 0 || ncbi4na == 15 {
        return (rng.get_rand() & 0x3) as u8;
    }
    // Multi-bit ambiguous code: pick one of the set bits uniformly.
    let bitcount = ncbi4na.count_ones();
    let mut pick  = rng.get_rand() % bitcount;
    for j in 0u8..4 {
        if ncbi4na & (1 << j) != 0 {
            if pick == 0 { return j; }
            pick -= 1;
        }
    }
    0
}

// ── blast_compress_blastna_sequence ───────────────────────────────────────────

/// Compress a BLASTNA sequence (1 byte/base) into 4-to-1 packed NCBI2NA.
///
/// Replicates makeblastdb's ambiguous-base encoding: the NCBI LFG RNG is seeded
/// with the sequence length, and each ambiguous base (BLASTNA >= 4) advances the
/// RNG to produce a deterministic pseudo-random 2-bit value — identical to what
/// makeblastdb stores in the .nsq packed sequence.
///
/// Format: packed[3+j] covers bases 4j..4j+3, base 4j in bits 7:6 (MSB),
/// base 4j+3 in bits 1:0 (LSB).  A partial final chunk is left-justified.
/// packed[0..2] are zeroed (pre-sequence padding).
pub fn blast_compress_blastna_sequence(seq_blk: &mut SeqBlk) {
    let len = seq_blk.length as usize;

    // 3 pre-sequence padding bytes + ceil(len/4) data bytes.
    seq_blk.packed = vec![0u8; len + 3];

    let old_seq = &seq_blk.blastna[1..1 + len];

    // Seed RNG with sequence length to match s_Ncbi4naToNcbi2na(ncbi4na, base_length, ...).
    let mut rng = NcbiRandom::new(len as u32);

    // Full 4-base chunks.
    for j in 0..(len / 4) {
        let b0 = blastna_to_ncbi2na(old_seq[4*j],   &mut rng);
        let b1 = blastna_to_ncbi2na(old_seq[4*j+1], &mut rng);
        let b2 = blastna_to_ncbi2na(old_seq[4*j+2], &mut rng);
        let b3 = blastna_to_ncbi2na(old_seq[4*j+3], &mut rng);
        seq_blk.packed[3 + j] = (b0 << 6) | (b1 << 4) | (b2 << 2) | b3;
    }
    // Partial last chunk: left-justify, right-pad with 0.
    let rem = len % 4;
    if rem > 0 {
        let base = len / 4;
        let mut byte = 0u8;
        for k in 0..rem {
            byte = (byte << 2) | blastna_to_ncbi2na(old_seq[4 * base + k], &mut rng);
        }
        seq_blk.packed[3 + base] = byte << (2 * (4 - rem) as u32);
    }
}

// ── SeqBlk construction ──────────────────────────────────────────────────────

/// Build a SeqBlk from a BLASTNA slice (1 byte/base, no sentinels).
///
/// Prepends a zero sentinel byte (matches `sequence_start[0] = 0` in NCBI),
/// so `blk.sequence()[i]` = base at position i.
pub fn seqblk_from_blastna(bases: &[u8]) -> SeqBlk {
    let length = bases.len() as i32;
    let mut blastna = Vec::with_capacity(bases.len() + 1);
    blastna.push(0u8); // sentinel
    blastna.extend_from_slice(bases);
    SeqBlk { blastna, packed: Vec::new(), length }
}

/// Build a SeqBlk, also compressing to NCBI2NA packed form.
pub fn seqblk_from_blastna_with_compress(bases: &[u8]) -> SeqBlk {
    let mut blk = seqblk_from_blastna(bases);
    blast_compress_blastna_sequence(&mut blk);
    blk
}
