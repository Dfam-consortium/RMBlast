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

// The RNG is CRandom (seq::twobit::NcbiRandom), seeded on first use.
use crate::seq::twobit::LazyRandom;

// ── BLASTNA → NCBI4NA table ───────────────────────────────────────────────────
//
// BLASTNA: A=0 C=1 G=2 T=3 R=4 Y=5 M=6 K=7 W=8 S=9 B=10 D=11 H=12 V=13 N=14 gap=15
//          (encoding::BLASTNA_TO_IUPAC order)
// NCBI4NA: bit0=A bit1=C bit2=G bit3=T; N=15 gap=0
// Copied from blast_encoding.c.  3.0.12's table assumed another BLASTNA order
// (W/S swapped, B/D/H/V permuted).  The CLI never reached it, since both
// database decoders resolve ambiguity codes before packing, but a library
// caller packing raw IUPAC codes drew from the wrong base sets.
const BLASTNA_TO_NCBI4NA: [u8; 16] = [1, 2, 4, 8, 5, 10, 3, 12, 9, 6, 14, 13, 11, 7, 15, 0];

/// Convert a single BLASTNA base to NCBI2NA, advancing the RNG for ambiguous bases.
///
/// Mirrors s_Ncbi4naToNcbi2na (blast_objmgr_tools.cpp):
///   - Unambiguous (A/C/G/T): returns the NCBI2NA code directly (no RNG call).
///   - N (NCBI4NA=15) or gap (NCBI4NA=0): returns `rng.GetRand() & 0x3`.
///   - Other ambiguous (2–3 bits set): returns `rng.GetRand() % bitcount`-th set-bit index.
#[inline]
fn blastna_to_ncbi2na(b: u8, rng: &mut LazyRandom) -> u8 {
    if b < 4 { return b; }  // A=0, C=1, G=2, T=3: BLASTNA == NCBI2NA
    let rng = rng.get();
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
/// Zero bytes after the last packed base.  The subject scanners (nascan.rs,
/// mblookup.rs) read a word as up to four bytes from the byte holding the
/// scan offset (`s_byte + 3`).  Since the last offset scanned leaves a whole
/// word before the end, 1 byte is enough today; 3 covers a read of
/// `s_byte + 3` from the last data byte, whatever the word length.  Tested by
/// match_at_subject_end_found_for_every_length (fails at 0).
pub const PACKED_TAIL_PAD: usize = 3;

pub fn blast_compress_blastna_sequence(seq_blk: &mut SeqBlk) {
    let len = seq_blk.length as usize;

    // 3 pre-sequence padding bytes + ceil(len/4) data bytes + PACKED_TAIL_PAD.
    // 3.0.12 allocated len + 3 bytes, four times what the bases need.
    seq_blk.packed = vec![0u8; 3 + len.div_ceil(4) + PACKED_TAIL_PAD];

    let old_seq = &seq_blk.blastna[1..1 + len];

    // Seed RNG with sequence length to match s_Ncbi4naToNcbi2na(ncbi4na, base_length, ...).
    let mut rng = LazyRandom::new(len as u32);

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encoding::BLASTNA_TO_IUPAC;

    /// Each BLASTNA code maps to the NCBI4NA bit set of its IUPAC meaning.
    #[test]
    fn blastna_to_ncbi4na_matches_iupac() {
        for (code, &letter) in BLASTNA_TO_IUPAC.iter().enumerate() {
            let bases: &[u8] = match letter {
                b'A' => b"A", b'C' => b"C", b'G' => b"G", b'T' => b"T",
                b'R' => b"AG", b'Y' => b"CT", b'M' => b"AC", b'K' => b"GT",
                b'W' => b"AT", b'S' => b"CG", b'B' => b"CGT", b'D' => b"AGT",
                b'H' => b"ACT", b'V' => b"ACG", b'N' => b"ACGT", b'-' => b"",
                _ => unreachable!(),
            };
            let bits = bases.iter().fold(0u8, |m, b| m | 1 << b"ACGT".iter().position(|x| x == b).unwrap());
            assert_eq!(BLASTNA_TO_NCBI4NA[code], bits, "BLASTNA {code} ({})", letter as char);
        }
    }
}
