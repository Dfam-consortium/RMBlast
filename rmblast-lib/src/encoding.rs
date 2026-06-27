//! BLASTNA sequence encoding used throughout the alignment engine.
//!
//! BLASTNA (16-symbol alphabet):
//!   A=0 C=1 G=2 T=3 R=4 Y=5 M=6 K=7 W=8 S=9 B=10 D=11 H=12 V=13 N=14 GAP=15
//!
//! The first four symbols (0-3) match ncbi2na, so the 2-bit packed lookup
//! table and the 16×16 scoring matrix share the same encoding for ACGT.
//!
//! UCSC 2bit packing: T=0b00 C=0b01 A=0b10 G=0b11 (high bits = first base).
//! Mapping to BLASTNA: T→3, C→1, A→0, G→2.

/// BLASTNA alphabet size (including gap/N/ambiguous codes).
pub const BLASTNA_SIZE: usize = 16;
/// Unambiguous base count (A, C, G, T).
pub const BLAST2NA_SIZE: usize = 4;
/// Sentinel value for sequence ends (same as GAP in BLASTNA).
pub const NUCL_SENTINEL: u8 = 15;

/// BLASTNA indices for unambiguous bases.
pub const BLASTNA_A: u8 = 0;
pub const BLASTNA_C: u8 = 1;
pub const BLASTNA_G: u8 = 2;
pub const BLASTNA_T: u8 = 3;

/// Convert IUPAC ASCII character to BLASTNA index.
/// Unknown characters map to N (14).
pub const IUPAC_TO_BLASTNA: [u8; 128] = {
    let mut t = [14u8; 128]; // default = N
    // Gap
    t[b'-' as usize] = 15;
    // Unambiguous
    t[b'A' as usize] = 0;  t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;  t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;  t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;  t[b't' as usize] = 3;
    t[b'U' as usize] = 3;  t[b'u' as usize] = 3; // RNA U → T
    // Ambiguous (IUPAC)
    t[b'R' as usize] = 4;  t[b'r' as usize] = 4;
    t[b'Y' as usize] = 5;  t[b'y' as usize] = 5;
    t[b'M' as usize] = 6;  t[b'm' as usize] = 6;
    t[b'K' as usize] = 7;  t[b'k' as usize] = 7;
    t[b'W' as usize] = 8;  t[b'w' as usize] = 8;
    t[b'S' as usize] = 9;  t[b's' as usize] = 9;
    t[b'B' as usize] = 10; t[b'b' as usize] = 10;
    t[b'D' as usize] = 11; t[b'd' as usize] = 11;
    t[b'H' as usize] = 12; t[b'h' as usize] = 12;
    t[b'V' as usize] = 13; t[b'v' as usize] = 13;
    t[b'N' as usize] = 14; t[b'n' as usize] = 14;
    // In a *sequence*, 'X' is treated as N: NCBI's FASTA reader encodes a query/
    // subject 'X' to BLASTNA N (it is rendered 'N' and aligned as an ordinary
    // ambiguous base).  NOTE: this is DIFFERENT from how 'X' is mapped when it
    // appears as a *matrix* column/row header -- there NCBI's IUPACNA_TO_BLASTNA
    // sends 'X'->15 (the gap slot) so it cannot clobber the N matrix row/column.
    // The matrix parser (matrix.rs) special-cases 'X'->15 itself; this table is
    // only for sequence encoding, so here 'X' must map to N (14).  (See bug #33
    // for the matrix side and the regression it caused when these two were
    // conflated into one X->15 mapping.)
    t[b'X' as usize] = 14; t[b'x' as usize] = 14;
    t
};

/// Convert BLASTNA index to IUPAC ASCII character.
pub const BLASTNA_TO_IUPAC: [u8; 16] = *b"ACGTRYMKWSBDHVN-";

/// Convert UCSC 2bit 2-bit base value (T=0,C=1,A=2,G=3) to BLASTNA.
pub const UCSC2BIT_TO_BLASTNA: [u8; 4] = [
    3, // 0b00 T → BLASTNA T=3
    1, // 0b01 C → BLASTNA C=1
    0, // 0b10 A → BLASTNA A=0
    2, // 0b11 G → BLASTNA G=2
];

/// Complement a BLASTNA base.
pub const BLASTNA_COMPLEMENT: [u8; 16] = [
    3,  // A → T
    2,  // C → G
    1,  // G → C
    0,  // T → A
    5,  // R(AG) → Y(CT)
    4,  // Y(CT) → R(AG)
    7,  // M(AC) → K(GT)
    6,  // K(GT) → M(AC)
    8,  // W(AT) → W(AT)
    9,  // S(CG) → S(CG)
    13, // B(CGT) → V(ACG): complement of {C,G,T} = {G,C,A} = V
    12, // D(AGT) → H(ACT): complement of {A,G,T} = {T,C,A} = H
    11, // H(ACT) → D(AGT): complement of {A,C,T} = {T,G,A} = D
    10, // V(ACG) → B(CGT): complement of {A,C,G} = {T,G,C} = B
    14, // N → N
    15, // GAP → GAP
];


/// Allowed bases for each BLASTNA ambiguity code, in NCBI4NA bit order (A,C,G,T),
/// for makeblastdb-compatible constrained-random base assignment.
///
/// Mirrors `CWriteDB_Impl::CMaskedRanges`-side `x_Random` in NCBI
/// `objtools/blast/seqdb_writer/writedb_convert.cpp`:
///   N (0xF):    `GetRand() & 0x3`            → any of {A,C,G,T}
///   otherwise:  `pick = GetRand() % bitcount`, then return the `pick`-th set bit,
///               scanning bits 0..3 (= A,C,G,T).
///
/// Layout per code: `[count, base0, base1, base2]`.  The caller uses
/// `bases[1 + (GetRand() % count)]` (and handles N=14 specially as `GetRand() & 3`).
///
/// NOTE: a previous table indexed `[GetRand() & 3]` with wraparound, which is WRONG
/// for the 3-bit codes B/D/H/V because `(x & 3) % 3 != x % 3`.  That mis-assigned the
/// random base at B/D/H/V positions (e.g. `D`), desyncing seed words from NCBI's
/// blastdb and dropping seeds in ambiguity-dense subjects.  `% count` matches NCBI.
pub const NCBI_AMBIG_ALLOWED: [[u8; 4]; 16] = [
    [0, 0, 0, 0], //  0: A  (unambiguous, unused)
    [0, 0, 0, 0], //  1: C  (unambiguous, unused)
    [0, 0, 0, 0], //  2: G  (unambiguous, unused)
    [0, 0, 0, 0], //  3: T  (unambiguous, unused)
    [2, 0, 2, 0], //  4: R (A,G)
    [2, 1, 3, 0], //  5: Y (C,T)
    [2, 0, 1, 0], //  6: M (A,C)
    [2, 2, 3, 0], //  7: K (G,T)
    [2, 0, 3, 0], //  8: W (A,T)
    [2, 1, 2, 0], //  9: S (C,G)
    [3, 1, 2, 3], // 10: B (C,G,T)
    [3, 0, 2, 3], // 11: D (A,G,T)
    [3, 0, 1, 3], // 12: H (A,C,T)
    [3, 0, 1, 2], // 13: V (A,C,G)
    [4, 0, 0, 0], // 14: N (special-cased: GetRand() & 3)
    [0, 0, 0, 0], // 15: GAP/sentinel (unused)
];

/// Encode a byte slice of IUPAC ASCII into a BLASTNA Vec.
/// Prepends and appends a sentinel (15) for the BLAST engine.
pub fn encode_iupac(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len() + 2);
    out.push(NUCL_SENTINEL);
    for &b in seq {
        out.push(if (b as usize) < 128 { IUPAC_TO_BLASTNA[b as usize] } else { 14 });
    }
    out.push(NUCL_SENTINEL);
    out
}

/// Decode a BLASTNA slice back to IUPAC ASCII (no sentinels).
pub fn decode_blastna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .filter(|&&b| b != NUCL_SENTINEL)
        .map(|&b| BLASTNA_TO_IUPAC[(b & 15) as usize])
        .collect()
}

/// Reverse-complement a BLASTNA encoded sequence (excluding sentinels).
pub fn revcomp_blastna(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| BLASTNA_COMPLEMENT[(b & 15) as usize])
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_roundtrip() {
        let input = b"ACGTRYN";
        let encoded = encode_iupac(input);
        // encoded[0] and encoded[last] are sentinels
        assert_eq!(encoded[0], NUCL_SENTINEL);
        assert_eq!(*encoded.last().unwrap(), NUCL_SENTINEL);
        let inner: Vec<u8> = encoded[1..encoded.len() - 1].to_vec();
        assert_eq!(inner, vec![0, 1, 2, 3, 4, 5, 14]);
    }

    #[test]
    fn test_revcomp() {
        // ACGT → revcomp → ACGT
        let seq = vec![0u8, 1, 2, 3]; // A C G T
        let rc = revcomp_blastna(&seq);
        assert_eq!(rc, vec![0, 1, 2, 3]); // ACGT revcomp = ACGT
    }
}
