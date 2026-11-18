// blast_types.rs — ports of the NCBI BLAST core data structures still in use.
//
// Primary source (NCBI BLAST 2.17.0+):
//   blast_def.h        → BlastOffsetPair, BLAST_SequenceBlk (SeqBlk here)

// ── blast_def.h ──────────────────────────────────────────────────────────────

/// COMPRESSION_RATIO: 4 nucleotide bases packed into 1 byte.
pub const COMPRESSION_RATIO: i32 = 4;

/// Analogous to BLAST_SequenceBlk.
///
/// NCBI maintains two forms of a sequence in the same struct:
///   sequence       – primary form (BLASTNA 1-byte/base for query; packed
///                    NCBI2NA 4-bases/byte for database subjects in prelim phase)
///   sequence_start – one position before `sequence`; holds the sentinel byte
///
/// For the query, a third form is also maintained:
///   compressed_nuc_seq       – 4-to-1 packed NCBI2NA built from the BLASTNA
///                              query; used by the ungapped extension scanner
///   compressed_nuc_seq_start – allocation root (3 bytes before compressed_nuc_seq)
///
/// In this Rust port we store the BLASTNA form and the packed form as separate
/// Vecs so that lifetime / ownership is explicit.  The `sequence()` helper
/// returns a slice starting at base-index 0 (i.e., past the sentinel byte),
/// matching C's `seq_blk->sequence`.
pub struct SeqBlk {
    /// BLASTNA (1-byte-per-base) storage.
    /// Layout: data[0] = sentinel (0), data[1..=length] = bases.
    /// `sequence` in NCBI C code points to data[1].
    pub blastna: Vec<u8>,

    /// Packed NCBI2NA storage for this sequence (built by
    /// `blast_compress_blastna_sequence`).
    ///
    /// Layout mirrors NCBI's compressed_nuc_seq:
    ///   packed[0..3] = zero padding before the sequence
    ///   packed[3..]  = packed bases; packed[3+j] covers bases 4j..4j+3 with
    ///                  base 4j in bits 7:6 (MSB) and base 4j+3 in bits 1:0;
    ///                  a partial last byte is left-justified
    ///   then PACKED_TAIL_PAD zero bytes (blast::util).
    ///
    /// `compressed_nuc_seq` in C points to packed[3].
    pub packed: Vec<u8>,

    /// Sequence length in bases (not including sentinel).
    pub length: i32,
}

impl SeqBlk {
    /// Slice of BLASTNA bases, analogous to `seq_blk->sequence` in NCBI.
    /// Index 0 corresponds to the first base.
    #[inline]
    pub fn sequence(&self) -> &[u8] {
        &self.blastna[1..] // skip sentinel
    }

    /// Mutable BLASTNA bases.
    #[inline]
    pub fn sequence_mut(&mut self) -> &mut [u8] {
        &mut self.blastna[1..]
    }

    /// `compressed_nuc_seq` pointer equivalent: slice starting at the
    /// packed base-0 byte (index 3 in `self.packed`).
    #[inline]
    pub fn compressed_nuc_seq(&self) -> &[u8] {
        &self.packed[3..]
    }
}

/// BlastOffsetPair (blast_def.h).
/// Holds the query/subject offsets of an initial word match.
#[derive(Clone, Copy, Debug, Default)]
pub struct BlastOffsetPair {
    pub q_off: u32,
    pub s_off: u32,
}
