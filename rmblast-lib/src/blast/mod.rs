// blast/mod.rs — faithful port of NCBI BLAST core algorithms.
//
// Module layout mirrors NCBI's C source file structure:
//   types.rs     → blast_def.h (BlastOffsetPair, BLAST_SequenceBlk)
//   util.rs      → blast_util.c  (BlastCompressBlastnaSequence, NCBI2NA_UNPACK_BASE)
//   nalookup.rs  → blast_lookup.c  (BlastLookupAddWordHit,
//                                   BlastLookupIndexQueryExactMatches)
//                  blast_nalookup.c (s_BlastNaLookupFinalize,
//                                    BlastNaLookupTableNew)
//                  blast_lookup.h / blast_nalookup.h (structs, PV macros)
//   nascan.rs    → blast_nascan.c  (s_BlastLookupGetNumHits,
//                                   s_BlastLookupRetrieve,
//                                   s_BlastNaScanSubject_8_4,
//                                   s_BlastNaScanSubject_Any)

pub mod types;
pub mod util;
pub mod nalookup;
pub mod nascan;
pub mod mblookup;
