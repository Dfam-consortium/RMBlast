//! Nucleotide word lookup table for seeding alignments.
//!
//! Builds a hash-based lookup from `word_size`-mer (ncbi2na, 2-bit packed)
//! to list of query positions.  The scan function then slides over the subject
//! sequence to find seed matches.
//!
//! TODO: implement the full lookup table + scan logic (port of
//! blast_nalookup.c / blast_nascan.c).

pub mod na_lookup;
pub use na_lookup::{NaLookupTable, WordHit};
