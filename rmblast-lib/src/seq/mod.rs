pub mod fasta;
pub mod fasta_db;
pub mod twobit;

pub use fasta::{FastaRecord, FastaReader};
pub use fasta_db::FastaDb;
pub use twobit::{SeqInfo, TwoBitFile, TwoBitSeqIter};

use std::io::Read;
use std::collections::HashMap;
use std::sync::Arc;
use anyhow::Result;

/// Unified subject-database handle that accepts either a FASTA or a 2bit file.
///
/// Format detection order:
///   1. Extension (case-insensitive): `.fa`, `.fasta`, `.fas`, `.fna`, `.fn`
///      → FASTA; `.2bit` → 2bit.
///   2. For files with no recognised extension, peek at the first byte:
///      `>` → FASTA; otherwise 2bit.
///
/// When a .2bit database is opened a note is printed to stderr explaining
/// that non-N IUB codes are not preserved and FASTA is preferred for full
/// IUB fidelity.
enum SubjectBackend {
    TwoBit(TwoBitFile),
    Fasta(FastaDb),
}

impl SubjectBackend {
    fn sequences(&self) -> &[SeqInfo] {
        match self {
            SubjectBackend::TwoBit(db) => &db.sequences,
            SubjectBackend::Fasta(db) => &db.sequences,
        }
    }
    fn decode(&self, name: &str) -> Result<Vec<u8>> {
        match self {
            SubjectBackend::TwoBit(db) => db
                .get_full_sequence_blastna(name)
                .map_err(|e| anyhow::anyhow!("{}", e)),
            SubjectBackend::Fasta(db) => db
                .get_full_sequence_blastna(name)
                .map_err(|e| anyhow::anyhow!("{}", e)),
        }
    }
    fn nmask(&self, name: &str) -> Vec<u8> {
        match self {
            SubjectBackend::TwoBit(db) => db.get_n_mask(name).unwrap_or_default(),
            SubjectBackend::Fasta(db) => db.get_n_mask(name),
        }
    }
}

/// Unified subject-database handle (FASTA or 2bit) with a **decoded-once cache**.
///
/// The backend stores the raw/packed sequences; at open time every subject is
/// decoded to BLASTNA (ambiguities resolved deterministically) and its n_mask is
/// computed once, then stored as `Arc<[u8]>`.  `get_full_sequence_blastna` /
/// `get_n_mask` hand out cheap `Arc` clones of the cached buffers instead of
/// re-decoding per call.  This matters hugely for the many-queries × large-subject
/// shape (e.g. TE-library query × genome DB), where the old per-call decode
/// re-cloned + re-RNG'd the whole chromosome for every query (and, under
/// `--mt-mode 1`, held one copy per worker).  The decode is deterministic
/// (`NcbiRandom` seeded only by `dna_size`), so caching is bit-for-bit identical
/// to re-decoding.
pub struct SubjectDb {
    backend: SubjectBackend,
    name_to_index: HashMap<String, usize>,
    decoded: Vec<Arc<[u8]>>,
    n_masks: Vec<Arc<[u8]>>,
}

impl SubjectDb {
    /// Open a database file, autodetecting format by extension then magic byte.
    pub fn open(path: &str) -> Result<Self> {
        let backend = Self::open_backend(path)?;
        Self::from_backend(backend)
    }

    fn open_backend(path: &str) -> Result<SubjectBackend> {
        if Self::path_looks_like_fasta(path) {
            let db = FastaDb::open(path)
                .map_err(|e| anyhow::anyhow!("opening FASTA database '{}': {}", path, e))?;
            return Ok(SubjectBackend::Fasta(db));
        }
        if Self::path_looks_like_twobit(path) {
            return Self::open_twobit(path);
        }
        // Unknown extension — peek at the first byte to decide.
        let first_byte = std::fs::File::open(path)
            .and_then(|mut f| {
                let mut buf = [0u8; 1];
                f.read_exact(&mut buf).map(|_| buf[0])
            })
            .map_err(|e| anyhow::anyhow!("cannot read database '{}': {}", path, e))?;
        if first_byte == b'>' {
            let db = FastaDb::open(path)
                .map_err(|e| anyhow::anyhow!("opening FASTA database '{}': {}", path, e))?;
            Ok(SubjectBackend::Fasta(db))
        } else {
            Self::open_twobit(path)
        }
    }

    /// Decode every subject once and cache it (BLASTNA + n_mask) as `Arc<[u8]>`.
    fn from_backend(backend: SubjectBackend) -> Result<Self> {
        let n = backend.sequences().len();
        let mut name_to_index = HashMap::with_capacity(n);
        let mut decoded: Vec<Arc<[u8]>> = Vec::with_capacity(n);
        let mut n_masks: Vec<Arc<[u8]>> = Vec::with_capacity(n);
        // Collect names first so we don't hold a borrow of backend.sequences()
        // across the decode calls.
        let names: Vec<String> = backend.sequences().iter().map(|s| s.name.clone()).collect();
        for (i, name) in names.iter().enumerate() {
            name_to_index.insert(name.clone(), i);
            decoded.push(Arc::from(backend.decode(name)?));
            n_masks.push(Arc::from(backend.nmask(name)));
        }
        Ok(SubjectDb { backend, name_to_index, decoded, n_masks })
    }

    fn path_looks_like_fasta(path: &str) -> bool {
        let lower = path.to_ascii_lowercase();
        lower.ends_with(".fa")
            || lower.ends_with(".fasta")
            || lower.ends_with(".fas")
            || lower.ends_with(".fna")
            || lower.ends_with(".fn")
    }

    fn path_looks_like_twobit(path: &str) -> bool {
        path.to_ascii_lowercase().ends_with(".2bit")
    }

    fn open_twobit(path: &str) -> Result<SubjectBackend> {
        eprintln!(
            "Note: .2bit database — non-N IUB ambiguity codes (R/Y/K/M/W/S/B/D/H/V) \
             are not preserved in the .2bit format and will be treated as N. \
             For full IUB fidelity supply the database in FASTA format."
        );
        let db = TwoBitFile::open(path)
            .map_err(|e| anyhow::anyhow!("opening 2bit database '{}': {}", path, e))?;
        Ok(SubjectBackend::TwoBit(db))
    }

    /// Metadata for all sequences in the database.
    pub fn sequences(&self) -> &[SeqInfo] {
        self.backend.sequences()
    }

    /// Cheap `Arc` clone of the cached, already-decoded BLASTNA sequence
    /// (sentinels included, ambiguities resolved).
    pub fn get_full_sequence_blastna(&self, name: &str) -> Result<Arc<[u8]>> {
        match self.name_to_index.get(name) {
            Some(&i) => Ok(self.decoded[i].clone()),
            None => Err(anyhow::anyhow!("sequence '{}' not found in database", name)),
        }
    }

    /// Cheap `Arc` clone of the cached ambiguity mask (empty if no ambiguity).
    pub fn get_n_mask(&self, name: &str) -> Arc<[u8]> {
        match self.name_to_index.get(name) {
            Some(&i) => self.n_masks[i].clone(),
            None => Arc::from(Vec::new()),
        }
    }
}
