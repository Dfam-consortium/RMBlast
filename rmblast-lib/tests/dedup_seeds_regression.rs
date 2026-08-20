//! Regression test: post-deduplication seeds on chr22.fa must match the NCBI-derived reference.
//!
//! Calls the library directly (no CLI instrumentation) so the test remains valid even after
//! debug flags like --dump-seeds are removed.  The golden reference was produced by NCBI
//! rmblastn 2.17.0+ with BLAST_DUMP_SEEDS=1 using the same parameters, validated against
//! the Rust implementation (diff = 0).
//!
//! Reference: tests/data/chr22_dedup_seeds.txt  (247,629 seeds)
//!
//! Run with: cargo test -p rmblast-lib --test dedup_seeds_regression -- --ignored --nocapture

use rmblast_lib::matrix::ScoreMatrix;
use rmblast_lib::options::{SearchParams, SeedMode};
use rmblast_lib::search::{build_query_lookup, search_with_query_lookup_seeds, SeedRecord};
use rmblast_lib::seq::{FastaReader, TwoBitFile};
use std::io::BufReader;
use std::path::PathBuf;

fn data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data")
}

fn workspace_test_data_dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("oxyblast-lib has a parent workspace dir")
        .join("tests/data")
}

/// Mirror the NCBI BlastSplitQuery chunking used in search_db_parallel (main.rs).
fn compute_chunks(query_len: usize) -> Vec<(usize, usize)> {
    const INITIAL_CHUNK_SIZE: usize = 1_000_000;
    const OVERLAP: usize = 100;
    let initial_step = INITIAL_CHUNK_SIZE - OVERLAP;
    let num_chunks = if INITIAL_CHUNK_SIZE > OVERLAP { query_len / initial_step } else { 0 };
    let (chunk_size, effective_num_chunks) = if num_chunks <= 1 {
        (query_len, 1usize)
    } else {
        let mut cs = (query_len + (num_chunks - 1) * OVERLAP) / num_chunks;
        if num_chunks < cs.saturating_sub(OVERLAP) { cs += 1; }
        (cs, num_chunks)
    };
    let step = chunk_size.saturating_sub(OVERLAP);
    let mut chunks = Vec::new();
    for k in 0..effective_num_chunks {
        let start = k * step;
        if start >= query_len { break; }
        let end = (start + chunk_size).min(query_len);
        chunks.push((start, end));
    }
    chunks
}

fn seed_text(s: &SeedRecord) -> String {
    let strand = if s.strand == rmblast_lib::hits::Strand::Plus { '+' } else { '-' };
    format!("SEED strand={} q={} s={}", strand, s.q_off, s.s_off)
}

#[test]
#[ignore]
fn dedup_seed_regression_chr22() {
    let ws_data = workspace_test_data_dir();
    let chr22    = ws_data.join("chr22.fa");
    let alu_2bit = ws_data.join("alu.2bit");
    let matrix_path = data_dir().join("comparison.matrix");
    let reference   = data_dir().join("chr22_dedup_seeds.txt");

    for p in &[&chr22, &alu_2bit, &matrix_path, &reference] {
        if !p.exists() {
            eprintln!("skipping: {} not found", p.display());
            return;
        }
    }

    let matrix = ScoreMatrix::from_file(matrix_path.to_str().unwrap()).expect("load matrix");

    let params = SearchParams {
        gap_open: 20,
        gap_extend: 5,
        word_size: 14,
        xdrop_ungap: 400,
        xdrop_gap: 100,
        xdrop_gap_final: 200,
        min_raw_gapped_score: 200,
        complexity_adjust: true,
        dust: false,
        mask_level: 80,
        num_threads: 1,
        seed_mode: SeedMode::Combined,
        matrix_name: matrix_path.to_string_lossy().into_owned(),
        ..Default::default()
    };

    let db = TwoBitFile::open(&alu_2bit).expect("open alu.2bit");
    let subject_names: Vec<String> = db.sequences.iter().map(|s| s.name.clone()).collect();

    let chr22_file = std::fs::File::open(&chr22).expect("open chr22.fa");
    let mut reader = FastaReader::new(BufReader::new(chr22_file));
    let qrec = reader.next_record().expect("read chr22").expect("chr22 not empty");

    let query = &qrec.seq;
    let full_q_len = query.len().saturating_sub(2);
    let chunks = compute_chunks(full_q_len);

    let mut all_seeds: Vec<String> = Vec::new();

    for (chunk_start, chunk_end) in chunks {
        let chunk_len = chunk_end - chunk_start;
        if chunk_len < params.word_size { continue; }

        let mut chunk: Vec<u8> = Vec::with_capacity(chunk_len + 2);
        chunk.push(query[0]);
        chunk.extend_from_slice(&query[1 + chunk_start..1 + chunk_end]);
        chunk.push(*query.last().unwrap_or(&14));

        let (ql, _) = build_query_lookup(&chunk, &params);

        for name in &subject_names {
            let seq = db.get_full_sequence_blastna(name).expect("get subject");
            let mut seeds: Vec<SeedRecord> = Vec::new();
            search_with_query_lookup_seeds(
                &ql, &chunk, &qrec.id, &seq, name,
                &params, &matrix, chunk_start as u32,
                &mut seeds, &[], None,
            );
            for s in &seeds {
                all_seeds.push(seed_text(s));
            }
        }
    }

    all_seeds.sort_unstable();

    let reference_text = std::fs::read_to_string(&reference).expect("read reference");
    let ref_seeds: Vec<&str> = reference_text.lines().collect();

    assert_eq!(
        all_seeds.len(), ref_seeds.len(),
        "dedup seed count: Rust={} reference={}", all_seeds.len(), ref_seeds.len()
    );

    let mut first_diff = None;
    for (i, (r, n)) in all_seeds.iter().zip(ref_seeds.iter()).enumerate() {
        if r.as_str() != *n {
            first_diff = Some((i, r.as_str(), *n));
            break;
        }
    }

    if let Some((i, r, n)) = first_diff {
        panic!("dedup seed {} differs:\n  Rust: {}\n  ref:  {}", i, r, n);
    }

    eprintln!("dedup_seed_regression_chr22: {} seeds matched", all_seeds.len());
}
