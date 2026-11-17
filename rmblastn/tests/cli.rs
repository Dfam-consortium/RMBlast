//! CLI-level tests: run the built `rmblastn` binary on small fixtures and check
//! stdout, stderr and the exit status.  Every run has a time limit, so a hang
//! fails the test instead of stalling `cargo test`.
//!
//! Fixtures (`tests/data/`):
//! - `query.fa`: q1 and q2 are 8 kb and 4 kb of `human-1mb.fa` (368001-376000
//!   and 72001-76000) with hits on both strands; q3 is 1 kb (22001-23000) with
//!   no hits.
//! - `db.fa`: the 17 `shortlib` records that hit q1 or q2.
//! - `comparison.matrix`: NCBI's RMBlast `comparison.matrix`.
//! - `expected_*.ncbi-2.17.1.tsv`: output of NCBI RMBlast 2.17.1 on these
//!   fixtures, made with `makeblastdb -in db.fa -dbtype nucl` and
//!
//!   ```text
//!   BLASTMAT=tests/data rmblastn -num_alignments 9999999 -db db -query query.fa \
//!     -gapopen 20 -gapextend 5 -word_size 14 -xdrop_ungap 400 -xdrop_gap_final 200 \
//!     -xdrop_gap 100 -min_raw_gapped_score 200 -matrix comparison.matrix \
//!     -num_threads 1 -outfmt "6 <RM_FIELDS>" \
//!     rm_defaults:    -mask_level 80 -complexity_adjust -dust yes
//!     nomask_nodust:  -mask_level 101 -dust no
//!   ```

use std::io::Read;
use std::path::PathBuf;
use std::process::{Command, Stdio};
use std::time::{Duration, Instant};

const RM_FIELDS: &str = "6 score perc_sub perc_query_gap perc_db_gap qseqid qstart qend qlen \
    sstrand sseqid sstart send slen kdiv cpg_kdiv transi transv cpg_sites qseq sseq";

const TIME_LIMIT: Duration = Duration::from_secs(60);

fn data(name: &str) -> String {
    let p: PathBuf = [env!("CARGO_MANIFEST_DIR"), "tests", "data", name].iter().collect();
    p.to_str().unwrap().to_owned()
}

struct Run {
    code: Option<i32>,
    stdout: String,
    stderr: String,
}

/// Run rmblastn with `args`; panic if it is still running after TIME_LIMIT.
fn rmblastn(args: &[&str]) -> Run {
    let mut child = Command::new(env!("CARGO_BIN_EXE_rmblastn"))
        .args(args)
        .env_remove("BLASTMAT")
        .stdin(Stdio::null())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .expect("spawn rmblastn");
    // Drain both pipes on their own threads so a large output can't block the child.
    let drain = |mut r: Box<dyn Read + Send>| {
        std::thread::spawn(move || {
            let mut v = Vec::new();
            r.read_to_end(&mut v).unwrap();
            String::from_utf8_lossy(&v).into_owned()
        })
    };
    let out = drain(Box::new(child.stdout.take().unwrap()));
    let err = drain(Box::new(child.stderr.take().unwrap()));
    let start = Instant::now();
    let status = loop {
        if let Some(s) = child.try_wait().unwrap() {
            break s;
        }
        if start.elapsed() > TIME_LIMIT {
            child.kill().ok();
            child.wait().ok();
            panic!("rmblastn {:?} still running after {:?}: killed", args, TIME_LIMIT);
        }
        std::thread::sleep(Duration::from_millis(20));
    };
    Run { code: status.code(), stdout: out.join().unwrap(), stderr: err.join().unwrap() }
}

/// The search parameters shared by every fixture run (RepeatMasker's).
fn search_args(db: &str, query: &str, matrix: &str) -> Vec<String> {
    [
        "--db", db, "--query", query, "--matrix", matrix,
        "--gapopen", "20", "--gapextend", "5", "--word-size", "14",
        "--xdrop-ungap", "400", "--xdrop-gap-final", "200", "--xdrop-gap", "100",
        "--min-raw-gapped-score", "200", "--outfmt", RM_FIELDS,
    ]
    .iter()
    .map(|s| s.to_string())
    .collect()
}

fn run_search(extra: &[&str]) -> Run {
    let (db, q, m) = (data("db.fa"), data("query.fa"), data("comparison.matrix"));
    let mut args = search_args(&db, &q, &m);
    args.extend(extra.iter().map(|s| s.to_string()));
    let refs: Vec<&str> = args.iter().map(|s| s.as_str()).collect();
    rmblastn(&refs)
}

fn expected(name: &str) -> String {
    std::fs::read_to_string(data(name)).unwrap_or_else(|e| panic!("{name}: {e}"))
}

/// Compare line by line so a failure names the first differing line.
fn assert_same_output(got: &Run, want: &str) {
    assert_eq!(got.code, Some(0), "exit status; stderr:\n{}", got.stderr);
    if got.stdout != want {
        let (g, w): (Vec<_>, Vec<_>) = (got.stdout.lines().collect(), want.lines().collect());
        let first = g.iter().zip(&w).position(|(a, b)| a != b).unwrap_or(g.len().min(w.len()));
        panic!(
            "output differs from NCBI 2.17.1 ({} lines vs {} expected); first difference at line {}:\n got: {:?}\nwant: {:?}",
            g.len(), w.len(), first + 1, g.get(first), w.get(first)
        );
    }
}

// ── output matches NCBI ─────────────────────────────────────────────────────

#[test]
fn rm_defaults_match_ncbi() {
    let r = run_search(&["--mask-level", "80", "--complexity-adjust", "--dust", "yes", "--num-threads", "1"]);
    assert_same_output(&r, &expected("expected_rm_defaults.ncbi-2.17.1.tsv"));
}

#[test]
fn rm_defaults_match_ncbi_threaded() {
    let r = run_search(&["--mask-level", "80", "--complexity-adjust", "--dust", "yes", "--num-threads", "4"]);
    assert_same_output(&r, &expected("expected_rm_defaults.ncbi-2.17.1.tsv"));
}

#[test]
fn rm_defaults_match_ncbi_mt_mode_1() {
    let r = run_search(&[
        "--mask-level", "80", "--complexity-adjust", "--dust", "yes",
        "--num-threads", "4", "--mt-mode", "1",
    ]);
    assert_same_output(&r, &expected("expected_rm_defaults.ncbi-2.17.1.tsv"));
}

#[test]
fn no_mask_level_no_dust_match_ncbi() {
    let r = run_search(&["--mask-level", "101", "--dust", "no", "--num-threads", "4"]);
    assert_same_output(&r, &expected("expected_nomask_nodust.ncbi-2.17.1.tsv"));
}

/// NCBI's single-dash, underscore spellings give the same result as ours.
#[test]
fn ncbi_spelling_matches_native_spelling() {
    let (db, q, m) = (data("db.fa"), data("query.fa"), data("comparison.matrix"));
    let r = rmblastn(&[
        "-db", &db, "-query", &q, "-matrix", &m,
        "-gapopen", "20", "-gapextend", "5", "-word_size", "14",
        "-xdrop_ungap", "400", "-xdrop_gap_final", "200", "-xdrop_gap", "100",
        "-min_raw_gapped_score", "200", "-outfmt", RM_FIELDS,
        "-mask_level", "80", "-complexity_adjust", "-dust", "yes", "-num_threads", "2",
    ]);
    assert_same_output(&r, &expected("expected_rm_defaults.ncbi-2.17.1.tsv"));
}

// ── errors ──────────────────────────────────────────────────────────────────

#[test]
fn missing_database_is_an_error() {
    let (q, m) = (data("query.fa"), data("comparison.matrix"));
    let r = rmblastn(&["--db", "/nonexistent/db", "--query", &q, "--matrix", &m]);
    assert_ne!(r.code, Some(0));
    assert!(r.stderr.contains("database not found"), "stderr: {}", r.stderr);
    assert!(r.stdout.is_empty());
}

#[test]
fn missing_matrix_is_an_error() {
    let (db, q) = (data("db.fa"), data("query.fa"));
    let r = rmblastn(&["--db", &db, "--query", &q, "--matrix", "no-such.matrix"]);
    assert_ne!(r.code, Some(0));
    assert!(!r.stderr.is_empty());
    assert!(r.stdout.is_empty());
}

#[test]
fn unknown_option_is_rejected() {
    let r = rmblastn(&["--no-such-option"]);
    assert_ne!(r.code, Some(0));
    assert!(r.stderr.contains("--no-such-option"), "stderr: {}", r.stderr);
}

#[test]
fn version_prints_crate_version() {
    let r = rmblastn(&["--version"]);
    assert_eq!(r.code, Some(0));
    assert_eq!(r.stdout.trim(), format!("rmblastn {}", env!("CARGO_PKG_VERSION")));
}
