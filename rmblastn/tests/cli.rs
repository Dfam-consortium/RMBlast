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
//!
//!   `expected_defaults` is the same but with only `-gapopen 20 -gapextend 5
//!   -num_threads 1 -outfmt "6 <RM_FIELDS>"`: every other option at NCBI's default.

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
    rmblastn_env(args, &[])
}

fn rmblastn_env(args: &[&str], env: &[(&str, &str)]) -> Run {
    let mut child = Command::new(env!("CARGO_BIN_EXE_rmblastn"))
        .args(args)
        .envs(env.iter().copied())
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

/// With only the (required) gap costs given, every other option takes NCBI's
/// default.  3.0.12 defaulted --word-size to 8 where NCBI uses 11.
#[test]
fn defaults_match_ncbi() {
    let (db, q, m) = (data("db.fa"), data("query.fa"), data("comparison.matrix"));
    let r = rmblastn(&[
        "--db", &db, "--query", &q, "--matrix", &m,
        "--gapopen", "20", "--gapextend", "5", "--outfmt", RM_FIELDS,
    ]);
    assert_same_output(&r, &expected("expected_defaults.ncbi-2.17.1.tsv"));
}

// ── errors ──────────────────────────────────────────────────────────────────

/// Each bad option value is a usage error that names the option, and nothing
/// is searched.  3.0.12 panicked (word size 0-2), hung (word size 3, negative
/// or huge gap costs) or silently substituted a value for most of these.
#[test]
fn bad_option_values_are_rejected() {
    let (db, q, m) = (data("db.fa"), data("query.fa"), data("comparison.matrix"));
    let cases: &[(&[&str], &str)] = &[
        (&["--word-size", "0"], "--word-size"),
        (&["--word-size", "3"], "--word-size"),
        (&["--word-size", "4.5"], "--word-size"),
        (&["--gapopen", "-60"], "--gapopen"),
        (&["--gapopen", "2000000000"], "--gapopen"),
        (&["--gapextend", "0"], "--gapextend"),
        (&["--gapextend", "-1"], "--gapextend"),
        (&["--num-threads", "0"], "--num-threads"),
        (&["--mt-mode", "2"], "--mt-mode"),
        (&["--seed-mode", "junk"], "--seed-mode"),
        (&["--evalue", "0"], "--evalue"),
        (&["--evalue", "-1"], "--evalue"),
        (&["--evalue", "nan"], "--evalue"),
        (&["--dust", "junk"], "--dust"),
        (&["--dust", "YES"], "--dust"),
        (&["--dust", "20 64 1"], "level window linker"),
        (&["--xdrop-gap", "-1"], "--xdrop-gap"),
        (&["--min-raw-gapped-score", "-5"], "--min-raw-gapped-score"),
        (&["--matrix-lambda", "-1"], "--matrix-lambda"),
        (&["-num_threads", "0"], "--num-threads"),
    ];
    for (bad, needle) in cases {
        // the bad value last, so it overrides nothing and duplicates nothing
        let mut args = vec!["--db", &db, "--query", &q, "--matrix", &m, "--outfmt", "6 score"];
        if !bad[0].contains("gap") || bad[0].contains("xdrop") {
            args.extend(["--gapopen", "20", "--gapextend", "5"]);
        } else if bad[0] == "--gapopen" {
            args.extend(["--gapextend", "5"]);
        } else {
            args.extend(["--gapopen", "20"]);
        }
        args.extend_from_slice(bad);
        let r = rmblastn(&args);
        assert!(matches!(r.code, Some(1) | Some(2)), "{bad:?}: exit {:?}, stderr:\n{}", r.code, r.stderr);
        assert!(r.stderr.contains(needle), "{bad:?}: stderr lacks {needle:?}:\n{}", r.stderr);
        assert!(r.stdout.is_empty(), "{bad:?}: printed output");
    }
}

/// Gap costs have no usable default (NCBI rmblastn's 0/0 is rejected).
#[test]
fn gap_costs_are_required() {
    let (db, q, m) = (data("db.fa"), data("query.fa"), data("comparison.matrix"));
    let r = rmblastn(&["--db", &db, "--query", &q, "--matrix", &m, "--gapopen", "20"]);
    assert_eq!(r.code, Some(2));
    assert!(r.stderr.contains("--gapextend"), "stderr:\n{}", r.stderr);
}

/// The edge values that are allowed still run.
#[test]
fn boundary_option_values_are_accepted() {
    let (db, q, m) = (data("db.fa"), data("query.fa"), data("comparison.matrix"));
    for extra in [
        &["--gapopen", "0", "--gapextend", "1"][..],
        &["--gapopen", "20", "--gapextend", "5", "--word-size", "4"][..],
        &["--gapopen", "1000000", "--gapextend", "1000000"][..],
        &["--gapopen", "20", "--gapextend", "5", "--dust", "false", "--evalue", "1e-5",
          "-matrix_beta", "-4.9", "-mask_level", "-1"][..],
    ] {
        let mut args = vec!["--db", &db, "--query", &q, "--matrix", &m, "--outfmt", "6 score"];
        args.extend_from_slice(extra);
        let r = rmblastn(&args);
        assert_eq!(r.code, Some(0), "{extra:?}: stderr:\n{}", r.stderr);
    }
}

#[test]
fn missing_database_is_an_error() {
    let (q, m) = (data("query.fa"), data("comparison.matrix"));
    let r = rmblastn(&["--db", "/nonexistent/db", "--query", &q, "--matrix", &m,
                       "--gapopen", "20", "--gapextend", "5"]);
    assert_ne!(r.code, Some(0));
    assert!(r.stderr.contains("database not found"), "stderr: {}", r.stderr);
    assert!(r.stdout.is_empty());
}

#[test]
fn missing_matrix_is_an_error() {
    let (db, q) = (data("db.fa"), data("query.fa"));
    let r = rmblastn(&["--db", &db, "--query", &q, "--matrix", "no-such.matrix",
                       "--gapopen", "20", "--gapextend", "5"]);
    assert_ne!(r.code, Some(0));
    assert!(r.stderr.contains("Cannot find substitution matrix"), "stderr: {}", r.stderr);
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

// ── mt-mode 1 failure handling (needs `--features diagnostics`) ─────────────

/// A query file with `n` records: the three fixture queries, repeated.
#[cfg(feature = "diagnostics")]
fn repeated_queries(n: usize) -> tempfile_path::TempFile {
    let text = std::fs::read_to_string(data("query.fa")).unwrap();
    let records: Vec<&str> = text.split('>').filter(|r| !r.is_empty()).collect();
    let mut out = String::new();
    for i in 0..n {
        // ">q1 desc\nACGT..." -> ">q1_<i> desc\nACGT...": unique ids, same sequences
        let (id, rest) = records[i % records.len()].split_once(' ').unwrap();
        out.push_str(&format!(">{id}_{i} {rest}"));
    }
    tempfile_path::TempFile::with_contents(&out)
}

#[cfg(feature = "diagnostics")]
mod tempfile_path {
    pub struct TempFile(pub std::path::PathBuf);
    impl TempFile {
        pub fn with_contents(s: &str) -> Self {
            let p = std::env::temp_dir().join(format!(
                "rmblastn-cli-test-{}-{:?}.fa", std::process::id(), std::thread::current().id()));
            std::fs::write(&p, s).unwrap();
            TempFile(p)
        }
        pub fn path(&self) -> &str { self.0.to_str().unwrap() }
    }
    impl Drop for TempFile {
        fn drop(&mut self) { std::fs::remove_file(&self.0).ok(); }
    }
}

/// A query search that panics in `--mt-mode 1` must end the run with an error,
/// whether the panic comes before the writer's window fills or after.  3.0.12
/// hung forever: the writer waited for the panicked query's results.
#[cfg(feature = "diagnostics")]
fn mt_mode_1_panic_exits(n_queries: usize, panic_at: usize) {
    let q = repeated_queries(n_queries);
    let (db, m) = (data("db.fa"), data("comparison.matrix"));
    let at = panic_at.to_string();
    let r = rmblastn_env(
        &["--db", &db, "--query", q.path(), "--matrix", &m,
          "--gapopen", "20", "--gapextend", "5",
          "--mt-mode", "1", "--num-threads", "2", "--outfmt", RM_FIELDS],
        &[("RMBLAST_TEST_PANIC_AT_QUERY", &at)],
    );
    assert_eq!(r.code, Some(1), "stderr:\n{}", r.stderr);
    assert!(r.stderr.contains(&format!("failed: RMBLAST_TEST_PANIC_AT_QUERY={panic_at}")),
        "stderr:\n{}", r.stderr);
}

#[cfg(feature = "diagnostics")]
#[test]
fn mt_mode_1_panic_past_window_exits() {
    // 2 threads -> window 16: query 30 is handed out only after the writer
    // has written 14 queries, so a hang would leave the dispatcher parked.
    mt_mode_1_panic_exits(60, 30);
}

#[cfg(feature = "diagnostics")]
#[test]
fn mt_mode_1_panic_within_window_exits() {
    mt_mode_1_panic_exits(3, 0);
}
