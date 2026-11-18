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
//!   `expected_query_dash` / `expected_db_dash` are rm_defaults with
//!   `query_dash.fa` (q1 with '-' at 8 positions) against `db.fa`, and
//!   `query.fa` against `db_dash.fa` (4 records with '-' at 4 positions).
//!
//!   `expected_pairwise_{defaults,rm}` are NCBI's pairwise output (no -outfmt)
//!   from a database built in this directory with `makeblastdb -in db.fa
//!   -dbtype nucl -parse_seqids -out db`, so its title is "db.fa"; the
//!   "Posted date" line is removed.  defaults: `-gapopen 20 -gapextend 5
//!   -num_threads 1`; rm: the rm_defaults options plus -xdrop_* as above,
//!   with dust on (the default).
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
    rmblastn_in(args, env, None)
}

fn rmblastn_in(args: &[&str], env: &[(&str, &str)], dir: Option<&str>) -> Run {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_rmblastn"));
    if let Some(d) = dir {
        cmd.current_dir(d);
    }
    let mut child = cmd
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

/// A '-' in a sequence is N, on either side (3.0.12 encoded it as the gap
/// sentinel: `GGC--GGG` instead of NCBI's `GGCNNGGG`, and lower scores).
#[test]
fn dash_in_sequence_is_n_like_ncbi() {
    let m = data("comparison.matrix");
    for (query, db, want) in [
        ("query_dash.fa", "db.fa", "expected_query_dash.ncbi-2.17.1.tsv"),
        ("query.fa", "db_dash.fa", "expected_db_dash.ncbi-2.17.1.tsv"),
    ] {
        let (q, d) = (data(query), data(db));
        let mut args = search_args(&d, &q, &m);
        args.extend(["--mask-level", "80", "--complexity-adjust", "--dust", "yes", "--num-threads", "1"]
            .iter().map(|s| s.to_string()));
        let refs: Vec<&str> = args.iter().map(|s| s.as_str()).collect();
        assert_same_output(&rmblastn(&refs), &expected(want));
    }
}

/// Pairwise output is NCBI's, line for line, apart from "Posted date": the
/// version line, thousands separators, the full subject defline, DUST-masked
/// query bases in lower case, the no-hits block and the blank lines between
/// queries all differed in 3.0.12.  Run from the fixture directory with
/// --db db.fa so the database title matches NCBI's ("db.fa").
#[test]
fn pairwise_matches_ncbi() {
    let dir = data("");
    let rm: &[&str] = &["--word-size", "14", "--mask-level", "80", "--complexity-adjust",
        "--min-raw-gapped-score", "200", "--xdrop-ungap", "400", "--xdrop-gap-final", "200",
        "--xdrop-gap", "100"];
    for (extra, want) in [(&[][..], "expected_pairwise_defaults.ncbi-2.17.1.txt"),
                          (rm, "expected_pairwise_rm.ncbi-2.17.1.txt")] {
        let mut args = vec!["--db", "db.fa", "--query", "query.fa", "--matrix", "comparison.matrix",
                            "--gapopen", "20", "--gapextend", "5", "--num-threads", "2"];
        args.extend_from_slice(extra);
        let mut r = rmblastn_in(&args, &[], Some(&dir));
        r.stdout = r.stdout.lines().filter(|l| !l.starts_with("    Posted date:"))
            .map(|l| format!("{l}\n")).collect();
        assert_same_output(&r, &expected(want));
    }
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

/// When E-values come from the built-in table, stderr says the entry was
/// chosen by the matrix file name; without statistics nothing is printed.
#[test]
fn table_statistics_say_they_are_chosen_by_name() {
    let (db, q, m) = (data("db.fa"), data("query.fa"), data("comparison.matrix"));
    let base = ["--db", &db, "--query", &q, "--matrix", &m, "--gapopen", "20", "--gapextend", "5"];
    let with = rmblastn(&[&base[..], &["--outfmt", "6 score evalue"]].concat());
    assert_eq!(with.code, Some(0));
    assert!(with.stderr.contains(
        "parameters for 'comparison.matrix' at gap costs 20/5 taken from the built-in table by matrix file name"),
        "stderr:\n{}", with.stderr);
    let without = rmblastn(&[&base[..], &["--outfmt", "6 score"]].concat());
    assert!(!without.stderr.contains("built-in table"), "stderr:\n{}", without.stderr);
}

// ── errors ──────────────────────────────────────────────────────────────────

/// A database with two records of the same id is rejected.  3.0.12 searched
/// the last one twice and the first never: with an Alu first and random
/// sequence second it printed nothing, where NCBI reports the Alu hit.
#[test]
fn duplicate_subject_ids_are_rejected() {
    let text = std::fs::read_to_string(data("db.fa")).unwrap();
    let first = &text[..text[1..].find('>').unwrap() + 1]; // first record
    let db = tempfile_path::TempFile::with_contents(&format!("{text}{first}"));
    let (q, m) = (data("query.fa"), data("comparison.matrix"));
    let r = rmblastn(&["--db", db.path(), "--query", &q, "--matrix", &m,
                       "--gapopen", "20", "--gapextend", "5", "--outfmt", "6 score"]);
    assert_eq!(r.code, Some(1), "stderr:\n{}", r.stderr);
    assert!(r.stderr.contains("duplicate sequence id 'AluJb_short_#SINE/Alu' (records 1 and 18)"),
        "stderr:\n{}", r.stderr);
    assert!(r.stdout.is_empty());
}

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

/// A malformed FASTA query or database is an error naming the line.  3.0.12
/// counted '*' as N in a database, skipped sequence before the first '>',
/// and so on.
#[test]
fn malformed_fasta_is_rejected() {
    let bad = tempfile_path::TempFile::with_contents(">x\nACGT\nAC*GT\n");
    let (db, q, m) = (data("db.fa"), data("query.fa"), data("comparison.matrix"));
    for (db, q) in [(db.as_str(), bad.path()), (bad.path(), q.as_str())] {
        let r = rmblastn(&["--db", db, "--query", q, "--matrix", &m, "--gapopen", "20",
                           "--gapextend", "5", "--outfmt", "6 score"]);
        assert_eq!(r.code, Some(1), "stderr:\n{}", r.stderr);
        assert!(r.stderr.contains("FASTA line 3: '*' is not a nucleotide code (in record 'x')"),
            "stderr:\n{}", r.stderr);
    }
}

/// An empty database is an error, FASTA or 2bit (3.0.12 exited 0 silently).
#[test]
fn empty_database_is_an_error() {
    let mut zero_2bit = Vec::new();
    for w in [0x1A41_2743u32, 0, 0, 0] {
        zero_2bit.extend_from_slice(&w.to_le_bytes()); // magic, version, seqCount 0, reserved
    }
    let (q, m) = (data("query.fa"), data("comparison.matrix"));
    for db in [tempfile_path::TempFile::with_contents(""), tempfile_path::TempFile::with_bytes(&zero_2bit, "2bit")] {
        let r = rmblastn(&["--db", db.path(), "--query", &q, "--matrix", &m,
                           "--gapopen", "20", "--gapextend", "5", "--outfmt", "6 score"]);
        assert_eq!(r.code, Some(1), "{}: stderr:\n{}", db.path(), r.stderr);
        assert!(r.stderr.contains("contains no sequences"), "stderr:\n{}", r.stderr);
    }
}

/// An empty query file warns as NCBI does and exits 0, in both mt-modes.
#[test]
fn empty_query_warns_like_ncbi() {
    let q = tempfile_path::TempFile::with_contents("");
    let (db, m) = (data("db.fa"), data("comparison.matrix"));
    for mode in ["0", "1"] {
        let r = rmblastn(&["--db", &db, "--query", q.path(), "--matrix", &m, "--gapopen", "20",
                           "--gapextend", "5", "--outfmt", "6 score", "--mt-mode", mode, "--num-threads", "2"]);
        assert_eq!(r.code, Some(0), "mt-mode {mode}: stderr:\n{}", r.stderr);
        assert!(r.stderr.contains("Warning: [rmblastn] Query is Empty!"), "stderr:\n{}", r.stderr);
        assert!(r.stdout.is_empty());
    }
}

/// A gilist that selects nothing is a valid search with no hits, with a warning.
#[test]
fn gilist_matching_nothing_warns() {
    let gl = tempfile_path::TempFile::with_contents("no-such-id\n");
    let (db, q, m) = (data("db.fa"), data("query.fa"), data("comparison.matrix"));
    let r = rmblastn(&["--db", &db, "--query", &q, "--matrix", &m, "--gapopen", "20",
                       "--gapextend", "5", "--outfmt", "6 score", "--gilist", gl.path()]);
    assert_eq!(r.code, Some(0), "stderr:\n{}", r.stderr);
    assert!(r.stderr.contains("matches no sequence"), "stderr:\n{}", r.stderr);
    assert!(r.stdout.is_empty());
}

/// Output formats and fields that are not implemented are refused, naming
/// the problem.  3.0.12 dropped unknown fields silently (shifting the columns),
/// treated "7 ..." and "10 ..." as 6, and refused a bare "6" only by accident.
#[test]
fn unsupported_outfmt_is_rejected() {
    let (db, q, m) = (data("db.fa"), data("query.fa"), data("comparison.matrix"));
    for (fmt, needle) in [
        ("6", "needs a field list"),
        ("7 score", "output format '7' is not supported"),
        ("10 score qseqid", "output format '10' is not supported"),
        ("6 score pident qseqid", "unknown -outfmt field 'pident'"),
        ("6 std", "unknown -outfmt field 'std'"),
        ("0 score", "takes no fields"),
    ] {
        let r = rmblastn(&["--db", &db, "--query", &q, "--matrix", &m, "--gapopen", "20",
                           "--gapextend", "5", "--outfmt", fmt]);
        assert_eq!(r.code, Some(1), "{fmt:?}: stderr:\n{}", r.stderr);
        assert!(r.stderr.contains(needle), "{fmt:?}: stderr:\n{}", r.stderr);
        assert!(r.stdout.is_empty(), "{fmt:?}");
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

mod tempfile_path {
    use std::sync::atomic::{AtomicUsize, Ordering};
    static NEXT: AtomicUsize = AtomicUsize::new(0);

    /// A file in the temp dir, removed on drop.
    pub struct TempFile(pub std::path::PathBuf);
    impl TempFile {
        pub fn with_contents(s: &str) -> Self {
            Self::with_bytes(s.as_bytes(), "fa")
        }
        pub fn with_bytes(b: &[u8], ext: &str) -> Self {
            let p = std::env::temp_dir().join(format!(
                "rmblastn-cli-test-{}-{}.{ext}", std::process::id(), NEXT.fetch_add(1, Ordering::Relaxed)));
            std::fs::write(&p, b).unwrap();
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
