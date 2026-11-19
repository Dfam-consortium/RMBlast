//! CLI tests for the dustmasker front end.
//!
//! `tests/data/deflines.fa` holds 306 records whose deflines were generated at
//! random from spaces, tabs, punctuation, brackets, HTML entities, TPA/TSA/MAG
//! prefixes and non-ASCII text, plus edge cases ("...", "TPA: ,").  The
//! `expected_deflines.<fmt>.ncbi-2.17.0` files are NCBI dustmasker 2.17.0
//! output for `-in deflines.fa -outfmt <fmt>`.

use std::path::PathBuf;
use std::process::{Command, Output};

fn data(name: &str) -> String {
    let p: PathBuf = [env!("CARGO_MANIFEST_DIR"), "tests", "data", name].iter().collect();
    p.to_str().unwrap().to_owned()
}

fn dustmasker(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_dustmasker")).args(args).output().expect("run dustmasker")
}

/// A temp file removed on drop.
struct TempFile(PathBuf);
impl TempFile {
    fn new(tag: &str, contents: &[u8]) -> Self {
        let p = std::env::temp_dir().join(format!("dustmasker-cli-{}-{tag}.fa", std::process::id()));
        std::fs::write(&p, contents).unwrap();
        TempFile(p)
    }
    fn path(&self) -> &str {
        self.0.to_str().unwrap()
    }
}
impl Drop for TempFile {
    fn drop(&mut self) {
        std::fs::remove_file(&self.0).ok();
    }
}

/// Deflines print exactly as NCBI prints them (spaces collapsed, cut at a
/// tab, trailing periods stripped, entities decoded ...).  3.0.12 echoed them
/// verbatim.
#[test]
fn deflines_match_ncbi_in_every_format() {
    for fmt in ["interval", "acclist", "fasta"] {
        let r = dustmasker(&["-in", &data("deflines.fa"), "-outfmt", fmt]);
        assert_eq!(r.status.code(), Some(0), "{fmt}: {}", String::from_utf8_lossy(&r.stderr));
        let want = std::fs::read(data(&format!("expected_deflines.{fmt}.ncbi-2.17.0"))).unwrap();
        if r.stdout != want {
            let first = r.stdout.split(|&b| b == b'\n').zip(want.split(|&b| b == b'\n'))
                .position(|(a, b)| a != b);
            panic!("{fmt}: output differs from NCBI, first at line {:?}", first.map(|i| i + 1));
        }
    }
}

/// An empty record is an error (NCBI: "FASTA-Reader: No residues given");
/// 3.0.12 skipped it and exited 0.
#[test]
fn empty_record_is_an_error() {
    let f = TempFile::new("empty", b">a\nACGT\n>b\n>c\nACGT\n");
    let r = dustmasker(&["-in", f.path()]);
    assert_eq!(r.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&r.stderr).contains("record 'b' has no sequence"));
}

/// A character that is not a nucleotide code is an error; NCBI drops it and
/// 3.0.12 counted it as N, shifting the intervals after it.
#[test]
fn bad_sequence_character_is_an_error() {
    let f = TempFile::new("bad", b">a\nACGT*ACGT\n");
    let r = dustmasker(&["-in", f.path()]);
    assert_eq!(r.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&r.stderr).contains("'*' is not a nucleotide code"));
}

/// Usage errors exit 1, as NCBI's do (clap's default is 2); --help exits 0.
#[test]
fn usage_errors_exit_1() {
    for args in [&["-window", "-1"][..], &["-linker", "x"][..], &["-no_such_option"][..]] {
        let r = dustmasker(args);
        assert_eq!(r.status.code(), Some(1), "{args:?}");
    }
    assert_eq!(dustmasker(&["-help"]).status.code(), Some(0));
}
