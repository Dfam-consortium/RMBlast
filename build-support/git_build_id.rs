// Shared by rmblastn/build.rs and dustmasker/build.rs (include!d, so it is
// compiled into each build script).
//
// Sets RMBLAST_VERSION to "<package version> (<commit>)", or
// "<package version> (<commit>-dirty)" when a file that builds the binaries
// has uncommitted changes, so `--version` says which code a binary was built
// from.  3.0.11's target/release once reported 3.0.11 while holding older
// code; the commit in the string makes that visible.  Outside a git checkout
// (a source tarball) it is the package version alone.

use std::path::{Path, PathBuf};
use std::process::Command;

/// Paths, relative to the workspace root, whose content goes into the binaries.
const BUILD_INPUTS: &[&str] = &[
    "Cargo.toml", "Cargo.lock", "rmblast-lib", "rmblastn", "dustmasker", "rmblast-stats",
    "build-support",
];

fn git(root: &Path, args: &[&str]) -> Option<String> {
    let out = Command::new("git").arg("-C").arg(root).args(args).output().ok()?;
    out.status.success().then(|| String::from_utf8_lossy(&out.stdout).trim().to_string())
}

fn emit_build_id() {
    let manifest = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR"));
    let root = manifest.parent().expect("crate inside the workspace").to_path_buf();
    let pkg = std::env::var("CARGO_PKG_VERSION").expect("CARGO_PKG_VERSION");

    for p in BUILD_INPUTS {
        println!("cargo:rerun-if-changed={}", root.join(p).display());
    }
    let version = match git(&root, &["rev-parse", "--short=7", "HEAD"]) {
        None => pkg,
        Some(commit) => {
            // Rerun when HEAD moves: a commit, checkout or reset.
            if let Some(dir) = git(&root, &["rev-parse", "--absolute-git-dir"]) {
                let dir = PathBuf::from(dir);
                println!("cargo:rerun-if-changed={}", dir.join("HEAD").display());
                println!("cargo:rerun-if-changed={}", dir.join("packed-refs").display());
                if let Some(r) = git(&root, &["symbolic-ref", "-q", "HEAD"]) {
                    println!("cargo:rerun-if-changed={}", dir.join(r).display());
                }
            }
            let mut args = vec!["status", "--porcelain", "--"];
            args.extend_from_slice(BUILD_INPUTS);
            let dirty = git(&root, &args).map(|s| !s.is_empty()).unwrap_or(false);
            format!("{pkg} ({commit}{})", if dirty { "-dirty" } else { "" })
        }
    };
    println!("cargo:rustc-env=RMBLAST_VERSION={version}");
}
