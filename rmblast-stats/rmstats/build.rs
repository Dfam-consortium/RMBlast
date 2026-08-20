fn main() {
    #[cfg(feature = "alp-fit")]
    build_alp();
}

#[cfg(feature = "alp-fit")]
fn build_alp() {
    use std::path::PathBuf;

    // Location of the ALP library sources (sls_*.cpp / njn_*.cpp). Override
    // with ALP_SRC_DIR; defaults to the vendored copy in this project tree.
    let default_dir = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../ALP_1.98_LIB/cpp"
    );
    let alp_dir = std::env::var("ALP_SRC_DIR").unwrap_or_else(|_| default_dir.to_string());
    let alp_dir = PathBuf::from(alp_dir);
    println!("cargo:rerun-if-env-changed=ALP_SRC_DIR");
    println!("cargo:rerun-if-changed=csrc/alp_shim.cpp");

    let sources = [
        "njn_dynprogprob.cpp",
        "njn_dynprogproblim.cpp",
        "njn_dynprogprobproto.cpp",
        "njn_ioutil.cpp",
        "njn_localmaxstat.cpp",
        "njn_localmaxstatmatrix.cpp",
        "njn_localmaxstatutil.cpp",
        "njn_random.cpp",
        "sls_alignment_evaluer.cpp",
        "sls_alp.cpp",
        "sls_alp_data.cpp",
        "sls_alp_regression.cpp",
        "sls_alp_sim.cpp",
        "sls_basic.cpp",
        "sls_pvalues.cpp",
    ];

    let mut build = cc::Build::new();
    build.cpp(true).opt_level(3).include(&alp_dir);
    for s in &sources {
        let p = alp_dir.join(s);
        assert!(
            p.exists(),
            "ALP source {} not found; set ALP_SRC_DIR to the ALP_1.98_LIB/cpp directory",
            p.display()
        );
        build.file(p);
    }
    build.file("csrc/alp_shim.cpp");
    build.compile("alp_shim");
    // The C++ standard library the shim needs is named differently per target:
    // libstdc++ on GNU/Linux, libc++ on Apple platforms (macOS SDKs no longer
    // ship libstdc++ at all, so hardcoding stdc++ fails the link there).
    let cpp_stdlib = match std::env::var("CARGO_CFG_TARGET_OS").as_deref() {
        Ok("macos") | Ok("ios") | Ok("freebsd") | Ok("openbsd") => "c++",
        _ => "stdc++",
    };
    println!("cargo:rustc-link-lib={cpp_stdlib}");
}
