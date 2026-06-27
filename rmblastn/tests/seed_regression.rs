//! CLI-level seed regression tests.
//!
//! The authoritative seed regression tests now live in the library:
//!   rmblast-lib/tests/dedup_seeds_regression.rs  — calls library functions directly (no CLI instrumentation)
//!
//! The former seed_regression_chr22 test (pre-dedup, 734,886 seeds via --dump-seeds) is retired:
//! the --dump-seeds flag now emits post-dedup seeds and the reference is stale.
//! The dedup_seed_regression_chr22 test is superseded by the library test above.
