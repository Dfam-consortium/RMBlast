//! E-value / bit-score statistics context for rmblastn output.
//!
//! Resolves the gapped Karlin-Altschul parameters for a custom matrix through
//! a five-level hierarchy, then builds a per-query [`RmStats`] used to fill
//! the `evalue` / `bitscore` tabular output fields:
//!
//! 1. **Baked table** — the RMBlast `rmblast_*_values` table by matrix
//!    basename + gap costs (NCBI Mode 1).
//! 2. **CLI overrides** — `-matrix_lambda/-matrix_k/-matrix_alpha`
//!    (+ optional `-matrix_beta`), NCBI Mode 2.
//! 3. **`# KARLIN` matrix comment** — parameters declared in the matrix file
//!    itself (rmblastn-rs extension; see [`crate::matrix::KarlinComment`]).
//! 4. **ALP fit** — startup-time Gumbel fit of the matrix with the vendored
//!    ALP library (feature `alp-fit`; rmblastn-rs extension).
//! 5. **Sentinel** — NCBI Mode 3: every hit reports `evalue = 1.0`,
//!    `bit_score = 0.0`.
//!
//! Sources 3 and 4 are funneled through the Mode-2 machinery (synthesized
//! [`MatrixCliOverrides`]), so their E-value semantics are identical to a
//! user passing the same values on the command line (H = lambda/alpha, beta
//! feeds the length adjustment).
//!
//! ## Two resolution modes
//!
//! NCBI 2.17.1 does one result-affecting thing with these statistics: the
//! -evalue reap (`Blast_HSPListReapByEvalue`, applied per HSP list before
//! masklevel, default threshold 10).  Historically (pre-2.17, all-sentinel)
//! that reap was a no-op; 2.17's hardcoded ALP table silently re-enabled it
//! for the table matrices.  rmblastn-rs therefore has two modes:
//!
//! * **`NcbiCompat`** (opt-in, `--ncbi-compat`): emulate 2.17.1 exactly —
//!   statistics only from NCBI-native sources (the hardcoded table or
//!   `-matrix_*` overrides), sentinel rendering (`1.0`/` 0.0`) elsewhere, and
//!   the reap always on (caller defaults the threshold to NCBI's 10).  Two
//!   C-table warts are reproduced faithfully: `30p53g.matrix` @22/5 uses the
//!   C source's placeholder (lambda=K=H=0.1, alpha=1, beta=0) for BOTH
//!   reporting and reaping, and `comparison.matrix` @20/5 (absent from the C
//!   table) is sentinel.  Use this for 2.17.1-parity benchmarks.
//! * **`Default`**: the full five-level hierarchy above supplies the printed
//!   columns (real fits for 30p53g and comparison.matrix), and NOTHING is
//!   culled unless the user explicitly passes `--evalue`; when they do, the
//!   reap uses the same statistics that are reported.
//!
//! In both modes the ungapped seeding cutoff (`search::ka_cutoff`) keeps its
//! own NCBI-verbatim table.
//!
//! **Laziness:** the table/CLI/comment sources are cheap arithmetic.  The
//! expensive ALP fit is gated by `allow_alp`, which callers set from
//! [`crate::output::outfmt_needs_stats`] OR an explicit `--evalue` — so the
//! fit (and all per-HSP E-value math) is free when statistics are neither
//! printed nor used as a threshold.  `NcbiCompat` never runs the ALP fit.
//!
//! **Thread safety:** the ALP fitter uses process-global RNG state.  Call
//! [`resolve`] once, from a single thread, before spawning search threads
//! (its result is immutable afterwards and safe to share).

use crate::matrix::ScoreMatrix;

pub use rmstats::{MatrixCliOverrides, RmStats};

/// Statistics resolution mode (see module docs).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StatsMode {
    /// Full hierarchy for reporting; reap only on explicit `--evalue`.
    Default,
    /// Emulate NCBI 2.17.1: native sources only, C-table warts included,
    /// reap always on (threshold defaulted to 10 by the caller).
    NcbiCompat,
}

/// Which level of the hierarchy supplied the Karlin-Altschul parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KaSource {
    /// Baked rmblast_*_values table entry (NCBI Mode 1).
    BakedTable,
    /// -matrix_lambda/-matrix_k/-matrix_alpha CLI overrides (NCBI Mode 2).
    Cli,
    /// `# KARLIN` comment line in the matrix file.
    MatrixComment,
    /// Startup-time ALP Gumbel fit of the matrix.
    AlpFit,
    /// No source available: NCBI Mode-3 sentinel (evalue 1.0 / bits 0.0).
    Sentinel,
}

impl KaSource {
    pub fn describe(&self) -> &'static str {
        match self {
            KaSource::BakedTable    => "baked RMBlast table (NCBI Mode 1)",
            KaSource::Cli           => "-matrix_lambda/k/alpha CLI overrides (NCBI Mode 2)",
            KaSource::MatrixComment => "# KARLIN comment in matrix file",
            KaSource::AlpFit        => "startup ALP Gumbel fit",
            KaSource::Sentinel      => "none (sentinel: evalue=1.0, bitscore=0.0)",
        }
    }
}

/// How to build an [`RmStats`] for a query (query-independent part).
#[derive(Debug, Clone, Copy)]
enum StatsRecipe {
    /// Baked-table lookup by basename; `user_cli` passes through untouched
    /// because NCBI still consults it for the alpha/beta arm of the length
    /// adjustment even on a table hit.
    Table { user_cli: MatrixCliOverrides },
    /// Mode-2 semantics from (possibly synthesized) CLI overrides.
    CliLike(MatrixCliOverrides),
    /// Mode-3 sentinel.
    Sentinel,
}

/// Query-independent statistics context, resolved once per search setup.
#[derive(Debug, Clone)]
pub struct KaContext {
    matrix_basename: String,
    gap_open: i32,
    gap_extend: i32,
    pub source: KaSource,
    /// Recipe for the *reported* evalue/bitscore columns (full hierarchy).
    report: StatsRecipe,
    /// Recipe for the -evalue reap; `None` = NCBI would be in sentinel mode
    /// for this invocation, so nothing is culled.
    reap: Option<StatsRecipe>,
    /// The gapped Karlin block the reporting recipe resolves to (for logging).
    pub kbp_gap: rmstats::KarlinBlk,
}

/// Basename of the -matrix argument, matching the normalization used by
/// `search::ka_cutoff::lookup_ka` (RepeatMasker passes full matrix paths).
pub fn matrix_basename(matrix_name: &str) -> &str {
    std::path::Path::new(matrix_name)
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or(matrix_name)
}

/// NCBI's placeholder entry for 30p53g.matrix @22/5 in blast_stat.c,
/// expressed as Mode-2 overrides (identical semantics: H = lambda/alpha =
/// 0.1, alpha = 1, beta = 0 — exactly the placeholder row's values).
const NCBI_30P53G_PLACEHOLDER: MatrixCliOverrides = MatrixCliOverrides {
    matrix_lambda: 0.1,
    matrix_k: 0.1,
    matrix_alpha: 1.0,
    matrix_beta: 0.0,
};

/// Resolve the statistics context for this matrix + gap costs.
///
/// Call once at startup, single-threaded (the ALP fitter is not thread-safe).
/// `allow_alp` gates the expensive ALP-fit fallback: pass true only when
/// statistics will actually be consumed (evalue/bitscore in the outfmt, or an
/// explicit --evalue threshold).  Ignored in `NcbiCompat` mode, which never
/// fits.
pub fn resolve(
    matrix: &ScoreMatrix,
    matrix_name: &str,
    gap_open: i32,
    gap_extend: i32,
    user_cli: &MatrixCliOverrides,
    allow_alp: bool,
    mode: StatsMode,
) -> KaContext {
    #[cfg(not(feature = "alp-fit"))]
    let _ = allow_alp;

    let basename = matrix_basename(matrix_name).to_owned();
    let table_hit = rmstats::rmblast_tables::karlin_blk_gapped_load_from_tables(
        gap_open, gap_extend, &basename,
    )
    .is_ok();
    // The rmstats table adds comparison.matrix; NCBI's C table does not have
    // it, so it must not count as a table hit in NcbiCompat mode.
    let is_comparison = basename.eq_ignore_ascii_case("comparison.matrix");
    let is_30p53g = basename.eq_ignore_ascii_case("30p53g.matrix");
    let ncbi_table_hit = table_hit && !is_comparison;

    let mk = |source: KaSource, report: StatsRecipe, reap: Option<StatsRecipe>| {
        let kbp_gap = kbp_for_recipe(&report, &basename, gap_open, gap_extend);
        KaContext {
            matrix_basename: basename.clone(),
            gap_open,
            gap_extend,
            source,
            report,
            reap,
            kbp_gap,
        }
    };

    if mode == StatsMode::NcbiCompat {
        // Emulate 2.17.1 exactly: native sources only, C-table warts and all.
        // The reap recipe always equals the report recipe — that's what NCBI
        // does (one kbp_gap drives both the printed columns and the cull).
        if ncbi_table_hit {
            let recipe = if is_30p53g {
                // The C table's 30p53g row is a placeholder; NCBI both prints
                // and reaps with it.
                StatsRecipe::CliLike(NCBI_30P53G_PLACEHOLDER)
            } else {
                StatsRecipe::Table { user_cli: *user_cli }
            };
            return mk(KaSource::BakedTable, recipe, Some(recipe));
        }
        if user_cli.is_set() {
            let recipe = StatsRecipe::CliLike(*user_cli);
            return mk(KaSource::Cli, recipe, Some(recipe));
        }
        // Sentinel: evalue 1.0 <= any threshold, so NCBI's reap culls
        // nothing; reap None spares the caller the no-op work.
        return mk(KaSource::Sentinel, StatsRecipe::Sentinel, None);
    }

    // Default mode: full hierarchy for the printed columns; the reap recipe
    // mirrors the report recipe and only takes effect if the caller supplies
    // an explicit --evalue threshold.

    // Baked table (the improved one: real 30p53g fit, comparison.matrix
    // included) by basename + exact gap costs.
    if table_hit {
        let recipe = StatsRecipe::Table { user_cli: *user_cli };
        return mk(KaSource::BakedTable, recipe, Some(recipe));
    }

    // CLI overrides.
    if user_cli.is_set() {
        let recipe = StatsRecipe::CliLike(*user_cli);
        return mk(KaSource::Cli, recipe, Some(recipe));
    }

    // Mode 2.3: # KARLIN comment in the matrix file.  Usable when it gives
    // lambda, K, and one of alpha/H (alpha = lambda/H when only H is given).
    if let Some(kc) = &matrix.karlin {
        let alpha = if kc.alpha > 0.0 {
            kc.alpha
        } else if kc.h > 0.0 {
            kc.lambda / kc.h
        } else {
            0.0
        };
        if kc.lambda > 0.0 && kc.k > 0.0 && alpha > 0.0 {
            let cli_eff = MatrixCliOverrides {
                matrix_lambda: kc.lambda,
                matrix_k: kc.k,
                matrix_alpha: alpha,
                matrix_beta: kc.beta,
            };
            let recipe = StatsRecipe::CliLike(cli_eff);
            return mk(KaSource::MatrixComment, recipe, Some(recipe));
        }
    }

    // Mode 2.5: ALP fit at startup.  Needs the # FREQS background; fit
    // failure falls through to the sentinel.
    //
    // Deterministic (LAST-style) mode: the default wall-clock-budget mode
    // makes the realization count — and therefore the fitted values — depend
    // on machine load (observed: 0.10908 vs 0.10881 lambda for the same
    // matrix on cold vs warm runs).  Reported E-values must not change
    // between reruns, so trade bit-matching the (equally timing-dependent)
    // baked-table trajectories for run-to-run reproducibility.
    #[cfg(feature = "alp-fit")]
    {
        if allow_alp && matrix.freqs[..4].iter().sum::<f64>() > 0.0 {
            let opts = rmstats::alp::AlpFitOptions {
                deterministic: true,
                max_time: 60.0, // LAST's convention for this mode's parameter
                ..Default::default()
            };
            match rmstats::alp::fit_gumbel_blastna(
                &matrix.scores,
                &matrix.freqs,
                gap_open,
                gap_extend,
                &opts,
            ) {
                Ok(fit) => {
                    // NCBI table mapping: alpha = a_J + a_I, beta = 0 (the
                    // baked Mode-1 entries approximate beta = 0 too).
                    let cli_eff = MatrixCliOverrides {
                        matrix_lambda: fit.raw.lambda,
                        matrix_k: fit.raw.k,
                        matrix_alpha: fit.raw.a_j + fit.raw.a_i,
                        matrix_beta: 0.0,
                    };
                    let recipe = StatsRecipe::CliLike(cli_eff);
                    return mk(KaSource::AlpFit, recipe, Some(recipe));
                }
                Err(e) => {
                    eprintln!(
                        "# rmblastn: ALP Gumbel fit failed for matrix '{}' \
                         (gapopen={} gapextend={}): {}; E-values unavailable",
                        basename, gap_open, gap_extend, e
                    );
                }
            }
        }
    }

    // Mode 3: sentinel.
    mk(KaSource::Sentinel, StatsRecipe::Sentinel, None)
}

/// The query-independent gapped Karlin block a recipe resolves to.
fn kbp_for_recipe(
    recipe: &StatsRecipe,
    basename: &str,
    gap_open: i32,
    gap_extend: i32,
) -> rmstats::KarlinBlk {
    match recipe {
        StatsRecipe::Table { user_cli } => rmstats::rmblast_tables::kbp_gapped_calc_custom_matrix(
            basename, gap_open, gap_extend, user_cli,
        ),
        StatsRecipe::CliLike(cli) => rmstats::KarlinBlk {
            lambda: cli.matrix_lambda,
            k: cli.matrix_k,
            log_k: cli.matrix_k.ln(),
            h: cli.matrix_lambda / cli.matrix_alpha,
        },
        StatsRecipe::Sentinel => rmstats::KarlinBlk::sentinel(),
    }
}

impl KaContext {
    /// Build the per-query stats for a recipe: pure math, safe on any thread.
    fn stats_for_recipe(
        &self,
        recipe: &StatsRecipe,
        query_length: i32,
        db_length: i64,
        db_num_seqs: i32,
    ) -> RmStats {
        match recipe {
            StatsRecipe::Table { user_cli } => RmStats::new_custom_matrix(
                &self.matrix_basename,
                self.gap_open,
                self.gap_extend,
                user_cli,
                query_length,
                db_length,
                db_num_seqs,
                0,    // -searchsp override: not exposed by this port's CLI
                None, // kbp_std: unreachable on the read_in_matrix path (see rmstats docs)
            ),
            StatsRecipe::CliLike(cli) => {
                // Mode-2 semantics regardless of table contents (used both
                // for genuine CLI overrides and for synthesized sources).
                let kbp_gap = kbp_for_recipe(recipe, &self.matrix_basename, self.gap_open, self.gap_extend);
                let inputs = rmstats::SearchSpaceInputs {
                    query_length,
                    db_length,
                    db_num_seqs,
                    eff_searchsp_override: 0,
                };
                let eff = rmstats::calc_eff_lengths_custom_matrix(
                    &inputs, &kbp_gap, None, cli, self.gap_open, self.gap_extend, true,
                );
                RmStats { kbp_gap, eff, round_down: false }
            }
            StatsRecipe::Sentinel => {
                let kbp_gap = rmstats::KarlinBlk::sentinel();
                let inputs = rmstats::SearchSpaceInputs {
                    query_length,
                    db_length,
                    db_num_seqs,
                    eff_searchsp_override: 0,
                };
                let eff = rmstats::calc_eff_lengths_custom_matrix(
                    &inputs,
                    &kbp_gap,
                    None,
                    &MatrixCliOverrides::default(),
                    self.gap_open,
                    self.gap_extend,
                    true,
                );
                RmStats { kbp_gap, eff, round_down: false }
            }
        }
    }

    /// Per-query statistics for the *reported* evalue/bitscore columns.
    pub fn for_query(&self, query_length: i32, db_length: i64, db_num_seqs: i32) -> RmStats {
        self.stats_for_recipe(&self.report, query_length, db_length, db_num_seqs)
    }

    /// Per-query statistics for the -evalue reap, or `None` when no reap
    /// statistics exist (sentinel).  Whether the reap actually runs is the
    /// caller's decision: always in `NcbiCompat` mode (threshold defaulted to
    /// 10), only on an explicit `--evalue` in `Default` mode.
    pub fn reap_stats_for_query(
        &self,
        query_length: i32,
        db_length: i64,
        db_num_seqs: i32,
    ) -> Option<RmStats> {
        let recipe = self.reap.as_ref()?;
        let stats = self.stats_for_recipe(recipe, query_length, db_length, db_num_seqs);
        stats.has_stats().then_some(stats)
    }

    /// True when reap statistics exist for this context (query-independent).
    /// Combined by the caller with the threshold decision (see
    /// [`Self::reap_stats_for_query`]).
    pub fn reap_available(&self) -> bool {
        self.reap.is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::ScoreMatrix;
    use std::io::Cursor;

    const COMPARISON_LIKE: &str = r"# FREQS A 0.265 C 0.235 G 0.235 T 0.265
   A   C   G   T   N
A  9 -15  -6 -17  -1
C -15  10 -15  -6  -1
G  -6 -15  10 -15  -1
T -17  -6 -15   9  -1
N  -1  -1  -1  -1  -1
";

    #[test]
    fn baked_table_hit_wins() {
        let m = ScoreMatrix::from_reader("20p41g.matrix", Cursor::new(COMPARISON_LIKE)).unwrap();
        for mode in [StatsMode::Default, StatsMode::NcbiCompat] {
            let ctx = resolve(&m, "/some/path/20p41g.matrix", 25, 5, &MatrixCliOverrides::default(), true, mode);
            assert_eq!(ctx.source, KaSource::BakedTable);
            assert!(ctx.kbp_gap.is_valid());
            assert!(ctx.reap_available());
        }
        // Gap-cost mismatch -> table miss.  Default: falls through to ALP
        // (reap stats then exist but only engage on explicit --evalue);
        // NcbiCompat: sentinel, nothing to reap.
        let ctx2 = resolve(&m, "20p41g.matrix", 11, 1, &MatrixCliOverrides::default(), true, StatsMode::Default);
        assert_ne!(ctx2.source, KaSource::BakedTable);
        let ctx3 = resolve(&m, "20p41g.matrix", 11, 1, &MatrixCliOverrides::default(), true, StatsMode::NcbiCompat);
        assert_eq!(ctx3.source, KaSource::Sentinel);
        assert!(!ctx3.reap_available());
    }

    #[test]
    fn cli_overrides_used_on_table_miss() {
        let m = ScoreMatrix::from_reader("nosuch.matrix", Cursor::new(COMPARISON_LIKE)).unwrap();
        let cli = MatrixCliOverrides {
            matrix_lambda: 0.1,
            matrix_k: 0.02,
            matrix_alpha: 0.3,
            matrix_beta: -25.0,
        };
        for mode in [StatsMode::Default, StatsMode::NcbiCompat] {
            let ctx = resolve(&m, "nosuch.matrix", 20, 5, &cli, true, mode);
            assert_eq!(ctx.source, KaSource::Cli);
            assert!((ctx.kbp_gap.lambda - 0.1).abs() < 1e-15);
            assert!((ctx.kbp_gap.h - 0.1 / 0.3).abs() < 1e-15);
            assert!(ctx.reap_available());
        }
    }

    #[test]
    fn karlin_comment_used_when_table_and_cli_miss() {
        let text = format!("# KARLIN lambda 0.1276 k 0.0179 h 0.3562\n{}", COMPARISON_LIKE);
        let m = ScoreMatrix::from_reader("nosuch.matrix", Cursor::new(text)).unwrap();
        let ctx = resolve(&m, "nosuch.matrix", 20, 5, &MatrixCliOverrides::default(), true, StatsMode::Default);
        assert_eq!(ctx.source, KaSource::MatrixComment);
        assert!((ctx.kbp_gap.lambda - 0.1276).abs() < 1e-15);
        // alpha derived from H: alpha = lambda/H, then Mode-2 re-derives
        // H = lambda/alpha, recovering the declared H.
        assert!((ctx.kbp_gap.h - 0.3562).abs() < 1e-12);
        // Reap stats exist (engage only on explicit --evalue), and the
        // NcbiCompat mode ignores the comment entirely (sentinel).
        assert!(ctx.reap_available());
        let stats = ctx.for_query(1_000_000, 50_000, 25);
        assert!(stats.has_stats());
        assert!(stats.evalue(300) > 0.0);
        assert!(stats.bit_score(300) > 0.0);
        let compat = resolve(&m, "nosuch.matrix", 20, 5, &MatrixCliOverrides::default(), true, StatsMode::NcbiCompat);
        assert_eq!(compat.source, KaSource::Sentinel);
    }

    #[test]
    fn sentinel_when_no_source() {
        // No FREQS line -> ALP fit is also unavailable.
        const NO_FREQS: &str = r"   A   C   G   T
A  9 -15  -6 -17
C -15  10 -15  -6
G  -6 -15  10 -15
T -17  -6 -15   9
";
        let m = ScoreMatrix::from_reader("nosuch.matrix", Cursor::new(NO_FREQS)).unwrap();
        let ctx = resolve(&m, "nosuch.matrix", 20, 5, &MatrixCliOverrides::default(), true, StatsMode::Default);
        assert_eq!(ctx.source, KaSource::Sentinel);
        assert!(!ctx.reap_available());
        let stats = ctx.for_query(1000, 1000, 1);
        assert!(!stats.has_stats());
        assert_eq!(stats.evalue(300), 1.0);
        assert_eq!(stats.bit_score(300), 0.0);
    }

    #[test]
    fn comparison_matrix_modes() {
        // comparison.matrix @20/5 is in the rmstats table but NOT in NCBI's
        // C table.  Default: report (and, on explicit --evalue, reap with)
        // the baked fit.  NcbiCompat: sentinel, exactly like 2.17.1.
        let m = ScoreMatrix::from_reader("comparison.matrix", Cursor::new(COMPARISON_LIKE)).unwrap();
        let ctx = resolve(&m, "comparison.matrix", 20, 5, &MatrixCliOverrides::default(), true, StatsMode::Default);
        assert_eq!(ctx.source, KaSource::BakedTable);
        assert!(ctx.kbp_gap.is_valid());
        assert!((ctx.kbp_gap.lambda - 0.09836223).abs() < 1e-12);
        assert!(ctx.reap_available());
        let compat = resolve(&m, "comparison.matrix", 20, 5, &MatrixCliOverrides::default(), true, StatsMode::NcbiCompat);
        assert_eq!(compat.source, KaSource::Sentinel);
        assert!(!compat.reap_available());
        assert!(compat.reap_stats_for_query(1_000_000, 50_000, 25).is_none());
        let stats = compat.for_query(1_000_000, 50_000, 25);
        assert_eq!(stats.evalue(300), 1.0);
        assert_eq!(stats.bit_score(300), 0.0);
        // With full CLI overrides NCBI is in Mode 2 (its table misses
        // comparison.matrix) -> CLI wins and the reap stats exist.
        let cli = MatrixCliOverrides {
            matrix_lambda: 0.09836223,
            matrix_k: 0.0789918084,
            matrix_alpha: 0.881279665,
            matrix_beta: -26.8400935093,
        };
        let ctx2 = resolve(&m, "comparison.matrix", 20, 5, &cli, true, StatsMode::NcbiCompat);
        assert_eq!(ctx2.source, KaSource::Cli);
        assert!(ctx2.reap_available());
    }

    #[test]
    fn placeholder_for_30p53g_in_compat_mode() {
        // NCBI's C table carries a placeholder (lambda=K=H=0.1) for 30p53g
        // @22/5 and both prints and reaps with it; the default mode uses the
        // real ALP fit for both instead.
        let m = ScoreMatrix::from_reader("30p53g.matrix", Cursor::new(COMPARISON_LIKE)).unwrap();
        let compat = resolve(&m, "30p53g.matrix", 22, 5, &MatrixCliOverrides::default(), true, StatsMode::NcbiCompat);
        assert_eq!(compat.source, KaSource::BakedTable);
        assert!(compat.reap_available());
        assert!((compat.kbp_gap.lambda - 0.1).abs() < 1e-15);
        let reap = compat.reap_stats_for_query(1_000_000, 50_000, 25).unwrap();
        assert!((reap.kbp_gap.lambda - 0.1).abs() < 1e-15);
        assert!((reap.kbp_gap.k - 0.1).abs() < 1e-15);
        assert!((reap.kbp_gap.h - 0.1).abs() < 1e-15);

        let ctx = resolve(&m, "30p53g.matrix", 22, 5, &MatrixCliOverrides::default(), true, StatsMode::Default);
        assert!((ctx.kbp_gap.lambda - 0.1204393447).abs() < 1e-12);
        let report = ctx.for_query(1_000_000, 50_000, 25);
        assert!((report.kbp_gap.lambda - 0.1204393447).abs() < 1e-12);
        let reap = ctx.reap_stats_for_query(1_000_000, 50_000, 25).unwrap();
        assert!((reap.kbp_gap.lambda - 0.1204393447).abs() < 1e-12);
    }

    #[cfg(feature = "alp-fit")]
    #[test]
    fn alp_fit_close_to_baked_and_reproducible() {
        use crate::matrix::ScoreMatrix;
        // 20p41g with its FREQS background at its canonical 25/5 gap costs,
        // but under a name the table cannot match: the ALP fit must fire.
        // The deterministic-mode sampling schedule differs from the one the
        // baked table was fitted with, so values are only required to agree
        // to a few percent — but they MUST be identical across runs (that's
        // the property the deterministic mode buys; the wall-clock mode was
        // observed to drift with machine load).
        let baked =
            rmstats::rmblast_tables::karlin_blk_gapped_load_from_tables(25, 5, "20p41g.matrix")
                .unwrap();

        let matrix_path = std::env::var("BLASTMAT")
            .map(|d| format!("{}/20p41g.matrix", d))
            .ok()
            .filter(|p| std::path::Path::new(p).exists());
        let Some(matrix_path) = matrix_path else {
            eprintln!("skipping: BLASTMAT not set or 20p41g.matrix missing");
            return;
        };
        let m = ScoreMatrix::from_file(&matrix_path).unwrap();
        let cli = MatrixCliOverrides::default();
        let ctx = resolve(&m, "renamed-so-table-misses.matrix", 25, 5, &cli, true, StatsMode::Default);
        assert_eq!(ctx.source, KaSource::AlpFit);
        assert!(ctx.reap_available());
        let rel = |a: f64, b: f64| ((a - b) / b).abs();
        assert!(
            rel(ctx.kbp_gap.lambda, baked.lambda) < 0.02,
            "lambda {} vs baked {}",
            ctx.kbp_gap.lambda,
            baked.lambda
        );
        assert!(
            rel(ctx.kbp_gap.k, baked.k) < 0.15,
            "K {} vs baked {}",
            ctx.kbp_gap.k,
            baked.k
        );
        // Reproducibility: a second resolve must give bit-identical params.
        let ctx2 = resolve(&m, "renamed-so-table-misses.matrix", 25, 5, &cli, true, StatsMode::Default);
        assert_eq!(ctx.kbp_gap, ctx2.kbp_gap, "deterministic fit must be reproducible");
    }
}
