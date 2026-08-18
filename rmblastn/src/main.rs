//! rmblastn — Rust port of NCBI rmblastn (BLAST 2.17.0)
//!
//! Nucleotide BLAST with custom substitution matrices, UCSC 2bit database
//! support, and RepeatMasker-specific output fields.

use std::collections::HashSet;
use std::io::{BufRead, BufWriter, Write};

use anyhow::{Context, Result};
use clap::Parser;
use rayon::prelude::*;

use rmblast_lib::matrix::ScoreMatrix;
use rmblast_lib::hits::Strand;
use rmblast_lib::ka_stats::{self, KaContext, MatrixCliOverrides, RmStats};
use rmblast_lib::options::{MtMode, SearchParams, SeedMode};
use rmblast_lib::output::{
    outfmt_needs_stats, parse_outfmt, write_tabular, write_pairwise_program_header,
    write_pairwise_results, write_pairwise_footer, AlignResult, OutField,
};
use rmblast_lib::search::{apply_mask_level, build_query_lookup, build_query_lookup_premask, mask_query_for_alignment, search_with_query_lookup, search_phase2a, run_phase2b, PrelimHsp};
use rmblast_lib::search::engine::{cull_prelims, resurrect_prelims, PrelimCullParams};
// use rmblast_lib::search::engine::{COUNT_SEEDS, COUNT_UNGAPPED_HITS, COUNT_PRELIM_GAPPED, COUNT_FINAL_GAPPED, COUNT_FINAL_HITS};
// use rmblast_lib::search::gapped::TOTAL_DP_CELLS;
use rmblast_lib::seq::{FastaReader, SubjectDb};

// ──────────────────────────────────────────────────────────────────────────────
// Argument-normalization shim for NCBI single-dash compatibility
// ──────────────────────────────────────────────────────────────────────────────
//
// NCBI's C++ Toolkit argument parser accepts long options with a single leading
// dash (e.g. `-word_size 8`, `-gapopen 4`).  clap — the Rust argument parser
// used here — follows POSIX/GNU convention and requires double-dash for long
// options (`--word_size 8`, `--gapopen 4`).
//
// `normalize_args()` bridges the gap by rewriting any token of the form
// `-<letter><rest>` to `--<letter><rest>` before clap sees argv.  Short
// single-character options (`-h`, `-v`), the stdin marker (`-`), and the
// end-of-options sentinel (`--`) are left unchanged.
//
// IMPORTANT LIMITATION: Because normalization happens before clap parses, the
// POSIX convention of bundling single-character flags into one token is NOT
// supported.  `-a -d -r` must be written as three separate arguments; `-adr`
// will be rewritten to `--adr` and rejected as an unknown option.
//
// DEPRECATION NOTICE: Single-dash long options are a compatibility shim for
// users migrating from NCBI rmblastn.  New scripts should use double-dash
// (`--`) prefixes.  The single-dash form may be removed in a future release.
//
/// Rewrite NCBI-style single-dash long options to double-dash before clap parsing.
///
/// Rule: `-<letter><one-or-more-chars>` → `--<letter><one-or-more-chars>`
/// All other tokens (bare `-`, `--`, `--foo`, `-x`, values) pass through unchanged.
fn normalize_args() -> Vec<std::ffi::OsString> {
    std::env::args_os()
        .map(|arg| {
            let s = match arg.to_str() {
                Some(s) => s,
                None => return arg,
            };
            if s.starts_with("--") || s == "-" {
                return arg;
            }
            if let Some(rest) = s.strip_prefix('-') {
                let mut chars = rest.chars();
                if let Some(first) = chars.next() {
                    if first.is_ascii_alphabetic() && chars.next().is_some() {
                        return format!("--{}", rest).into();
                    }
                }
            }
            arg
        })
        .collect()
}
// ──────────────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug)]
#[command(
    name = "rmblastn",
    version,
    about = "RepeatMasker BLAST — Rust port of NCBI rmblastn",
    after_help = "\
ARGUMENT SYNTAX NOTES:
  Both single-dash and double-dash long options are accepted:
    -word_size 8        (NCBI rmblastn style, compatibility shim)
    --word_size 8       (preferred)

  Single-character flags CANNOT be bundled into one token.
  Write '-a -b -c' as three separate arguments, not '-abc'.
  ('-abc' will be rewritten to '--abc' and rejected as unknown.)

  The single-dash compatibility form may be removed in a future release.
"
)]
struct Args {
    /// UCSC 2bit database file
    #[arg(long)]
    db: String,

    /// FASTA query file (use '-' for stdin)
    #[arg(long)]
    query: String,

    /// Scoring matrix file (RMBlast FREQS format)
    #[arg(long, default_value = "")]
    matrix: String,

    /// Gap open penalty
    #[arg(long, default_value_t = 4)]
    gapopen: i32,

    /// Gap extend penalty
    #[arg(long, default_value_t = 4)]
    gapextend: i32,

    /// Word size for seeding
    #[arg(long, default_value_t = 8)]
    word_size: usize,

    /// X-dropoff for ungapped extension
    #[arg(long, default_value_t = 20)]
    xdrop_ungap: i32,

    /// X-dropoff for gapped extension
    #[arg(long, default_value_t = 30)]
    xdrop_gap: i32,

    /// X-dropoff for final gapped extension
    #[arg(long, default_value_t = 100)]
    xdrop_gap_final: i32,

    /// Cutoff for the preliminary gapped stage. NOT a floor on the reported
    /// score: traceback re-aligns under --xdrop-gap-final and may score lower,
    /// and that lower score is reported (matches NCBI rmblastn). Post-filter if
    /// you need a hard floor.
    #[arg(long, default_value_t = 0)]
    min_raw_gapped_score: i32,

    /// Use complexity adjusted scoring (bare flag, no value required)
    #[arg(long, action = clap::ArgAction::SetTrue)]
    complexity_adjust: bool,

    /// Filter query with DUST: 'yes' to enable, 'no' to disable.
    /// NCBI also accepts 'level window linker' (e.g. '20 64 1'); custom
    /// parameters are accepted but ignored — only yes/no is acted on.
    #[arg(long, default_value = "yes")]
    dust: String,

    /// Number of threads
    #[arg(long, default_value_t = 1)]
    num_threads: usize,

    /// Threading mode: 0 = split by DB, 1 = split by queries
    #[arg(long, default_value_t = 0)]
    mt_mode: u8,

    /// Output format: 0 = pairwise (default, matches NCBI rmblastn), 6 [fields...] = tabular.
    #[arg(long, default_value = "0")]
    outfmt: String,

    /// Masklevel: drop an HSP if >N% of its query span is covered by a better HSP (default 80; 101 = disabled)
    #[arg(long, default_value_t = 80)]
    mask_level: u32,

    /// File with subject sequence IDs to include (one per line)
    #[arg(long)]
    gilist: Option<String>,

    /// Lookup-table seeding strategy: combined (NCBI-faithful) or separate (efficient).
    ///
    /// combined: indexes k-mers from both forward and RC query in one table, matching
    ///   NCBI's BlastMBLookupTableNew exactly (including ~2× wasted context-1 hits).
    ///   DUST is applied to the query before LUT construction.
    ///
    /// separate: indexes only forward-query k-mers; minus-strand seeding uses the
    ///   same LUT against the revcomp subject.  No wasted context-1 hits.
    ///   DUST is applied to the query before LUT construction.
    #[arg(long, default_value = "combined")]
    seed_mode: String,

    /// Write all merged Phase 2a prelim HSPs to this file before Phase 2b runs.
    /// Format: sseqid TAB strand TAB qstart TAB qend TAB sstart TAB send TAB score
    #[arg(long)]
    dump_prelims: Option<String>,

    /// FAST MODE (not NCBI-faithful): skip Phase 2b traceback for preliminary
    /// HSPs dominated (mask_level-style) by higher-scoring prelims, then
    /// re-test against the survivors' final alignments and resurrect any no
    /// longer dominated.  Output differs from faithful mode in a small
    /// fraction of heavily-overlapped, mostly low-scoring hits (chr22
    /// benchmark: balanced ~0.2% of annotated bp).  Modes: off (default),
    /// conservative, balanced, aggressive.  Active only for multi-chunk
    /// queries (>1 Mbp) with mask_level < 100.
    #[arg(long, default_value = "off")]
    prelim_cull: String,

    /// Karlin-Altschul lambda override for E-value statistics of a matrix
    /// absent from the baked table (NCBI -matrix_lambda; needs _k and _alpha too)
    #[arg(long, alias = "matrix_lambda", default_value_t = 0.0)]
    matrix_lambda: f64,

    /// Karlin-Altschul K override (NCBI -matrix_k)
    #[arg(long, alias = "matrix_k", default_value_t = 0.0)]
    matrix_k: f64,

    /// Karlin-Altschul alpha override (NCBI -matrix_alpha; H = lambda/alpha)
    #[arg(long, alias = "matrix_alpha", default_value_t = 0.0)]
    matrix_alpha: f64,

    /// Karlin-Altschul beta override for the length adjustment (NCBI -matrix_beta)
    #[arg(long, alias = "matrix_beta", default_value_t = 0.0, allow_hyphen_values = true)]
    matrix_beta: f64,

    /// Expect-value threshold: drop HSPs with E-value above this.  By
    /// default NO E-value culling happens unless this flag is given (the
    /// statistics then come from the full source hierarchy: baked table,
    /// -matrix_* overrides, # KARLIN comment, or ALP fit).  Under
    /// --ncbi-compat this defaults to 10, matching NCBI.
    #[arg(long)]
    evalue: Option<f64>,

    /// Emulate NCBI rmblastn 2.17.1 E-value behavior exactly: statistics
    /// only from the hardcoded matrix table (warts included: 30p53g's
    /// placeholder values) or -matrix_* overrides, sentinel 1.0/0.0
    /// rendering otherwise, and HSPs with E-value above --evalue (default
    /// 10) silently dropped before masklevel — 2.17.1 has done this for
    /// table matrices since the table was introduced.  Use for parity
    /// benchmarks and legacy pipelines.
    #[arg(long, alias = "ncbi_compat", action = clap::ArgAction::SetTrue)]
    ncbi_compat: bool,
}

/// Resolve a `--db` argument to an existing file path.
///
/// Accepted forms (tried in order):
///   1. Exact path — used as-is if it exists.
///   2. With each of these extensions appended in order: `.2bit`, `.fa`,
///      `.fasta`, `.fas`, `.fna`.
///
/// The caller should pass the result to `SubjectDb::open`, which detects
/// the format automatically from the extension.
fn resolve_db(arg: &str) -> Result<String> {
    let p = std::path::Path::new(arg);
    if p.exists() {
        return Ok(arg.to_owned());
    }
    for ext in &[".2bit", ".fa", ".fasta", ".fas", ".fna"] {
        let candidate = format!("{}{}", arg, ext);
        if std::path::Path::new(&candidate).exists() {
            return Ok(candidate);
        }
    }
    anyhow::bail!(
        "database not found: '{}' does not exist; also tried suffixes \
         .2bit, .fa, .fasta, .fas, .fna",
        arg
    )
}

/// Output format selected by --outfmt.
enum OutputFormat {
    /// NCBI-style pairwise alignment output (outfmt 0, the default).
    Pairwise,
    /// Tab-separated tabular output (outfmt 6) with the given field list.
    Tabular(Vec<OutField>),
}

fn parse_output_format(s: &str) -> Result<OutputFormat> {
    let trimmed = s.trim();
    if trimmed == "0" || trimmed.is_empty() {
        return Ok(OutputFormat::Pairwise);
    }
    let first_token = trimmed.split_whitespace().next().unwrap_or("");
    if first_token == "0" {
        return Ok(OutputFormat::Pairwise);
    }
    // Attempt tabular parse (format 6 or bare field list).
    let fields = parse_outfmt(trimmed);
    if fields.is_empty() {
        anyhow::bail!(
            "--outfmt '{}' is not recognized; use '0' for pairwise or '6 <fields>' for tabular",
            s
        );
    }
    Ok(OutputFormat::Tabular(fields))
}

fn main() -> Result<()> {
    let args = Args::parse_from(normalize_args());

    if args.matrix.is_empty() {
        anyhow::bail!("--matrix is required");
    }

    let matrix_path = if std::path::Path::new(&args.matrix).exists() {
        args.matrix.clone()
    } else if let Ok(blastmat) = std::env::var("BLASTMAT") {
        let candidate = format!("{}/{}", blastmat, args.matrix);
        if std::path::Path::new(&candidate).exists() {
            candidate
        } else {
            anyhow::bail!("cannot find matrix '{}' (also tried BLASTMAT: '{}')", args.matrix, candidate);
        }
    } else {
        args.matrix.clone()
    };
    let matrix = ScoreMatrix::from_file(&matrix_path)
        .with_context(|| format!("loading matrix '{}'", matrix_path))?;

    let seed_mode = match args.seed_mode.as_str() {
        "separate" | "separate-strands" | "SeparateStrands" => SeedMode::SeparateStrands,
        _ => SeedMode::Combined,
    };

    let dust = !matches!(args.dust.to_lowercase().as_str(), "no" | "false" | "0");

    let params = SearchParams {
        gap_open: args.gapopen,
        gap_extend: args.gapextend,
        matrix_name: args.matrix.clone(),
        word_size: args.word_size,
        xdrop_ungap: args.xdrop_ungap,
        xdrop_gap: args.xdrop_gap,
        xdrop_gap_final: args.xdrop_gap_final,
        min_raw_gapped_score: args.min_raw_gapped_score,
        // Resolved per-search in search_db_parallel (KA cutoff when ALP params exist).
        ungapped_cutoff: None,
        complexity_adjust: args.complexity_adjust,
        dust,
        mask_level: args.mask_level,
        num_threads: args.num_threads,
        mt_mode: if args.mt_mode == 1 { MtMode::SplitByQueries } else { MtMode::SplitByDb },
        seed_mode,
    };

    let output_format = parse_output_format(&args.outfmt)?;

    rayon::ThreadPoolBuilder::new()
        .num_threads(params.num_threads)
        .build_global()
        .ok();

    let db_path = resolve_db(&args.db)?;
    let db = SubjectDb::open(&db_path)
        .with_context(|| format!("opening database '{}'", db_path))?;

    let gilist = load_gilist(args.gilist.as_deref())?;
    let subject_names: Vec<String> = db.sequences().iter()
        .map(|s| s.name.clone())
        .filter(|n| gilist.as_ref().map_or(true, |gl| gl.contains(n.as_str())))
        .collect();

    if subject_names.is_empty() {
        return Ok(());
    }

    // Compute database statistics for pairwise header/footer.
    let n_db_seqs = subject_names.len();
    let total_db_letters: u64 = db.sequences().iter()
        .filter(|s| subject_names.contains(&s.name))
        .map(|s| s.dna_size as u64)
        .sum();

    // Average DB sequence length (NCBI BlastSeqSrcGetAvgSeqLen: integer-truncated
    // total residues / number of sequences).  Feeds the Karlin-Altschul ungapped
    // cutoff for matrices with ALP params (see search::ka_cutoff).
    let avg_subj_length: u64 = if n_db_seqs > 0 {
        total_db_letters / n_db_seqs as u64
    } else {
        0
    };

    // Karlin-Altschul statistics context.
    //
    // Default mode: nothing is culled unless --evalue is given; the printed
    // evalue/bitscore columns (when requested) come from the full source
    // hierarchy.  Statistics are only resolved at all when they will be
    // consumed (columns requested or --evalue given) — the lazy-stats rule
    // that keeps the ALP fit free otherwise.
    //
    // --ncbi-compat: emulate 2.17.1 — native sources only, and the reap is
    // always on with NCBI's default threshold of 10.
    //
    // Resolution must happen here — single-threaded, before any searches —
    // because the ALP fitter uses process-global RNG state.
    let stats_requested = match &output_format {
        OutputFormat::Tabular(fields) => outfmt_needs_stats(fields),
        OutputFormat::Pairwise => false,
    };
    let stats_mode = if args.ncbi_compat { ka_stats::StatsMode::NcbiCompat } else { ka_stats::StatsMode::Default };
    let reap_threshold: Option<f64> = if args.ncbi_compat {
        Some(args.evalue.unwrap_or(10.0))
    } else {
        args.evalue
    };
    let user_cli = MatrixCliOverrides {
        matrix_lambda: args.matrix_lambda,
        matrix_k: args.matrix_k,
        matrix_alpha: args.matrix_alpha,
        matrix_beta: args.matrix_beta,
    };
    let ka_ctx: KaContext = ka_stats::resolve(
        &matrix, &args.matrix, args.gapopen, args.gapextend, &user_cli,
        stats_requested || reap_threshold.is_some(),
        stats_mode,
    );
    let reap_on = reap_threshold.is_some() && ka_ctx.reap_available();
    if stats_requested || reap_threshold.is_some() {
        eprintln!(
            "# rmblastn: E-value statistics{}: {}",
            if args.ncbi_compat { " (--ncbi-compat 2.17.1 emulation)" } else { "" },
            ka_ctx.source.describe()
        );
        if ka_ctx.kbp_gap.is_valid() {
            eprintln!(
                "#   lambda={:.10} K={:.10} H={:.10}",
                ka_ctx.kbp_gap.lambda, ka_ctx.kbp_gap.k, ka_ctx.kbp_gap.h
            );
        }
        match (reap_on, reap_threshold) {
            (true, Some(t)) => eprintln!("#   -evalue cutoff ACTIVE: threshold {}", t),
            (false, Some(t)) => eprintln!(
                "#   -evalue cutoff requested (threshold {}) but no statistics available — nothing culled",
                t
            ),
            (_, None) => eprintln!("#   -evalue cutoff: off (no --evalue given)"),
        }
    }
    // Reporting stats are built per query below; this closure builds the
    // reap context (None = nothing culled).
    let reap_for_query = |query_len: i32| -> Option<(RmStats, f64)> {
        let t = reap_threshold?;
        ka_ctx
            .reap_stats_for_query(query_len, total_db_letters as i64, n_db_seqs as i32)
            .map(|s| (s, t))
    };

    // Get the 2bit file modification time for the "Posted date" footer line.
    let db_mtime_unix: Option<u64> = std::fs::metadata(&db_path)
        .ok()
        .and_then(|m| m.modified().ok())
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map(|d| d.as_secs());

    // Fast-mode prelim cull presets (measured on chr22 x longlib; see
    // PrelimCullParams docs).  Requires an active mask_level.
    let prelim_cull: Option<PrelimCullParams> = match args.prelim_cull.as_str() {
        "off" => None,
        "conservative" => Some(PrelimCullParams {
            cull_margin_pct: 110, cull_coverage: 95,
            resurrect_margin_pct: 90, resurrect_slack: 100,
        }),
        "balanced" => Some(PrelimCullParams {
            cull_margin_pct: 110, cull_coverage: 95,
            resurrect_margin_pct: 100, resurrect_slack: 50,
        }),
        "aggressive" => Some(PrelimCullParams {
            cull_margin_pct: 100, cull_coverage: 95,
            resurrect_margin_pct: 100, resurrect_slack: 0,
        }),
        other => anyhow::bail!(
            "--prelim-cull: unknown mode '{}' (off|conservative|balanced|aggressive)", other
        ),
    };
    let prelim_cull = if prelim_cull.is_some() && params.mask_level >= 100 {
        eprintln!("Note: --prelim-cull requires mask_level < 100; running faithful pipeline.");
        None
    } else {
        prelim_cull
    };

    let query_reader: Box<dyn std::io::Read> = if args.query == "-" {
        Box::new(std::io::stdin())
    } else {
        Box::new(
            std::fs::File::open(&args.query)
                .with_context(|| format!("opening query '{}'", args.query))?,
        )
    };
    let query_reader = FastaReader::new(std::io::BufReader::new(query_reader));

    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    // Print the pairwise program/db header once (before any query output).
    if matches!(output_format, OutputFormat::Pairwise) {
        write_pairwise_program_header(&mut out, &args.db, n_db_seqs, total_db_letters)
            .context("writing pairwise header")?;
    }

    match params.mt_mode {
        MtMode::SplitByDb => {
            for qrec in query_reader {
                let qrec = qrec.context("reading query FASTA")?;
                let full_q_len = qrec.len() as u32;
                let mut results = search_db_parallel(
                    &qrec.seq, &qrec.id, &subject_names, &db, &params, &matrix,
                    avg_subj_length, args.dump_prelims.as_deref(), prelim_cull,
                    reap_for_query(full_q_len as i32),
                );
                results.sort_by(|a, b| {
                    let qa = a.hsp.q_len;
                    let qb = b.hsp.q_len;
                    let q_off_a = match a.hsp.strand { Strand::Plus => a.hsp.q_start, Strand::Minus => qa - a.hsp.q_end };
                    let q_off_b = match b.hsp.strand { Strand::Plus => b.hsp.q_start, Strand::Minus => qb - b.hsp.q_end };
                    let q_end_a = match a.hsp.strand { Strand::Plus => a.hsp.q_end, Strand::Minus => qa - a.hsp.q_start };
                    let q_end_b = match b.hsp.strand { Strand::Plus => b.hsp.q_end, Strand::Minus => qb - b.hsp.q_start };
                    b.hsp.score.cmp(&a.hsp.score)
                        .then_with(|| a.hsp.s_start.cmp(&b.hsp.s_start))
                        .then_with(|| b.hsp.s_end.cmp(&a.hsp.s_end))
                        .then_with(|| q_off_a.cmp(&q_off_b))
                        .then_with(|| q_end_b.cmp(&q_end_a))
                });
                match &output_format {
                    OutputFormat::Pairwise => {
                        write_pairwise_results(&mut out, &results, &qrec.defline, full_q_len)
                            .context("writing pairwise results")?;
                    }
                    OutputFormat::Tabular(fields) => {
                        let qstats = stats_requested.then(|| {
                            ka_ctx.for_query(full_q_len as i32, total_db_letters as i64, n_db_seqs as i32)
                        });
                        write_all(&mut out, &results, fields, qstats.as_ref())?;
                    }
                }
            }
        }
        MtMode::SplitByQueries => {
            let all_queries: Vec<_> = query_reader
                .collect::<Result<Vec<_>, _>>()
                .context("reading query FASTA")?;
            let all_results: Vec<Vec<AlignResult>> = all_queries
                .par_iter()
                .map(|q| {
                    let mut qr = search_db_parallel(
                        &q.seq, &q.id, &subject_names, &db, &params, &matrix,
                        avg_subj_length, None, prelim_cull,
                        reap_for_query(q.len() as i32),
                    );
                    qr.sort_by(|a, b| b.hsp.score.cmp(&a.hsp.score));
                    qr
                })
                .collect();
            for (q, results) in all_queries.iter().zip(all_results.into_iter()) {
                match &output_format {
                    OutputFormat::Pairwise => {
                        let qlen = q.len() as u32;
                        write_pairwise_results(&mut out, &results, &q.defline, qlen)
                            .context("writing pairwise results")?;
                    }
                    OutputFormat::Tabular(fields) => {
                        let qstats = stats_requested.then(|| {
                            ka_ctx.for_query(q.len() as i32, total_db_letters as i64, n_db_seqs as i32)
                        });
                        write_all(&mut out, &results, fields, qstats.as_ref())?;
                    }
                }
            }
        }
    }

    // Print the pairwise footer once (after all queries).
    if matches!(output_format, OutputFormat::Pairwise) {
        write_pairwise_footer(
            &mut out,
            &args.db,
            n_db_seqs,
            total_db_letters,
            db_mtime_unix,
            args.gapopen,
            args.gapextend,
        ).context("writing pairwise footer")?;
    }

    out.flush()?;

    if std::env::var_os("RMBLAST_DEBUG_COUNTERS").is_some() {
        eprintln!(
            "REVERSE_FBI_CLAMPED={} IMPROVE_SEED_NEGATIVE_OFFSET={}",
            rmblast_lib::search::gapped::REVERSE_FBI_CLAMPED
                .load(std::sync::atomic::Ordering::Relaxed),
            rmblast_lib::search::engine::IMPROVE_SEED_NEGATIVE_OFFSET
                .load(std::sync::atomic::Ordering::Relaxed),
        );
    }

    // use std::sync::atomic::Ordering;
    // let dp_cells = TOTAL_DP_CELLS.load(Ordering::Relaxed);
    // eprintln!("COUNTS seeds={} ungapped_hits={} prelim_gapped={} final_gapped={} final_hits={} dp_cells={}",
    //     COUNT_SEEDS.load(Ordering::Relaxed),
    //     COUNT_UNGAPPED_HITS.load(Ordering::Relaxed),
    //     COUNT_PRELIM_GAPPED.load(Ordering::Relaxed),
    //     COUNT_FINAL_GAPPED.load(Ordering::Relaxed),
    //     COUNT_FINAL_HITS.load(Ordering::Relaxed),
    //     dp_cells,
    // );

    if std::env::var_os("RMBLAST_DP_STATS").is_some() {
        use std::sync::atomic::Ordering;
        use rmblast_lib::search::gapped::{
            SCORE_ONLY_CELLS, SCORE_ONLY_ROWS, SCORE_ONLY_SIMD_ROWS, SCORE_ONLY_WIDE_CELLS,
            TOTAL_DP_CELLS,
        };
        use rmblast_lib::search::engine::{
            COUNT_FINAL_GAPPED, COUNT_PRELIM_GAPPED, COUNT_SEEDS, COUNT_UNGAPPED_HITS,
        };
        let so_cells = SCORE_ONLY_CELLS.load(Ordering::Relaxed);
        let total_cells = TOTAL_DP_CELLS.load(Ordering::Relaxed);
        eprintln!(
            "DP_STATS score_only: rows={} simd_rows={} cells={} wide_cells={}",
            SCORE_ONLY_ROWS.load(Ordering::Relaxed),
            SCORE_ONLY_SIMD_ROWS.load(Ordering::Relaxed),
            so_cells,
            SCORE_ONLY_WIDE_CELLS.load(Ordering::Relaxed),
        );
        eprintln!(
            "DP_STATS stages: seeds={} ungapped={} prelim_gapped={} final_gapped={} tb_cells={}",
            COUNT_SEEDS.load(Ordering::Relaxed),
            COUNT_UNGAPPED_HITS.load(Ordering::Relaxed),
            COUNT_PRELIM_GAPPED.load(Ordering::Relaxed),
            COUNT_FINAL_GAPPED.load(Ordering::Relaxed),
            total_cells - so_cells,
        );
    }

    Ok(())
}

/// Cross-chunk HSP merge, mirroring NCBI BlastHSPStreamMerge → Blast_HSPListsMerge
/// → s_BlastMergeTwoHSPs.
///
/// When chunk k+1 (starting at `split_point` in global coords) has been processed,
/// boundary candidates from `accumulated` (those extending past the split) are merged
/// with boundary candidates from `new_prelims` (those starting within the 100-base
/// overlap).  Merged results update `accumulated` in-place; unmerged new entries are
/// appended.
fn merge_chunk_prelims(
    accumulated: &mut Vec<PrelimHsp>,
    new_prelims: Vec<PrelimHsp>,
    split_point: u32,
) {
    const OVERLAP: u32 = 100;
    const OVERLAP_DIAG_CLOSE: i64 = 10;

    let mut to_add: Vec<PrelimHsp> = Vec::new();

    'new: for new in new_prelims {
        // New-chunk HSPs are boundary candidates iff they start within the overlap region.
        // (q_start >= split_point is guaranteed by chunk construction.)
        if new.q_start >= split_point + OVERLAP {
            to_add.push(new);
            continue;
        }

        for acc in accumulated.iter_mut() {
            if acc.strand != new.strand { continue; }
            // Accumulated boundary candidate: extends into the overlap.
            if acc.q_end <= split_point { continue; }

            // Diagonal proximity (same formula for both strands in Rust FWD-q + strand-s space):
            // |(acc.q_end - acc.s_end) - (new.q_start - new.s_start)| < OVERLAP_DIAG_CLOSE
            let diag_acc_end   = acc.q_end as i64 - acc.s_end as i64;
            let diag_new_start = new.q_start as i64 - new.s_start as i64;
            if (diag_acc_end - diag_new_start).abs() >= OVERLAP_DIAG_CLOSE { continue; }

            // CONTAINED_IN_HSP: at least one endpoint of `new` falls inside `acc`'s bbox.
            let start_in = new.q_start >= acc.q_start && new.q_start <= acc.q_end
                        && new.s_start >= acc.s_start && new.s_start <= acc.s_end;
            let end_in   = new.q_end   >= acc.q_start && new.q_end   <= acc.q_end
                        && new.s_end   >= acc.s_start && new.s_end   <= acc.s_end;
            if !start_in && !end_in { continue; }

            // Merge: expand bounding box; keep higher-scoring HSP's seed.
            let (q_seed, s_seed) = if new.score > acc.score {
                (new.q_seed, new.s_seed)
            } else {
                (acc.q_seed, acc.s_seed)
            };
            acc.q_start = acc.q_start.min(new.q_start);
            acc.q_end   = acc.q_end.max(new.q_end);
            acc.s_start = acc.s_start.min(new.s_start);
            acc.s_end   = acc.s_end.max(new.s_end);
            acc.score   = acc.score.max(new.score);
            acc.q_seed  = q_seed;
            acc.s_seed  = s_seed;
            continue 'new;
        }

        // No accumulated candidate matched: add as a standalone entry.
        to_add.push(new);
    }

    accumulated.extend(to_add);
}

/// Search one query against all subjects, chunking the query to match NCBI's
/// BlastSplitQuery behavior (chunk_size=1_000_000, overlap=100, step=999_900).
///
/// 3-phase pipeline (multi-chunk): Phase 1+2a per chunk (parallel over subjects) →
/// cross-chunk merge per subject (sequential over chunks) → Phase 2b with full query
/// (parallel over subjects).  Single-chunk queries use the original run_gapped_phase
/// path unchanged.
/// Emit a one-time note to stderr when a Karlin-Altschul ungapped cutoff is in
/// effect (i.e. the matrix+gap combo matched the ALP table), making it clear in
/// the run log that the ungapped pre-filter threshold is the KA-derived value
/// rather than the default `min_raw_gapped_score / 2`.
fn log_ungapped_cutoff(
    params: &SearchParams,
    cutoff: i32,
    ka: Option<&rmblast_lib::search::KaCutoffInfo>,
) {
    use std::sync::Once;
    static ONCE: Once = Once::new();
    if let Some(info) = ka {
        ONCE.call_once(|| {
            let base = std::path::Path::new(&params.matrix_name)
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or(&params.matrix_name);
            eprintln!(
                "# rmblastn: Karlin-Altschul ungapped cutoff IN EFFECT (matrix={} gapopen={} gapextend={})",
                base, params.gap_open, params.gap_extend
            );
            eprintln!(
                "#   lambda={:.10} K={:.10} E={} avg_db_seq_len={} searchsp={} -> KA score={}",
                info.lambda,
                info.k,
                rmblast_lib::search::ka_cutoff::CUTOFF_E_BLASTN,
                info.avg_subj_length,
                info.searchsp,
                info.es
            );
            eprintln!(
                "#   ungapped_cutoff = {} (overrides default min_raw_gapped_score/2 = {})",
                cutoff,
                params.min_raw_gapped_score / 2
            );
        });
    }
}

fn search_db_parallel(
    query: &[u8],
    query_id: &str,
    subject_names: &[String],
    db: &SubjectDb,
    params: &SearchParams,
    matrix: &ScoreMatrix,
    avg_subj_length: u64,
    dump_prelims: Option<&str>,
    prelim_cull: Option<PrelimCullParams>,
    // -evalue reap (NCBI Blast_HSPListReapByEvalue): drop HSPs whose E-value
    // exceeds the threshold, BEFORE masklevel — reaping changes which HSPs
    // compete there.  None when statistics are unavailable (NCBI sentinel
    // mode) or come from a reporting-only source.
    reap: Option<(rmblast_lib::ka_stats::RmStats, f64)>,
) -> Vec<AlignResult> {
    const INITIAL_CHUNK_SIZE: usize = 1_000_000;
    const OVERLAP: usize = 100;

    let full_q_len = query.len().saturating_sub(2);

    // Resolve the ungapped pre-filter cutoff for this search.  NCBI's
    // BlastInitialWordParametersNew derives it from Karlin-Altschul stats when
    // the matrix+gap combo has ALP params; otherwise it falls back to
    // min_raw_gapped_score/2.  Replicate that here and pass it down via params.
    let (ung_cutoff, ka_info) = rmblast_lib::search::ungapped_cutoff(
        &params.matrix_name,
        params.gap_open,
        params.gap_extend,
        params.min_raw_gapped_score,
        full_q_len as u64,
        avg_subj_length,
    );
    log_ungapped_cutoff(params, ung_cutoff, ka_info.as_ref());
    let params = &{
        let mut p = params.clone();
        p.ungapped_cutoff = Some(ung_cutoff);
        p
    };
    let full_q_len_u32 = full_q_len as u32;

    // Mirror NCBI's SplitQuery_CalculateNumChunks + x_ComputeChunkRanges.
    //
    // Step 1: num_chunks = floor(query_length / (initial_chunk_size - overlap))
    // Step 2: if num_chunks <= 1, use a single chunk (no splitting)
    // Step 3: re-adjust chunk_size to distribute load evenly:
    //           chunk_size = (query_length + (num_chunks-1)*overlap) / num_chunks
    //           if num_chunks < chunk_size - overlap: chunk_size++
    // Step 4: generate num_chunks chunks with step = chunk_size - overlap
    let initial_step = INITIAL_CHUNK_SIZE - OVERLAP;
    let num_chunks = if INITIAL_CHUNK_SIZE > OVERLAP { full_q_len / initial_step } else { 0 };

    let (chunk_size, effective_num_chunks) = if num_chunks <= 1 {
        (full_q_len, 1usize)  // no splitting; single chunk covers full query
    } else {
        let mut cs = (full_q_len + (num_chunks - 1) * OVERLAP) / num_chunks;
        if num_chunks < cs.saturating_sub(OVERLAP) { cs += 1; }
        (cs, num_chunks)
    };

    let step = chunk_size.saturating_sub(OVERLAP);

    // Single-chunk path: use original run_gapped_phase pipeline (unchanged).
    if effective_num_chunks <= 1 {
        if prelim_cull.is_some() {
            eprintln!(
                "Note: --prelim-cull applies only to multi-chunk queries (>1 Mbp); \
                 running the faithful pipeline for '{}'.",
                query_id
            );
        }
        let (chunk_lookup, _) = build_query_lookup(query, params);
        let mut results: Vec<AlignResult> = subject_names
            .par_iter()
            .flat_map(|name| {
                let seq = match db.get_full_sequence_blastna(name) {
                    Ok(s) => s,
                    Err(e) => { eprintln!("warning: skipping {}: {}", name, e); return vec![]; }
                };
                let n_mask = db.get_n_mask(name);
                // Shared per-subject RC/packed strands: computed once in the db
                // cache and reused across every query searched against `name`.
                let prep = db.get_prepared(name).ok();
                let mut r = search_with_query_lookup(
                    &chunk_lookup, query, query_id, &seq, name, params, matrix, 0, &n_mask,
                    prep.as_deref(),
                );
                for h in &mut r { h.hsp.q_len = full_q_len_u32; }
                r
            })
            .collect();
        if let Some((stats, expect)) = &reap {
            results.retain(|r| stats.evalue(r.hsp.score) <= *expect);
        }
        // Cross-subject masklevel: mirrors Blast_HSPResultsApplyMasklevel — all subjects
        // combined, sorted by (score DESC, oid DESC), then filtered globally.
        apply_mask_level(&mut results, params.mask_level, subject_names);
        results.sort_unstable_by(|a, b| {
            a.subject_id.cmp(&b.subject_id)
                .then_with(|| a.hsp.strand.as_str().cmp(b.hsp.strand.as_str()))
                .then_with(|| a.hsp.q_start.cmp(&b.hsp.q_start))
        });
        return results;
    }

    // Multi-chunk path.
    // Phase 1+2a (per chunk, parallel over subjects) → cross-chunk merge (sequential).
    // accumulated_prelims[i] = global-coord prelim HSPs for subject_names[i].
    //
    // NCBI applies DUST masking to the full query ONCE (before chunking).  Per-chunk
    // DUST produces different intervals near chunk boundaries (sentinel artifacts), so
    // we replicate NCBI's approach: mask the full query first, then extract masked
    // chunks from it.  build_query_lookup_premask builds the LUT from a pre-masked
    // chunk without re-applying DUST.  Masked chunk is used for seeding only (LUT +
    // word extension + ungapped scoring); Phase 2a and Phase 2b use unmasked query,
    // matching NCBI's mask_at_hash=TRUE behavior for blastn/rmblastn.
    let masked_full_query = mask_query_for_alignment(query, params);
    let mut accumulated_prelims: Vec<Vec<PrelimHsp>> = vec![Vec::new(); subject_names.len()];

    for k in 0..effective_num_chunks {
        let chunk_start = k * step;
        if chunk_start >= full_q_len { break; }
        let chunk_end = (chunk_start + chunk_size).min(full_q_len);
        let chunk_len = chunk_end - chunk_start;
        if chunk_len < params.word_size { break; }

        // Unmasked chunk (original bases) for Phase 2a gapped alignment scoring.
        let mut chunk: Vec<u8> = Vec::with_capacity(chunk_len + 2);
        chunk.push(query[0]);
        chunk.extend_from_slice(&query[1 + chunk_start..1 + chunk_end]);
        chunk.push(*query.last().unwrap_or(&14));

        // Masked chunk (DUST-masked bases) for seeding: LUT + word extension + ungapped scoring.
        let mut masked_chunk: Vec<u8> = Vec::with_capacity(chunk_len + 2);
        masked_chunk.push(masked_full_query[0]);
        masked_chunk.extend_from_slice(&masked_full_query[1 + chunk_start..1 + chunk_end]);
        masked_chunk.push(*masked_full_query.last().unwrap_or(&14));

        let chunk_lookup = build_query_lookup_premask(&masked_chunk, params);

        let chunk_prelims: Vec<Vec<PrelimHsp>> = subject_names
            .par_iter()
            .map(|name| {
                let seq = match db.get_full_sequence_blastna(name) {
                    Ok(s) => s,
                    Err(e) => { eprintln!("warning: skipping {}: {}", name, e); return vec![]; }
                };
                let prep = db.get_prepared(name).ok();
                search_phase2a(&chunk_lookup, &chunk, &masked_chunk, &seq, name, params, matrix, chunk_start as u32, prep.as_deref())
            })
            .collect();

        let split_point = chunk_start as u32;
        if k == 0 {
            for (acc, new) in accumulated_prelims.iter_mut().zip(chunk_prelims.into_iter()) {
                *acc = new;
            }
        } else {
            for (acc, new) in accumulated_prelims.iter_mut().zip(chunk_prelims.into_iter()) {
                merge_chunk_prelims(acc, new, split_point);
            }
        }
    }

    // Optional: dump all merged prelim HSPs before Phase 2b (for phase-attribution analysis).
    if let Some(path) = dump_prelims {
        use std::io::Write as IoWrite;
        let f = std::fs::File::create(path)
            .expect("cannot create dump-prelims file");
        let mut bw = std::io::BufWriter::new(f);
        for (name, prelims) in subject_names.iter().zip(accumulated_prelims.iter()) {
            for p in prelims {
                let strand = if p.strand == rmblast_lib::hits::Strand::Plus { "plus" } else { "minus" };
                writeln!(bw, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    name, strand, p.q_start, p.q_end, p.s_start, p.s_end, p.score, p.q_seed, p.s_seed)
                    .expect("write prelim");
            }
        }
    }

    // Phase 2b: full-query traceback, parallel over subjects.
    // Reverse-complement the full query ONCE and share it (by ref) across all
    // per-subject Phase 2b tasks; previously each parallel task recomputed it,
    // holding one ~full-query copy per concurrent thread.
    let query_rc = rmblast_lib::encoding::revcomp_blastna(query);
    let query_rc = &query_rc[..];

    // Fast-mode prelim cull, round 1: split the merged prelims into a Phase 2b
    // wave-1 set and a culled store (see PrelimCullParams).
    // Interior chunk-edge coordinates for the cull's boundary exemption:
    // prelims truncated at a chunk edge under-represent their final extent.
    let chunk_boundaries: Vec<u32> = if prelim_cull.is_some() {
        let mut b: Vec<u32> = Vec::new();
        for k in 0..effective_num_chunks {
            let cs = k * step;
            if cs >= full_q_len { break; }
            let ce = (cs + chunk_size).min(full_q_len);
            if cs > 0 { b.push(cs as u32); }
            if ce < full_q_len { b.push(ce as u32); }
        }
        b.sort_unstable();
        b
    } else {
        Vec::new()
    };

    let (wave1, culled_store): (Vec<Vec<PrelimHsp>>, Vec<(usize, PrelimHsp)>) =
        match prelim_cull {
            Some(cp) => {
                let keep = cull_prelims(
                    &accumulated_prelims, cp.cull_coverage, cp.cull_margin_pct,
                    &chunk_boundaries,
                );
                let mut w1 = Vec::with_capacity(accumulated_prelims.len());
                let mut culled = Vec::new();
                for (si, (v, k)) in accumulated_prelims.into_iter().zip(keep).enumerate() {
                    let mut kept = Vec::with_capacity(v.len());
                    for (p, kp) in v.into_iter().zip(k) {
                        if kp { kept.push(p); } else { culled.push((si, p)); }
                    }
                    w1.push(kept);
                }
                (w1, culled)
            }
            None => (accumulated_prelims, Vec::new()),
        };

    let phase2b_wave = |per_subject: Vec<Vec<PrelimHsp>>| -> Vec<AlignResult> {
        subject_names
            .par_iter()
            .zip(per_subject.into_par_iter())
            .flat_map(|(name, prelims)| {
                if prelims.is_empty() { return vec![]; }
                let seq = match db.get_full_sequence_blastna(name) {
                    Ok(s) => s,
                    Err(e) => { eprintln!("warning: skipping {}: {}", name, e); return vec![]; }
                };
                let n_mask = db.get_n_mask(name);
                let mut r = run_phase2b(query, query_rc, query_id, &seq, name, prelims, params, matrix, &n_mask);
                for h in &mut r { h.hsp.q_len = full_q_len_u32; }
                r
            })
            .collect()
    };

    let mut all_results = phase2b_wave(wave1);

    // Round 2: re-test culled prelims against the survivors' FINAL query spans
    // and scores; resurrect (and traceback) those no longer dominated.
    if let Some(cp) = prelim_cull {
        if !culled_store.is_empty() {
            let mut survivors: Vec<(u32, u32, i64)> = all_results
                .iter()
                .map(|r| {
                    let (qs, qe) = if r.hsp.strand == Strand::Minus {
                        (r.hsp.q_start + 1, r.hsp.q_end + 1)
                    } else {
                        (r.hsp.q_start, r.hsp.q_end)
                    };
                    (qs, qe, r.hsp.score as i64)
                })
                .collect();
            survivors.sort_unstable();
            let surv_max_span = survivors.iter().map(|s| s.1 - s.0).max().unwrap_or(0);
            let resurrected = resurrect_prelims(
                &culled_store, &survivors, surv_max_span,
                params.mask_level, cp.resurrect_margin_pct, cp.resurrect_slack,
            );
            if !resurrected.is_empty() {
                let mut by_subj: Vec<Vec<PrelimHsp>> = vec![Vec::new(); subject_names.len()];
                for &ci in &resurrected {
                    let (si, ref p) = culled_store[ci];
                    by_subj[si].push(p.clone());
                }
                let mut r2 = phase2b_wave(by_subj);
                all_results.append(&mut r2);
            }
        }
    }

    if let Some((stats, expect)) = &reap {
        all_results.retain(|r| stats.evalue(r.hsp.score) <= *expect);
    }
    // Cross-subject masklevel: mirrors Blast_HSPResultsApplyMasklevel — all subjects
    // combined, sorted by (score DESC, oid DESC), then filtered globally.
    apply_mask_level(&mut all_results, params.mask_level, subject_names);
    all_results.sort_unstable_by(|a, b| {
        a.subject_id.cmp(&b.subject_id)
            .then_with(|| a.hsp.strand.as_str().cmp(b.hsp.strand.as_str()))
            .then_with(|| a.hsp.q_start.cmp(&b.hsp.q_start))
    });
    all_results
}

fn write_all<W: Write>(
    w: &mut W,
    results: &[AlignResult],
    fields: &[OutField],
    ka: Option<&RmStats>,
) -> Result<()> {
    for r in results {
        write_tabular(w, r, '\t', fields, ka).context("writing output")?;
    }
    Ok(())
}

fn load_gilist(path: Option<&str>) -> Result<Option<HashSet<String>>> {
    let p = match path {
        None => return Ok(None),
        Some(p) => p,
    };
    let f = std::fs::File::open(p).with_context(|| format!("opening gilist '{}'", p))?;
    let mut set: HashSet<String> = HashSet::new();
    for line in std::io::BufReader::new(f).lines() {
        let entry = match line {
            Ok(l) => l.trim().to_string(),
            Err(_) => continue,
        };
        if entry.is_empty() {
            continue;
        }
        // A bare GI number also matches the sequence named "gi|<n>".  This is
        // the form RepeatModeler writes (RepeatModeler:1851 emits
        // `($startGID+1) .. $sampleDBSize`, i.e. plain integers) while its
        // databases name sequences "gi|N"; without this the all-vs-other
        // batches would filter every subject out and silently return no hits.
        if entry.chars().all(|c| c.is_ascii_digit()) {
            set.insert(format!("gi|{}", entry));
        }
        set.insert(entry);
    }
    Ok(Some(set))
}

#[cfg(test)]
mod gilist_tests {
    use super::load_gilist;
    use std::io::Write;

    /// RepeatModeler writes bare integers; the databases name sequences "gi|N".
    /// Both spellings must select the subject (see wrappers/blastdb_aliastool).
    #[test]
    fn bare_numbers_and_gi_pipe_both_match() {
        let dir = std::env::temp_dir().join("rmblastn_gilist_test");
        std::fs::create_dir_all(&dir).unwrap();

        let bare = dir.join("bare.txt");
        writeln!(std::fs::File::create(&bare).unwrap(), "2\n3\n").unwrap();
        let set = load_gilist(Some(bare.to_str().unwrap())).unwrap().unwrap();
        assert!(set.contains("gi|2") && set.contains("gi|3"));
        assert!(set.contains("2"), "the raw spelling must still match");
        assert!(!set.contains("gi|1"));

        let piped = dir.join("piped.txt");
        writeln!(std::fs::File::create(&piped).unwrap(), "gi|2\ngi|3\n").unwrap();
        let set = load_gilist(Some(piped.to_str().unwrap())).unwrap().unwrap();
        assert!(set.contains("gi|2") && set.contains("gi|3"));

        assert!(load_gilist(None).unwrap().is_none());
        std::fs::remove_dir_all(&dir).ok();
    }
}
