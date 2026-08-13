# rmblast-stats — E-values & bit scores for the rmblastn Rust port

**INTEGRATED (2026-08-12):** now a subdirectory of the `rmblast/` workspace;
`rmstats` is a workspace member and a dependency of `rmblast-lib`.  The
integration lives in `rmblast-lib/src/ka_stats.rs` (source hierarchy:
baked table → `-matrix_*` CLI → `# KARLIN` matrix comment → deterministic
ALP fit → sentinel), `output.rs` (`evalue`/`bitscore` outfmt fields with
NCBI tabular rendering), and `main.rs` (lazy resolution + E-value cutoff).
The ALP sources are vendored in `ALP_1.98_LIB/` and built by the default-on
`alp-fit` feature.  Two runtime modes (see `ka_stats.rs` module docs):
default = full hierarchy for reporting, no E-value culling unless
`--evalue` is given; `--ncbi-compat` = exact NCBI 2.17.1 emulation
(hardcoded-table/CLI sources only, 30p53g placeholder wart, sentinel
rendering, and the silent E>10 reap 2.17.1 performs for table matrices) —
required for 2.17.1-parity benchmarks.  The sections below are the
original hand-off notes.

Everything needed to add NCBI-faithful E-values and bit scores to the Rust
rmblastn:

```
rmblast-stats/
├── alp-fit/            C++ fitter: RepeatMasker matrix -> ALP Gumbel fit
├── fits/               Per-matrix fit outputs (.fit) + full ALP state (.par)
│   └── SUMMARY.tsv     All 42 fits, one row per matrix
├── gen_tables.py       Regenerates rmstats/src/tables_data.rs from
│                       blast_stat.c + fits/ (no hand transcription)
├── run_fits.sh         Reruns all fits at canonical gap costs
├── ncbi_table_update.patch   C-side sync for 30p53g + comparison.matrix
└── rmstats/            The Rust crate (see below)
```

## Validation status (2026-08-12)

* `alp-fit` reproduces all **40** fitted `rmblast_*_values` entries in the
  RMBlast-patched `blast_stat.c` **exactly (10 decimal places)** with ALP
  defaults (eps_lambda 0.01, eps_K 0.05, seed 1, iad=false). The C table's
  `30p53g` entry is an explicit placeholder; `comparison.matrix` is absent.
  Both now have real fits (see `ncbi_table_update.patch`).
* `rmstats` cross-validated against the patched NCBI `rmblastn`
  (`BLAST_KA_DEBUG=1`) on human-1mb × shortlib:
  * **Mode 1** (baked table, `20p41g.matrix` @25/5): 6558/6558 HSPs match
    (E-value within %.6g print precision, bits within %.4f).
  * **Mode 2** (CLI `-matrix_lambda/k/alpha/beta`, comparison.matrix @20/5):
    6638/6638 HSPs match, including the beta=-26.84 length adjustment.
  * **Mode 3** (sentinel): NCBI prints `evalue=1.0 bitscore=0.0` for every
    hit; `rmstats` reproduces exactly.
* 11 unit tests + 1 FFI test (`cargo test`, `cargo test --features alp-fit`),
  including Altschul's worked example for the ungapped Lambda/H/K machinery.

## The rmstats crate

Faithful ports (function-by-function, from
`ncbi-blast-2.17.0+-src/c++/src/algo/blast/core/`):

| rmstats module | NCBI source |
|---|---|
| `tables_data` (generated) | `blastn_values_*`, `rmblast_*_values` (blast_stat.c) |
| `rmblast_tables` | `Blast_KarlinBlkGappedLoadFromTables`, three-mode hierarchy of `Blast_ScoreBlkKbpGappedCalc` (blast_setup.c:89-126) |
| `nucl_tables` | `s_GetNuclValuesArray`, `s_SplitArrayOf8`, `s_AdjustGapParametersByGcd`, `Blast_KarlinBlkNuclGappedCalc`, `Blast_GetNuclAlphaBeta`, `s_GetUngappedBeta` |
| `karlin` | ungapped KA: `Blast_ResFreqStdComp/String`, `BlastScoreFreqCalc`, `Blast_KarlinLambdaNR`, `BlastKarlinLtoH`, `BlastKarlinLHtoK`, `Blast_KarlinBlkUngappedCalc` (blastna path) |
| `length_adjust` | `BLAST_ComputeLengthAdjustment` (blast_stat.c:5555) |
| `eff_lengths` | `BLAST_CalcEffLengths` blastn arm incl. -RMH- alpha/beta hierarchy and sentinel skip (blast_setup.c:733-906) |
| `evalue` | `BLAST_KarlinStoE_simple`, `BlastKarlinEtoS_simple`, `BLAST_GapDecayDivisor`, `BLAST_Cutoffs`, sentinel semantics of `Blast_HSPListGetEvalues`/`GetBitScores` (blast_hits.c) |
| `alp` (feature `alp-fit`) | LAST-style startup fitting via ALP FFI |

### High-level API

```rust
use rmstats::{RmStats, MatrixCliOverrides};

// Once per query x database (both blastn strand contexts share this):
let stats = RmStats::new_custom_matrix(
    "20p41g.matrix",          // basename passed to -matrix (case-insensitive)
    25, 5,                    // -gapopen / -gapextend
    &MatrixCliOverrides::default(), // or -matrix_lambda/k/alpha/beta values
    query_length as i32,      // this query's length
    total_db_letters as i64,  // sum of subject lengths (before masking)
    n_db_seqs as i32,
    0,                        // -searchsp override; 0 = unset
    None,                     // kbp_std; see note below
);

// Per HSP, using the FINAL score (i.e. after complexity adjustment):
let evalue = stats.evalue(hsp.score);      // 1.0 when stats unavailable
let bits   = stats.bit_score(hsp.score);   // 0.0 when stats unavailable
```

Semantics guaranteed to match NCBI rmblastn:

* **Mode 1**: matrix+gap-costs found in the baked table → table Lambda/K/H,
  `alpha = Lambda/H`, `beta = 0` for the length adjustment.
* **Mode 2**: table miss but `-matrix_lambda/-matrix_k/-matrix_alpha` all
  supplied → those values, `H = lambda/alpha`, and `-matrix_beta` (default
  0.0) feeds the length adjustment.
* **Mode 3**: neither → sentinel; every hit reports `evalue = 1.0`,
  `bit_score = 0.0` (do **not** print stats in that case, or print them
  as NCBI does — it prints the 1.0/0.0 values in tabular output).
* `round_down` is always false on the custom-matrix path (it belongs to the
  reward/penalty path, also ported in `nucl_tables` for completeness).
* Bit score always uses the raw (un-rounded) score; the even-score rounding
  applies to E-values only, and only when `round_down` (SB-2303).

`kbp_std` (the ungapped per-context Karlin block, computable with
`rmstats::kbp_ungapped_calc_blastna(&matrix16, query_context_bytes)`) is only
consulted by the Mode-2 *alpha/beta fallback branch* of
`BLAST_CalcEffLengths`, which is unreachable in `read_in_matrix` mode: when
CLI params are set they take the first branch, when the table hit is valid the
second, and when neither the sentinel skip fires before the value is used.
Passing `None` is faithful for all rmblastn configurations; the machinery is
ported anyway to keep the crate usable for plain blastn statistics.

## Integrating into the rmblast port (for the other session)

1. **Dependency**: `rmstats = { path = "../rmblast-stats/rmstats" }` in
   `rmblast-lib/Cargo.toml` (default features; no C++ toolchain needed).
2. **Thread the DB totals**: `main.rs` already computes `total_db_letters`
   and `n_db_seqs` (main.rs ~300-330) but doesn't pass them into
   `search_db_parallel`. Either extend `SearchParams` (options.rs) with
   `db_total_len: u64` / `db_num_seqs: u32`, or construct the `RmStats`
   in `main.rs` and pass it down alongside `avg_subj_length`.
3. **Construct once per query**: matrix name must be the basename as given
   to `-matrix`. Note the CLI additions if full parity is wanted:
   `-matrix_lambda`, `-matrix_k`, `-matrix_alpha`, `-matrix_beta`
   (cmdline_flags.cpp:81-84) — without them Mode 2 is unreachable, which is
   fine for a first pass.
4. **Fill per-HSP fields**: add `evalue: f64` / `bit_score: f64` to `Hsp`
   (hits.rs) or `AlignResult` (output.rs); fill at the two construction
   sites — `engine.rs:1694` (`run_phase2b`) and `engine.rs:2383`
   (`run_gapped_phase`) — or once in the output layer, since
   evalue/bit-score depend only on `hsp.score` and the per-query `RmStats`.
   IMPORTANT: use the score *after* complexity adjustment (NCBI computes
   E-values at traceback time from the adjusted score;
   blast_traceback.c:544-593 runs before Blast_HSPListGetEvalues).
5. **Output**: add `Evalue` / `Bitscore` to `OutField` + `from_str` +
   `write_tabular` arms (output.rs). To byte-match NCBI's rendering, port
   `CAlignFormatUtil::GetScoreString`
   (`src/objtools/align_format/align_format_util.cpp:940`); exact rules:
   * e-value: `< 1e-180` → `"0.0"`; `< 1e-99` → `%2.0le`; `< 0.0009` →
     `%3.0le`; `< 0.1` → `%4.3lf`; `< 1.0` → `%3.2lf`; `< 10.0` → `%2.1lf`;
     else `%2.0lf`. (The Mode-3 sentinel `evalue = 1.0` therefore prints as
     `1.0`.)
   * bit score: `> 99999` → `%5.3le`; `> 99.9` → `%3.0ld` (cast to long,
     i.e. truncated); else `%4.1lf` (note the space-padded width — the
     sentinel `0.0` prints as `" 0.0"`, as observed in the validation run).
6. **Cutoff reuse**: `ka_cutoff.rs` keeps its own transcription of
   lambda/K for the seeding cutoff. It can now delegate to
   `rmstats::rmblast_tables::find_matrix_entry` /
   `karlin_blk_gapped_load_from_tables` to avoid the duplicate table
   (the rmstats table also fixes 30p53g and adds comparison.matrix — the
   ka_cutoff copy lacks comparison.matrix entirely).

### Divergence from NCBI (intentional, documented)

The rmstats baked table differs from the C tree in exactly two entries:
real ALP values for `30p53g.matrix` (C has a placeholder) and a new
`comparison.matrix` @20/5 entry (C has none, so NCBI reports 1.0/0.0 for
the canonical RepeatModeler configuration!). Apply
`ncbi_table_update.patch` to the C tree and rebuild when you want both
sides identical — do not do that while comparison runs depend on the
current NCBI binary.

## Startup-time fitting (feature `alp-fit`)

```rust
use rmstats::alp::{fit_gumbel_blastna, AlpFitOptions};

let fit = fit_gumbel_blastna(&matrix16, &freqs16, gap_open, gap_extend,
                             &AlpFitOptions::default())?;
let kbp = fit.to_karlin_blk();   // NCBI mapping: alpha=a_J+a_I, H=lambda/alpha
```

* Builds the vendored ALP sources via `cc` (`ALP_SRC_DIR` overrides the
  default `../ALP_1.98_LIB/cpp` (vendored)).
* `AlpFitOptions::default()` = the exact settings of the baked fits, so a
  known matrix refits to bit-identical values (regression-tested).
* `deterministic: true` gives LAST's reproducible fixed-sample mode
  (`set_gapped_computation_parameters_simplified(max_time)` + retry loop
  `temperature = 1.07 + 0.01*attempt`, ≤21 attempts).
* **Not thread-safe** (ALP uses process-global RNG state): fit at startup,
  single-threaded, before spawning search threads. Typical fit time for
  these 4x4 DNA matrices: well under a second.
* Suggested engine flow: table lookup → CLI params → *ALP fit at startup*
  (new Mode 2.5) → sentinel. That keeps NCBI behavior for everything NCBI
  can do, and adds statistics for arbitrary custom matrices.

## Refitting / extending the table

```sh
cd rmblast-stats
(cd alp-fit && make)        # builds ALP lib + fitter
./run_fits.sh               # all 42 canonical fits -> fits/
./alp-fit/alp_fit -matrix path/to/foo.matrix -gapopen 20 -gapextend 5 \
    -max_time 10 -out fits/foo.par > fits/foo.fit   # one-off fit
python3 gen_tables.py       # regenerate rmstats/src/tables_data.rs
(cd rmstats && cargo test)  # verify
```
