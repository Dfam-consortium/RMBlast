# verify_wrappers/

Drop-in replacements for `rmblastn`, `makeblastdb`, and `blastdbcmd` that run
**both** the real NCBI engine and the Rust port serially, compare their results,
and save a debug bundle whenever they disagree.

These share the same file names as the real NCBI tools (and as the plain
`../wrappers/` set) so they slot into any caller that finds them first on `PATH`
or via `RMBLAST_DIR`. They must live in a *separate* directory from the plain
single-engine wrappers precisely because the names collide.

## How it differs from `../wrappers/`

| tool        | `../wrappers/` (plain)                       | `verify_wrappers/` (this dir)                                  |
|-------------|----------------------------------------------|----------------------------------------------------------------|
| `rmblastn`  | runs Rust only (db -> `.2bit`)               | runs **NCBI + Rust**, diffs, bundles mismatches; emits Rust    |
| `makeblastdb`| builds `.2bit` + dummy index files          | builds the **real** NCBI BLAST index (FASTA kept for Rust)     |
| `blastdbcmd`| reads `.2bit` via `twoBitInfo`               | defers to the **real** NCBI `blastdbcmd`                        |

Key point: the Rust engine is pointed at the **FASTA file** as its database,
never a `.2bit` — only the FASTA reproduces NCBI's ambiguity-code handling.

Because of that, both `makeblastdb` wrappers copy the `-in` FASTA to the `-out`
base path whenever the two differ, so a FASTA is always sitting where
`rmblastn` looks for `-db <base>`.  RepeatMasker builds its databases in place
(`-in x.fa -out x.fa`) and is unaffected, but RepeatModeler's
`RepeatUtil::ncbiMaskDatabaseNativeMT` builds
`makeblastdb -out <round-N>/tmpConsDB-1 -in <root>/consensi.fa` and then
searches `-db <round-N>/tmpConsDB-1`; without the copy the `rmblastn` wrappers
abort with `database file not found`.  It is a copy rather than a symlink
because `consensi.fa` keeps growing between rounds and the database has to stay
paired with the sequences it was built from.

## `--ncbi-compat` is forced on the Rust side

This wrapper appends `--ncbi-compat` to every Rust invocation (both the
comparison run and the caller-format run). It is *required* for a meaningful
comparison:

NCBI 2.17.1 silently drops HSPs with E-value > 10 (`Blast_HSPListReapByEvalue`,
applied before masklevel) whenever it has valid Karlin-Altschul parameters —
i.e. for every matrix in its hardcoded table, which is all the RepeatMasker
`p##g` matrices. The Rust port applies **no** E-value cutoff by default, so
without the flag it legitimately keeps a few marginal hits (E ≈ 10–30, typically
scores just above `-min_raw_gapped_score`) that NCBI discards, and the wrapper
would report them as mismatches. Measured example: human-1mb × shortlib with
`20p41g.matrix` @25/5 — 691 hits without the flag vs NCBI's 686.

`comparison.matrix` is absent from NCBI's table, so both engines run in sentinel
mode there and cull nothing; those runs are unaffected either way.

The flag is also stripped from the **NCBI** argument vector if a caller supplies
it (NCBI has no such option and would abort with a usage error), and it is not
added twice if the caller already passed it in either spelling
(`-ncbi_compat` / `--ncbi-compat`).

## Usage

Put this directory ahead of the real NCBI tools on `PATH` (or set it as the
RepeatMasker/RepeatModeler `RMBLAST_DIR`):

```sh
export PATH=<rmblast-root>/verify_wrappers:$PATH
# callers then invoke `makeblastdb`/`rmblastn`/`blastdbcmd`/`dustmasker` unchanged
```

`<rmblast-root>` is either an installed release (`/usr/local/rmblast-<version>`)
or the `rmblast/` directory of a source checkout — the wrappers find the Rust
binaries relative to themselves, so neither needs editing.

The caller must still set `BLASTMAT` to the matrix directory (as usual for NCBI
rmblastn); both engines inherit it and resolve the matrix the same way.

## What gets compared

The cross-check **always compares tabular alignment data** (`-outfmt 6`), never
the caller's possibly-pairwise output. Pairwise/default output carries headers
(engine version, DB path, build dates, number formatting) that differ
cosmetically between the two engines and are *not* what we are verifying —
diffing them produced spurious "failures" with no actual alignment difference.

- If the caller already requested tabular output (`-outfmt "6 ..."`), that output
  is compared directly (and the single Rust run also serves as the caller's
  output).
- Otherwise (pairwise or any non-tabular format) the wrapper runs a dedicated
  canonical `-outfmt 6` pass on each engine for the comparison, plus a separate
  caller-format Rust run for the output handed back to the caller. The canonical
  field list is the same one the project's 16-combo validation scripts use
  (score, positions, strand, `qseq`/`sseq`, etc.).

## `-version` / `-help`

Informational invocations do **not** trigger a cross-check. `rmblastn -version`
(or `-help`/`-h`) simply runs the primary (Rust) engine and reports its output.
Note this reports the Rust port's version string (`rmblastn 0.1.0`); if a caller
parses the version and needs an NCBI-style string, change the `-version` branch
near the top of `rmblastn` to exec `$NCBI_RMBLASTN --version` instead.

## Output and exit status

- The **Rust** output is what gets handed back to the caller (this directory is
  a drop-in for the Rust deployment); the wrapper exits with the Rust engine's
  return code.
- A one-line `# verify rmblastn: Rust MATCHES NCBI (...)` note is printed to
  STDERR on a clean match.
- On any mismatch (different output lines, the same lines in a different
  order, or either engine exiting non-zero) the wrapper prints a preview to
  STDERR and saves a full bundle.  The log status is `MISMATCH`, `ORDER`, or
  `ERROR`.

## Mismatch bundles

Written under `$RMBLAST_VERIFY_DIR` (default `./rmblast_verify_failures`,
relative to the caller's working directory). Each failure gets its own
subdirectory containing:

- the query FASTA, the database FASTA + its BLAST index files, the matrix
- `ncbi.tab` / `rust.tab` — raw outputs, and `*.stderr` captures
- `diff.txt` — `<` lines are NCBI-only, `>` lines are Rust-only; for an
  order-only mismatch, a plain diff of the two raw outputs
- `commands.txt` — the exact invocations and `BLASTMAT`
- `reproduce.sh` — reruns both engines against the copied-in files and diffs

## Configuration

Environment variables (no editing required):

- `NCBI_DIR` — where the real NCBI tools live.
  Default `/usr/local/rmblast-2.17.1/bin`.
- `RMBLAST_BIN_DIR` — where the Rust binaries live.  By default the wrappers
  resolve `../bin/<name>` (installed release) and then
  `../target/release/<name>` (source checkout), relative to the wrapper itself
  with symlinks followed.
- `RMBLAST_VERIFY_DIR` — overrides the bundle directory.
- `RMBLAST_VERIFY_DUSTMASKER` — `ncbi` (default) or `rust`.  Only `rmblastn` is
  A/B-compared here; `makeblastdb`, `blastdbcmd` and `dustmasker` defer to the
  real NCBI binaries so that any reported diff is attributable to `rmblastn`
  alone.  Set this to `rust` to exercise the port's dustmasker instead.
- `blastdb_aliastool` is the same pass-through copy as in `../wrappers/`: it
  writes the GI list out as plain text, which both engines accept.
- `-version` / `--version` / `-v` is answered by the wrapper, not passed
  through.  The Rust binary prints clap's plain `rmblastn 3.0.3`; both
  `rmblastn` wrappers reshape that into NCBI's form:

  ```
  rmblastn: 3.0.3+
   Package: rmblast 3.0.3, Rust port
  ```

  The colon and the trailing `+` are load-bearing.  Every caller-side version
  probe needs them: `NCBIBlastSearchEngine::setPathToEngine` matches
  `/(\d+)\.(\d+)\.(\d+)\+/` before setting `hasTabFormat`, and without it the
  caller omits `-outfmt "6 score perc_sub ..."` and silently falls back to the
  legacy pairwise parser (losing kdiv / cpg_kdiv / transi / transv /
  cpg_sites).  `RepModelConfig.pm` and `RepeatModeler/configure` additionally
  need the colon.  The reshaping lives in the wrappers on purpose: `+` is an
  NCBI packaging convention, not this project's version, so the binary keeps
  reporting its own `rmblastn <ver>`.  The number is always read back from the
  binary, never hardcoded.
