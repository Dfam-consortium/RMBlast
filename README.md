# RMBlast — Rust port

A Rust reimplementation of `rmblastn`, the RepeatMasker nucleotide alignment
program.

***THIS IS AN EXPERIMENTAL PORT and should not be used with the current releases
   of RepeatMasker or RepeatModeler in production settings at this time. For the
   current supported version please see: http://www.repeatmasker.org/rmblast***

Given the same FASTA formatted database, it reproduces NCBI RMBlast 2.17.1 output 
faithfully. In addition, we have switched out the NCBI database format for the UCSC
twobit format (currently without IUB support).  Finally, in this port E-value
culling is off by default unless you pass the `--ncbi-compat` flag.

The package also ships a `dustmasker` port that is a drop-in for the NCBI application
of the same name, and wrapper scripts that support testing of this port with 
the current version of RepeatMasker (nucleotide tools only).  RepeatModeler can als
be tested using the wrapper scripts, however, there is a dependency on blastx in
the last classification step that this package does not provide. 

## What this is a port of

The algorithms here were translated from two sources, both in the public domain:

- NCBI BLAST 2.17.0 C++ toolkit: the lookup table, ungapped and gapped
  extension, traceback, DUST filtering, and Karlin-Altschul statistics.
- The rmblastn extensions in RMBlast 2.17.1, by Robert Hubley
  at the Institute for Systems Biology. These are what separate
  `rmblastn` from stock `blastn`: custom scoring matrices without
  Karlin-Altschul statistics, cross_match-style complexity-adjusted scoring, and
  cross_match-style masklevel filtering.

We ported from the C sources rather than from NCBI's output. Where NCBI has a
quirk that changes results, the port reproduces it on purpose. 

## Requirements

- Rust stable (2021 edition).
- A C++ toolchain. The `alp-fit` feature is on by default and compiles the
  vendored ALP 1.98 sources, which supply Gumbel parameters for matrices absent
  from NCBI's hardcoded table. Build `rmblast-lib` with `--no-default-features`
  to skip it; E-value statistics then fall back to the sentinel.
- For the wrapper scripts only: UCSC's `faToTwoBit`, `twoBitInfo`, and
  `twoBitToFa`, found in `$UCSCTOOLS_DIR`, then on `PATH`, then in
  `/usr/local/ucscTools`. `makeblastdb` builds a usable database without
  `faToTwoBit` (searches read the FASTA); `blastdbcmd` needs the other two only
  for a `.2bit` database.

## Build and install

```sh
make                      # cargo build --release
make test                 # cargo test --release
make install PREFIX=/opt  # -> /opt/rmblast-<version>/
```

The default build leaves out the stage counters and dump hooks
(`RMBLAST_DP_STATS`, `RMBLAST_DEBUG_COUNTERS`, `BLAST_DUMP_SEEDS`,
`RMBLAST_DUMP_IMPROVE`, `RMBLAST_DUMP_DP_ROWS`). To get them, build with the `diagnostics` feature:

```sh
cargo build --release --features diagnostics
```

`make install` produces a self-contained, relocatable tree:

```
rmblast-<version>/
├── bin/{rmblastn,dustmasker}
├── wrappers/{rmblastn,makeblastdb,blastdbcmd,blastdb_aliastool,dustmasker}
├── verify_wrappers/{...,README.md}
├── README.md
└── LICENSE
```

Both wrapper directories find the binaries relative to themselves, so nothing
needs editing after unpacking. `make dist` packages that same tree as a
versioned tarball, which is what the tagged GitHub Actions release publishes;
it stages the install under `./stage` and never writes to `PREFIX`.

`rmblastn --version` names the commit it was built from, for example
`rmblastn 3.1.0 (01cc4cc)`, with `-dirty` added when the tree had uncommitted
changes. The wrappers print NCBI's `rmblastn: 3.1.0+`, which is what
RepeatMasker and RepeatModeler parse.

To bundle scoring matrices alongside the binaries, point `MATRIX_SRC` at a
directory of `*.matrix` files: `make install PREFIX=/opt MATRIX_SRC=../matrices`.

## Running it directly

```sh
rmblastn --query genome.fa --db library.fa \
         --matrix comparison.matrix \
         --gapopen 20 --gapextend 5 --word-size 7 \
         --complexity-adjust --mask-level 80 --dust no \
         --outfmt "6 score perc_sub qseqid qstart qend sseqid sstart send"
```

`--help` lists every option. A few worth knowing about:

| Option | Notes |
|---|---|
| `--db` | FASTA or UCSC `.2bit`, autodetected. Prefer FASTA. `.2bit` cannot store non-N IUB codes, so `faToTwoBit` folds them all to N and output diverges from NCBI on any sequence containing them. |
| `--matrix` | Required. Resolved the way NCBI resolves it, including `$BLASTMAT/nt` (see below). |
| `--mask-level` | Defaults to `-1` (off), matching NCBI. RepeatMasker passes `80` explicitly; RepeatModeler does not pass it at all. |
| `--min-raw-gapped-score` | A cutoff for the preliminary gapped stage, not a floor on the reported score. Traceback can score lower, and that lower score is what gets reported, as in NCBI. Post-filter if you need a hard floor. |
| `--evalue` / `--ncbi-compat` | By default nothing is culled on E-value. NCBI 2.17.1 silently drops hits with E > 10 whenever the matrix is in its baked table, and `--ncbi-compat` reproduces that. Pass it for parity benchmarks. |
| `--outfmt` | `0` for pairwise, or `6` followed by field names. Other formats (`7`, `10`, ...), a bare `6` (NCBI's default 12 columns) and field names not listed below are refused. |
| `--mt-mode` | Which axis `--num-threads` splits on. See below. |

Tabular fields: `score`, `perc_sub`, `perc_query_gap`, `perc_db_gap`, `qseqid`,
`qstart`, `qend`, `qlen`, `sstrand`, `sseqid`, `sstart`, `send`, `slen`, `kdiv`,
`cpg_kdiv`, `transi`, `transv`, `cpg_sites`, `qseq`, `sseq`, `evalue`,
`bitscore`.

NCBI's argument spellings work as-is, so `-word_size 7` is accepted alongside
`--word-size 7`. Single-character flags cannot be bundled.

### Choosing `--mt-mode`

`--mt-mode` picks the axis threads are split across, and the right choice
depends on the shape of the search:

- `0` (default) parallelises across database sequences. This suits the standard
  RepeatMasker shape, a genome query against a many-sequence TE library, where
  it scales close to linearly.
- `1` parallelises across query sequences. Use it for the reverse shape, many TE
  queries against a single-sequence genome database, where mode 0 has nothing to
  split across.

Output is bit-identical across both modes and every thread count.

### Finding the scoring matrix

`--matrix comparison.matrix` names a file, and the port looks for it where NCBI
does: `$NCBI_DATA_PATH`, then `$BLASTMAT` and `$BLASTMAT/nt`, then `./data/`.
RepeatMasker sets `BLASTMAT` to the directory holding the matrix, so this
usually takes care of itself.

An absolute path also works, which stock `rmblastn` rejects. That is the one
place the lookup diverges, and the port tries it last, after every directory
NCBI would have searched, so it can never shadow one of them.
`find_matrix_path` in `rmblastn/src/main.rs` has the full order and the NCBI
quirks inside it.

## Using it with RepeatMasker and RepeatModeler

`wrappers/` holds stand-ins for the NCBI tools those pipelines call: `rmblastn`,
`makeblastdb`, `blastdbcmd`, `blastdb_aliastool`, and `dustmasker`. They
translate NCBI-style arguments, build the database representation the port
expects, and answer the database queries the callers make. Point RepeatMasker or
RepeatModeler at that directory as its rmblast location, or put it first on
`PATH`:

```sh
export PATH=/opt/rmblast-<version>/wrappers:$PATH
```

The `makeblastdb` here never calls NCBI's. It converts the FASTA to `.2bit`,
leaves a FASTA at the database base path for the engine to read, and writes
placeholder index files so that callers checking for a BLAST database are
satisfied. Every file is written under a temporary name and renamed into place,
so a search during a rebuild never reads a half-written database. It accepts the
options RepeatMasker and RepeatModeler pass (`-in`, `-out`, `-dbtype nucl`,
`-parse_seqids`, `-blastdb_version`, `-title`). It refuses `-dbtype prot` (the
port is nucleotide only), `-mask_data` (the engine would not apply the masks)
and any other option.

The `rmblastn` wrapper passes every argument through unchanged, values
included, except that it drops `-num_alignments`, points `-db` at the FASTA,
and prints `-version` in NCBI's `rmblastn: X.Y.Z+` shape, which RepeatMasker
and RepeatModeler parse. `bash tests/wrappers.sh` tests all of the wrappers,
and compares against NCBI where it is installed.

### verify_wrappers/

`verify_wrappers/` holds a parallel set that runs both engines on every
invocation, diffs the results, and saves a bundle of inputs and commands
whenever they disagree. It exercises the port against real RepeatMasker and
RepeatModeler runs instead of a fixed test set, which catches what the suites
structurally cannot: the four combo suites pin `-mask_level` on every
invocation, so none of them could ever have caught the wrong default that a
RepeatModeler run exposed immediately. Using it requires NCBI RMBlast installed
as well. See `verify_wrappers/README.md`.

Bundles land in `./rmblast_verify_failures` by default; `RMBLAST_VERIFY_DIR`
moves them. Every invocation also appends a row to `verify_log.tsv` with
timings, hit counts, and a PASS/MISMATCH/ERROR status.

## Known differences from NCBI

On the same FASTA database the port matches NCBI's results. Below are the
exceptions, and the places where the tools around it stop short of NCBI's:

- `.2bit` databases lose ambiguity codes. `faToTwoBit` stores only A/C/G/T
  plus an N-block table, so R/Y/K/M/W/S/B/D/H/V all collapse to N and the seeds
  differ from NCBI's on any sequence containing them. Supply the database as
  FASTA for full IUB fidelity. The port prints a note to stderr when it opens a
  `.2bit`.
- E-value culling is off unless you ask for it. NCBI 2.17.1 silently drops hits
  with E > 10 whenever the matrix is in its baked table; the port keeps them
  unless you pass `--evalue` or `--ncbi-compat`. On a matrix absent from that
  table, such as `comparison.matrix`, both engines cull nothing.
- Pairwise output (`-outfmt 0`) is NCBI's line for line, for a database built
  with `-parse_seqids` (as RepeatMasker and RepeatModeler build theirs), with
  two exceptions. The "Database:" title is the `-db` path given to rmblastn,
  where NCBI prints the `-in` path given to makeblastdb. "Posted date" is the
  database file's modification time in UTC, where NCBI prints the build time
  in local time. Without `-parse_seqids`, NCBI wraps a long subject defline
  slightly differently, and an id-only defline has no trailing space.
- The bare binary rejects `-num_alignments`; `wrappers/rmblastn` drops it.
- NCBI 2.17.1 rejects `-xdrop_ungap` ("Attempt to cast to a wrong (Integer)
  type") unless `-min_raw_gapped_score` is also given; the port accepts it on
  its own.
- Three small fidelity gaps remain, all far from RepeatMasker's parameters
  (which always set a minimum score of 200 or more):
  - At `-mask_level 0`, which HSPs a touching neighbour masks can differ,
    because NCBI's choice depends on its interval-tree layout.
  - With `identity.matrix` and no `--min-raw-gapped-score`, the port keeps a
    minus-strand seed that NCBI's exact extension drops, so a few hits end a
    little differently.
  - With no minimum score at all, both programs report single-base HSPs; at
    `-word_size 4`, one of those can carry the other strand label (1 line in
    about 679,000 on the test fixtures).
- Not ported, with no known effect on results:
  - NCBI's small lookup table (`eSmallNaLookupTable`) and its specialised
    megablast scanners: the port uses the standard table and the general
    scanner for every width. Results match on every comparison suite; the
    `identity.matrix` seed difference above may come from here.
  - `AdjustSubjectRange`, which narrows the subject range NCBI extends over.
- Not supported: `-searchsp` (the effective search space is always computed),
  gzip-compressed input, and these NCBI spellings, all refused with an error:
  negative gap costs, `-mt_mode 2`, and `-dust 'level window linker'`.
- The built-in Karlin-Altschul table is looked up by the matrix file's name
  and gap costs, as NCBI does; the matrix contents are not checked against it.
  The same matrix under another name falls back to a startup ALP fit, whose
  lambda can differ by about 1%. When E-values or bit scores are computed,
  stderr says which source was used.
- The whole database stays resident. NCBI streams database volumes; the port
  loads and decodes all of it, and caches the reverse complement and packed
  strands of every subject it touches at about 1.5 times that subject's size.
  Both are negligible for a typical TE library of a megabyte or so, but they
  grow with database size, so a genome-sized subject needs more memory here than
  under NCBI.
- `dustmasker` covers part of the NCBI interface. It handles `-in`, `-out`,
  `-window`, `-level`, `-linker`, `-infmt fasta`, `-outfmt` of
  `interval`/`fasta`/`acclist`, and `-hard_masking`. Everything else, including
  `-infmt blastdb`, the ASN.1 and XML `-outfmt` variants, and `-parse_seqids`,
  exits non-zero with a "not ported" message rather than emitting near-miss
  output. Deflines are printed as NCBI prints them, through a port of NCBI's
  title cleanup: spaces collapsed, the line cut at a tab, trailing periods
  stripped, HTML entities decoded. Like NCBI, a final non-ASCII character is cut
  to its first byte. The output matched NCBI dustmasker 2.17.0 byte for byte
  across all three output formats, on 13,000 random deflines and on the
  repository's own tests (`dustmasker/tests/`).
  Where NCBI is lenient, this port refuses: an empty record (NCBI also aborts),
  any character that is not a nucleotide code (NCBI drops it), and a negative
  `-window` or `-linker` (NCBI clamps them).
- The `blastdbcmd` and `blastdb_aliastool` wrappers cover what RepeatMasker and
  RepeatModeler call, not the full tools. Unsupported options exit non-zero, as
  do a missing `-entry` and a batch with entries skipped (as in NCBI). A batch
  interval that starts past the end of its sequence is skipped with an error;
  NCBI silently returns the whole sequence for it.
- `blastdbcmd` reads the `.2bit` our `makeblastdb` writes or a FASTA beside the
  database, never an NCBI BLAST index. For a database BuildDatabase built with
  the real NCBI `makeblastdb`, it therefore falls back to BuildDatabase's input
  FASTA, and reads `<db>.translation` so that it still reports the `gi|N` names
  the index holds. It also uppercases the bases, the way a nucleotide database
  does. RepeatModeler's Refiner needs the names; RAMExtend and the round
  searches see the case.

## Repository layout

| Path | What it is |
|---|---|
| `rmblast-lib/` | The engine: lookup tables, extension, traceback, DUST, statistics |
| `rmblastn/` | The `rmblastn` command-line front end |
| `dustmasker/` | The `dustmasker` command-line front end |
| `rmblast-stats/` | E-value and bit-score statistics, the vendored ALP sources, and the fitting tools that generated the baked tables |
| `build-support/` | Shared build-script code: puts the git commit in `--version` |
| `wrappers/`, `verify_wrappers/` | Caller-facing shims described above |

## Credits

The RepeatMasker extensions this port reimplements are the work of Arian Smit,
Robert Hubley, and Jeb Rosen at the Institute for Systems Biology, built on
NCBI's BLAST toolkit. ALP is by John Spouge at NCBI. Complexity-adjusted scoring
and masklevel filtering follow Phil Green's cross_match (www.phrap.org).

Upstream RMBlast: <https://www.repeatmasker.org/rmblast/>

## License

CC0 1.0 Universal. See `LICENSE`.

The NCBI and ALP sources this port derives from and vendors are "United States
Government Work" in the public domain, distributed under NCBI's PUBLIC DOMAIN
NOTICE, which imposes no terms CC0 could conflict with. That notice asks that
you cite the author in any work based on the material. CC0 requires no
attribution, but the request costs nothing to honour: if you publish work using
this software, please cite NCBI BLAST and RMBlast. `LICENSE` reproduces both
texts.
