# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## Unreleased
### Behaviour changes
Input that used to run is now refused with an error that names the problem.
Each case below either gave wrong or incomplete output, or hung or crashed
- Bad option values are usage errors: --word-size below 4, --gapopen below 0,
  --gapextend below 1, gap costs over 1000000, --num-threads 0, --mt-mode other
  than 0 or 1, an unknown --seed-mode, --evalue that is not a positive number,
  negative --xdrop-* or --min-raw-gapped-score, negative -matrix_lambda/_k/_alpha,
  and any --dust value except yes, no, true or false. NCBI also accepts negative
  gap costs and -mt_mode 2; we don't. NCBI's -dust 'level window linker' form is
  refused, where it used to run silently with the default parameters
- --seed-mode is removed. Its non-default "separate" mode had known
  minus-strand bugs; the default seeding, which matches NCBI, is unchanged
- --gapopen and --gapextend are required. Their defaults were 4/4; NCBI's are
  0/0, which NCBI itself rejects
- A database with two sequences of the same id is refused. Before, the first
  was never searched and the last was searched twice
- An empty database is refused
- A query or database sequence longer than 2^31 - 1 bases is refused instead of
  silently giving no hits
- Malformed matrix files are refused with the line number: a score that is not
  an integer, a short or long row, a matrix that is not square, a row label that
  does not match its header symbol, a repeated header symbol, a missing or
  malformed FREQS line, and a matrix whose expected score is not negative.
  NCBI also refuses the missing FREQS line, short rows and non-square matrices;
  it silently gives no hits for a positive expected score (RepeatModeler's
  at.matrix.4.2)
- Malformed FASTA (query, database, or dustmasker input) is refused with the
  line number: sequence before the first '>', a '>' with no id, ';' comment
  lines, a defline that is not UTF-8, whitespace inside a sequence line, and any
  character that is not a nucleotide code. Before, these were skipped, counted
  as N, or aborted with "stream did not contain valid UTF-8"
- dustmasker refuses an empty record (NCBI also aborts), and its usage errors
  exit 1 as NCBI's do
- -outfmt accepts only 0 and "6 <fields>" with the 22 implemented fields. Other
  formats were treated as 6, and unknown field names were dropped, shifting the
  columns after them. Now both are refused, as is a bare "6"
- The makeblastdb wrapper refuses options it does not implement, -dbtype prot
  and -mask_data. All three used to "succeed": prot wrote nucleotide dummy files
- The blastdbcmd wrapper exits 1 for an unknown -entry, for a batch with
  skipped entries (both as NCBI does), and for a batch interval that starts
  past the end of its sequence

### Fixed
- --word-size defaulted to 8; it now defaults to 11, as in NCBI. RepeatMasker
  and RepeatModeler always pass it
- --mt-mode 1 hung forever if the search of one query panicked; it now exits
  with an error naming the query
- Word size 0 to 2 panicked and word size 3 hung; negative gap costs, or gap
  costs near 2^30, hung
- A '-' in a query or database sequence was treated as a gap code: it printed
  as '-' and scored differently from NCBI. It is now N, as in NCBI
- A truncated or corrupt .2bit file panicked, or loaded with fewer sequences;
  it is now an error that says what is wrong and where
- Matrices with a Z column (ctools.matrix, symmetricNscores1matrix) scored N as
  +3 instead of the file's -1: Z took N's slot. Header symbols outside the IUPAC
  codes are now discarded as NCBI discards them, and all four such matrices give
  output identical to NCBI
- A FREQS line written as `#  FREQS ...` was missed, which left lambda at 0 and
  turned -complexity_adjust off silently. NCBI finds FREQS anywhere on a comment
  line, and so does the port now
- Leading spaces on a FASTA sequence line were counted as N; they are now
  ignored, as in NCBI
- dustmasker printed deflines verbatim; it now prints them as NCBI does (spaces
  collapsed, cut at a tab, trailing periods stripped, HTML entities decoded),
  byte-identical to NCBI on 13,000 random deflines
- When E-values come from the built-in table, stderr says the entry was chosen
  by the matrix file name
- Pairwise output (-outfmt 0) now matches NCBI 2.17.1 line for line, apart from
  "Posted date" and the database title (see README). It now prints the
  2.17.1 version line, thousands separators in database sizes, the full
  subject defline wrapped as NCBI wraps it, DUST-masked query bases in lower
  case, a "No hits found" block for queries without hits, and NCBI's blank
  lines between queries. Tabular output is unchanged
- Queries of 2 Mb or more (searched in 1 Mb chunks) used an older minus-strand formula
  for the start of the preliminary gapped alignment. At word size 7 it gave
  different hits from NCBI (33 of about 4,200 lines on a 3 Mb chr22 slice); it
  now uses the same formula as single-chunk queries and matches NCBI
- Complexity adjustment (-complexity_adjust) of a minus-strand hit counted the
  forward query's bases; NCBI counts the reverse complement. This matters only
  for matrices whose FREQS are not complement-symmetric, which the RepeatMasker
  matrices are
- With --min-raw-gapped-score 0 (the default), zero-score, zero-length hits
  could be printed, and the ungapped cutoff was 0 where NCBI's is at least 1
- -mask_level used integer arithmetic where NCBI divides in floating point
  (e.g. 29 of 50 bases is 57% in NCBI, 58% before), and at -mask_level 0 an
  HSP only touching a better one was not masked
- The rmblastn wrappers mangled negative option values: `-matrix_beta -4.9`
  became `--4.9` and failed. Values now pass through unchanged
- The makeblastdb wrapper needed faToTwoBit at a fixed path, so it failed on any
  other host. It now looks in $UCSCTOOLS_DIR, then on PATH, and without
  faToTwoBit it still builds a usable database. It writes every file under a
  temporary name and renames it into place
- The blastdbcmd wrapper's 2bit commands exited 0 when twoBitToFa failed;
  -ctrl_a, -target_only and -long_seqids swallowed the next argument; a
  "> name" defline lost its record; and -out was truncated before the arguments
  were checked
- The verify wrappers failed every search under RepeatMasker 4.2.6, which
  drives a 3.x rmblastn with GNU-style options and a matrix path; NCBI rejected
  both. They now give NCBI its own spelling. They also share a query on stdin
  between the two engines, clean up their temporary directory, run without GNU
  time, and take NCBI_DIR from the environment
- A gilist that matches no sequence and an empty query file used to exit with
  no output and no message. Both now warn on stderr (the empty query with
  NCBI's own warning) and give the normal empty output

### Performance
- Startup took 28 s on a database of 100,000 sequences; it takes 0.4 s
- Per-subject work for short queries against many database sequences: the
  query was DUST-masked and reverse-complemented again for every database
  sequence, a 4096-entry buffer and a hash table were allocated and zeroed for
  each, and each was given a seeded random generator whether or not it had an
  ambiguous base. A 60 bp query against 100,000 sequences took 14-42 s; it
  takes 0.3 s
- The gapped-alignment stage copied each database sequence, and its reverse
  complement, for every query instead of using the copies already cached: two
  chromosome-sized copies per thread per query when searching against a
  genome. 200 TE queries against chr20: peak memory 360 MB -> 294 MB
- -mask_level filtering (and --prelim-cull) was quadratic in the number of
  HSPs kept, and slowed further after one long HSP was kept: 1 million HSPs
  took 12.6 s at -mask_level 99; it takes 1.1 s
- --mt-mode 1 read the whole query file before searching; it now reads each
  query as it is searched, so memory holds only the queries in flight. A
  malformed query record is now reported after the output of the queries
  before it, as in --mt-mode 0

## 3.0.12 - 2026-11-04
### Fixed
- rmblastn could take many hours on a search producing millions of HSPs per
  query, such as a low-complexity library family with -dust no. The step that
  purges HSPs with common endpoints moved whole HSP records once per removal,
  which is quadratic; it is now one linear pass with the same results. A Xenopus
  laevis search that was killed after 5.5 hours with no output now finishes in
  about 2 hours 10 minutes on 20 threads
- Long queries could give different hits from NCBI 2.17.1 next to the 1 Mb
  query chunk boundaries: 15 of about 3 million hits at 7 sites on Xenopus
  laevis, while the existing comparison suites showed no differences. The merge
  of HSPs across chunks now matches NCBI in three ways:
  - HSPs are paired in NCBI's order, best score first, and one HSP can absorb
    several others
  - The merged score is the score density over the merged length, or the better
    of the two scores if that is higher
  - Only plus-strand HSPs merge; NCBI never merges minus-strand HSPs across chunks

  This search now gives output identical to NCBI 2.17.1 on a FASTA database.
  Results near chunk boundaries can differ from 3.0.11
- Prelim HSPs now use a stable sort before traceback, so equal HSPs keep their
  order

## 3.0.11 - 2026-10-21
### Fixed
- rmblastn -mt_mode 1 could deadlock (all threads idle, no output) when a
  query job waiting for the writer was run inside another query's parallel
  search; queries are now handed out by a thread outside the worker pool and
  query jobs never wait

## 3.0.10 - 2026-10-20
### Breaking Changes
- None
### Added
- None
### Fixed
- Results are now streamed rather than buffered in -mt_mode 1
- The entire project including dustmasker and rmblast-lib now share
  the same project version number
- Silenced harmless compiler warnings from the vendored ALP library
  (-Woverloaded-virtual) that were shown on every build
### Removed
- None
### Known Issues
- None

## 3.0.0 - 2020-09-02
- Initial release of standalone Rust port
