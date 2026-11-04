# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

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
