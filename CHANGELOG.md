# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

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
