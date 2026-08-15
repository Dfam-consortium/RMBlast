# Design Notes — rmblastn Rust Port

## Subject Database Format

### Decision (2026-06): Support FASTA and 2bit as Database Inputs

The Rust port accepts either FASTA or UCSC `.2bit` files as the subject database
(`--db`).  Format is autodetected:

1. **Extension** (case-insensitive): `.fa`, `.fasta`, `.fas`, `.fna`, `.fn` → FASTA.
   `.2bit` → 2bit.
2. **Magic-byte fallback** for files with unrecognised or absent extensions: peek at
   the first byte.  `>` → FASTA; otherwise treated as 2bit.

#### FASTA database (recommended)

Full IUB/IUPAC support.  Ambiguous codes (R/Y/K/M/W/S/B/D/H/V/N) are preserved at
load time and resolved to specific bases at decode time using NCBI's constrained-random
algorithm (`CAmbigDataBuilder::x_Random`): seed = sequence length, one `GetRand()` call
per ambiguous position in 5′→3′ order, base = `NCBI_IUPAC_EXPAND[code][rand & 3]`.

This exactly matches what `makeblastdb` bakes into an NCBI BLAST database, so the Rust
binary reproduces NCBI rmblastn output bit-for-bit on the same FASTA input.

#### 2bit database

UCSC `.2bit` encodes only the four unambiguous bases (T/C/A/G) and a separate N-block
table.  All non-N IUB codes (`R/Y/K/M/W/S/B/D/H/V`) in the original FASTA are
converted to `N` by `faToTwoBit`, and that information is permanently lost.

When a 2bit database is opened, a note is printed to stderr:

> Note: .2bit database — non-N IUB ambiguity codes (R/Y/K/M/W/S/B/D/H/V) are not
> preserved in the .2bit format and will be treated as N.  For full IUB fidelity
> supply the database in FASTA format.

All ambiguous positions are randomized as N (code 14), meaning every ambiguous site
draws from `{A, C, G, T}` with equal probability rather than the constrained subset
implied by the original IUB code.  This produces different seeds than NCBI for
sequences that contain non-N IUB ambiguities, so 2bit-mode output will diverge from
NCBI on such databases (e.g. the `longlib` repeat library, which contains K/Y/S/R
codes in several families).

**Use FASTA for validation against NCBI.  Use 2bit only when space is the binding
constraint and the database contains only A/C/G/T/N.**

---

## Future Proposal: Extended 2bit Format with IUB Auxiliary Table

*For discussion with the UCSC Genome Browser team.*

The UCSC `.2bit` format is built around two tables per sequence:

- **Sequence block table** — 2-bit-packed bases (T=0, C=1, A=2, G=3).
- **N-block table** — run-length-encoded ranges of N positions.

### Proposed extension

Add an optional third table per sequence: an **IUB-block table** in the same
run-length style as the N-block table, but carrying the original IUPAC code for
each ambiguous run:

```
iup_count   uint32   number of IUB blocks
For each block:
  start     uint32   0-based start position in the decoded sequence
  size      uint32   block length in bases
  iub_code  uint8    BLASTNA code (4–13: R/Y/M/K/W/S/B/D/H/V)
  _pad      uint8[3] alignment padding
```

Backward compatibility: readers that do not know about the IUB table ignore it
(they already use the N-block table to identify ambiguous positions and randomize as N).
New readers that understand the table can apply constrained-random resolution per IUPAC
code, matching NCBI makeblastdb behavior.

### Benefits

- **8× space savings** of 2bit over FASTA are retained for the main sequence data.
- **Full IUB fidelity** for tools that need it (rmblastn, makeblastdb-compatible seeders).
- **Backward compatible** — old `.2bit` readers continue to work.
- **Random access** — large genomes with gilist-style subset access remain efficient;
  auxiliary tables are small relative to sequence data.
- **Single file** — no sidecars.

### Comparison of formats

| Format         | Size (vs FASTA) | IUB support | Random access |
|----------------|-----------------|-------------|---------------|
| FASTA          | 1×              | Full        | Sequential    |
| bgzf FASTA     | ~0.25–0.33×     | Full        | Block-level   |
| UCSC 2bit      | ~0.125×         | N only      | Full          |
| 2bit + IUB tbl | ~0.125× + ε     | Full        | Full          |

The IUB table adds at most a few KB per sequence (only ambiguous runs are stored),
so the overall file size is essentially identical to plain 2bit for typical repeat
libraries or genome assemblies.

### Status

Proposal only — not implemented.  Raise with UCSC Genome Browser team before
committing to a format extension, as backward compatibility and format versioning
need coordination.

---

## Why not bgzf-compressed FASTA?

bgzf (blocked gzip used by samtools/htslib) gives ~3–4× compression over plain FASTA
because it is compressing ASCII text (one byte/base).  UCSC 2bit achieves ~8× because
it packs four bases into one byte at the binary level.  For large databases where space
or I/O bandwidth dominates, 2bit is significantly more efficient.  bgzf does support
random access (via a `.gzi` index) but does not improve per-base density.

---

## Minus-strand representation (ungapped uses FWD_q/RC_s — PROVEN FAITHFUL 2026-06-24)

This records the two different internal frames the port uses for minus-strand
processing, and the 2026-06-24 finding that the ungapped `FWD_q/RC_s` frame — long
suspected as the source of the minus-strand bug class — is in fact **bit-for-bit
equivalent** to NCBI's `RC_q/FWD_s` frame and is NOT the cause of bug #36.

### TL;DR (2026-06-24)
- **Ungapped extension** runs in `FWD_q/RC_s` with complement+left/right-swap.
  Unit test `search::ungapped::tests::test_minus_frame_matches_ncbi_frame_asym`
  proves it produces byte-identical scores AND endpoints to NCBI's `RC_q/FWD_s`
  frame, even on an asymmetric, non-complement-symmetric matrix (the 20p39g case).
  It is faithful; it was **kept**, not refactored.
- **Gapped Phase 2a/2b** already run in NCBI's `RC_q/FWD_s` frame (`query_rc` +
  forward subject), with a DP tie-break byte-identical to `blast_gapalign.c`.
- Bug #36's real root is therefore NOT the ungapped frame. It is downstream:
  containment / `mask_level` filtering (drops HAL1ME/L1PA12 hits) and minus-strand
  traceback gap-justification (the one-base gap-shift ties at engine.rs:1596).

The historical record below is kept because the `FWD_q/RC_s` choice is still the
reason earlier bugs (#30, #32) needed minus-specific patches, and because the
performance note at the end is still relevant if the frame is ever revisited.

### NCBI's representation (the faithful one we now mirror)

For a minus-strand hit NCBI reverse-complements the **query** into a separate "minus
context" (`RC_q`), keeps the **subject forward** (`FWD_s`), and runs **one uniform code
path** for both strands: same seeding, same ungapped extension (`matrix[query][subject]`),
same left-then-right Xdrop with the same fixed/adaptive assignment, same diagonal dedup,
same traceback tie-breaks.  Minus strand is literally "plus strand, computed on the
reverse-complement of the query."  Only the query is duplicated (both strands live in the
lookup table); the subject is scanned forward exactly once.

### The port's original optimisation (FWD_q / RC_s)

To avoid building/storing a second (RC) query buffer and to reuse the **forward-query**
lookup table for both strands, the port instead kept the **query forward** and
reverse-complemented the **subject** (`FWD_q, RC_s`) during the scan and ungapped
extension.  This is the mirror image of NCBI's frame.  Because the substitution math is
complement-symmetric the *total* alignment score is identical, so the shortcut "worked"
for most cases.

To make the mirror frame produce NCBI-equivalent results it had to be patched in three
places:
  1. **base-complement before the matrix lookup** in ungapped extension
     (`matrix.score(COMPLEMENT[qb], COMPLEMENT[sb])`) — bug #30's fix;
  2. a **swapped left/right extension order** with swapped fixed/adaptive Xdrop
     assignment in the minus path of `extend_ungapped`;
  3. a different matrix index order between the (FWD_q/RC_s) ungapped stage and the
     gapped stage — note the **gapped stage already used NCBI's `RC_q/FWD_s` frame**
     (it builds `query_rc`), so the two stages disagreed on frame.

### Why it was abandoned

The reframing is **score-equivalent but NOT decision-equivalent**.  The intermediate
steps that decide *which* alignment survives are order- and direction-sensitive:
Xdrop termination (partial-sum profile depends on walk direction), diagonal dedup
(depends on which seed extends first and how far), and traceback tie-breaks.  Mirroring
the frame silently changes these decisions.  Every minus-strand bug we chased traces to
this single cause:
  - **#30** ungapped ambiguity scored in the wrong frame;
  - **#32** `purge_common_start` tiebreak (gap-in-query vs gap-in-subject) inverted;
  - **#36** minus ungapped extension over-extended (FWD-s 1767 vs NCBI 1396), which
    dedup-suppressed a seed NCBI keeps, anchoring the gapped alignment at the wrong place.

**Correction (2026-06-24):** #30 and #32 were genuine, but they were fixed in place
(the complement patch and the tiebreak ordering). The remaining suspicion that the
`FWD_q/RC_s` ungapped frame still silently diverged from NCBI was **tested and refuted**
— see the TL;DR above. No frame refactor was performed; the equivalence test now guards
the behaviour. Bug #36 is being pursued in the containment/traceback code instead.

**Resolution (2026-06-26):** the residual minus-strand divergences (#36 over/under-
extension, #38 seed-anchor gap-ties, and the new `cc` false-positive) were ultimately
**all one bug — the preliminary bidirectional *split convention* (#42), see below** —
not the ungapped frame, not seed-anchoring, and not ambiguity randomization (those were
all dead ends, instrumented and ruled out).  After the #42 fix the entire mirs + cc
residual class is closed and all four validation suites are 16/16 with zero diffs across
all verify bundles.

### If revisiting for performance

The FWD_q/RC_s approach's only advantage was avoiding a second query buffer and reusing
the forward-query lookup table.  Any future attempt to reclaim that must reproduce NCBI's
**decisions**, not just its scores — i.e. it must replay Xdrop termination, diagonal
dedup, and tie-breaks identically in the mirrored frame.  That is the hard part and the
reason the optimisation was not worth its correctness cost.  A cheaper win is to share the
single RC query buffer across all subjects of a search (build once) rather than to avoid
building it at all.  **(Update 2026-06-26: this "share the single RC query buffer"
idea was partly realised on the subject side — see "Subject decode cache" below.)**

---

## Gapped alignment: split convention, orientation, traceback, purge (2026-06-26)

This section documents the core gapped-alignment design after the 2026-06-25/26 work,
which corrected several long-standing minus-strand and asymmetric-matrix divergences.

### Bidirectional split convention — the pivot is NCBI's `q_length`/`s_length` (#42)

The preliminary and traceback gapped extensions are **bidirectional**: from a pivot
they extend LEFT (reverse) and RIGHT (forward) and sum the two scores.  The subtle
design point is **which side owns the pivot base**.

NCBI's `s_BlastDynProgNtGappedAlignment` rounds the seed up to the next 4-base subject
boundary (`offset_adjustment = 4 - s_off % 4`) to get `q_length`/`s_length`, then calls
`s_BlastAlignPackedNucl(query, subject, q_length, s_length, …)` for the LEFT extension
and `(query+q_length-1, subject+…, …)` for the RIGHT.  So **`q_length`/`s_length` is the
EXCLUSIVE left bound and the FIRST base of the right extension**: left covers `[0, q_length)`,
right covers `[q_length, end)`.

The port's prelim (`gapped_extend_score_only`, gapped.rs) must mirror this exactly:

```
left  = align_ex_score_only(&qa[..q_seed],  &sa[..s_seed],  REVERSE)   // pivot EXCLUDED
right = align_ex_score_only(&qa[q_seed..],  &sa[s_seed..],  FORWARD)   // pivot = first cell
```

A one-base error here (including the pivot in the LEFT extension, `..q_seed+1`) leaves
the span and endpoints unchanged but scores the pivot base in a different running/xdrop
context, so under a tight prelim x-drop the band is pruned differently and the prelim
*filter* score drifts by a few points.  That is invisible except for hits sitting right
at `min_raw_gapped_score`, where it flips keep↔drop.  **This single off-by-one was the
root of the entire mirs (plus, missing-hit) and `cc` (minus, false-positive + gap-tie)
residual class** — superseding the earlier #36/#38/randomization theories.

The pivot must be `first-of-right` for **both** strands.  The plus branch already passes
`q_length` directly.  The minus branches compute it in NCBI's `(RC_q, FWD_s)` frame and
map back: `run_gapped_phase` uses `pre_q = qa_len-1-q_len_piv` (so the caller's
`qa_len-1-pre_q` yields the pivot `= q_len_piv`), and `run_phase2a_inner`'s RC-subtract
formula already yields `first-of-right`.  All three paths are now consistent.

### Preliminary score orientation is subject-outer (#41 — NOT a bug)

The prelim scores `matrix[subject][query]` (subject indexes the matrix row).  This is
correct and matches NCBI: `s_BlastAlignPackedNucl(Uint1* B, Uint1* A, …)` is called as
`(query, subject, …)`, so **B = query, A = subject**, and the outer loop runs over
`M = s_length` with `matrix_row = matrix[A[a_index]] = matrix[subject]`.  The
final/reported score instead comes from the query-outer **traceback**
(`align_ex(query, subject)` = `matrix[query][subject]`), which is also faithful.  A
2026-06-25 hypothesis that the prelim should be query-outer was a misread of that arg
order; transposing it is catastrophic (forces the wrong orientation) and must not be
retried.  The prelim is only a *filter*; the reported score is the traceback's.

### Traceback DP structure — single-pass `ALIGN_EX` (#2)

`align_ex_inner` (gapped.rs) mirrors NCBI's `ALIGN_EX`:

- **No score-only pre-pass.**  The score array starts small (`num_extra+101`), is grown
  by **doubling realloc** when the band's right edge nears the allocation
  (`last_b_index + num_extra + 3 >= dp_cap`, capped at `n+2`), via a raw `dp_ptr`
  refreshed after each `Vec::resize` (resizing `ws.dp` does not touch `ws.flat_edit`).
- **Traceback grows incrementally**, not as the old O(dp_cap²) pre-reserve.  Each row
  reserves its worst case (`(b_size-first_b_index) + num_extra + slack`) and writes the
  per-cell script bytes through a raw `flat_ptr` (refreshed per row after the reserve);
  `row_info[a] = (flat offset of row a, first_b_index)`.  Total traceback memory is
  O(total band cells), matching NCBI's `s_GapGetState` chunked scheme.  The raw-pointer
  write is the literal analogue of NCBI's `edit_script_row[b_index] = script` and avoids
  a per-cell capacity check on this hot path.

### Common-endpoint purge — faithful two-pass cut/delete (#3)

`purge_hsps_with_common_endpoints(arr: &mut [Option<AlignResult>], purge: bool)` is a
direct port of `Blast_HSPListPurgeHSPsWithCommonEndpoints` (blast_hits.c): `None` models
a freed/NULL slot; two sub-passes (common START via `cmp_query_offset`, common END via
`cmp_query_end`) each MOVE removed HSPs (cut → `Some`, freed → `None`) to the end of the
shrinking array.  The traceback callers run it twice exactly as NCBI does: `purge(false)`
(cut) → re-trace the cut remainders `arr[extra_start..]` in array order
(`reevaluate_gapped`) → flatten (PurgeNull) → `purge(true)` (delete).  A stable `sort_by`
supplies glibc-qsort's quasi-stable ordering for the small per-subject arrays.

---

## Subject decode cache (2026-06-26)

`SubjectDb` (seq/mod.rs) wraps a private `SubjectBackend` (FASTA or 2bit) plus a
**decode-once cache**: at `open`, every subject is decoded to BLASTNA (ambiguities
resolved by `NcbiRandom`) and its n_mask computed once, stored as `Arc<[u8]>`.
`get_full_sequence_blastna`/`get_n_mask` hand out cheap `Arc` clones.

**Why:** previously each call re-`clone()`d the raw sequence and re-ran the ambiguity
decode into a fresh `Vec`.  For the many-query × large-subject shape (TE-library query ×
genome DB), the genome was re-decoded *per query* (hundreds of times under `--mt-mode 0`;
one copy per worker under `--mt-mode 1`).  The decode is deterministic (`NcbiRandom`
seeded only by `dna_size`), so caching is bit-for-bit identical to re-decoding.  Backend
per-call methods are unchanged (used to populate the cache and by their own unit tests).
Callers are unchanged thanks to deref coercion (`&seq`/`&n_mask` coerce from `Arc<[u8]>`
exactly as from `Vec`).

---

## Threading model (`--mt-mode`, mirrors NCBI)

- **`--mt-mode 0` (default, SplitByDb):** queries processed one at a time; within each,
  `search_db_parallel` does `subjects.par_iter()`.  Scales with the number of DB
  sequences — the right model for the standard RepeatMasker shape (genome query ×
  many-sequence TE library), where it scales near-linearly.
- **`--mt-mode 1` (SplitByQueries):** `all_queries.par_iter()`, one query per thread
  (with the subject `par_iter` nesting underneath rayon's shared pool).  Scales with the
  number of query sequences — the right model for the swap shape (many TE queries ×
  single-sequence genome DB), where mode 0 cannot parallelize.

Output is bit-identical across both modes and thread counts (the final per-query score
sort makes ordering deterministic).

---

## Memory architecture & known divergences from NCBI

These are deliberate trade-offs, fine for the typical RepeatMasker workload but worth
revisiting for large-DB use:

- **Whole DB resident, no volume streaming.**  NCBI streams DB volumes; the port loads
  and decodes the whole DB into memory.  Negligible for the typical small TE library
  (≤~1 MB; ~0.1% of footprint, loads in ms, search is compute-bound).  It only matters
  for large-subject cases (genome-as-subject / nt-sized DB) — there NCBI's streaming wins
  and the port's footprint grows with DB size.  Revisit only if large-DB use cases arise.
- **Query copies.**  The port holds the unmasked query, a masked copy (lookup), and —
  for large (multi-chunk) queries — a reverse-complement (`query_rc`).  NCBI stores the
  query as both strands × (masked `sequence` + `sequence_nomask`) ≈ 4× query length; the
  port is at or below that (single-chunk keeps essentially just the unmasked query +
  lookup; the masked copy is transient and minus uses the small per-subject subject-RC).
  The one genuinely port-only allocation is the transient `[fwd | sep | rc]` buffer that
  `build_query_lookup` concatenates to feed the lookup constructor; NCBI indexes its
  already-both-strands query in place.  Only matters for very large single queries.
- **Swap-orientation performance.**  In the standard orientation the port is ~10–20%
  faster than NCBI single-threaded (~30–40% at 4 threads).  The swap orientation
  (TE query × genome DB) was historically ~1.7× *slower*; callgrind showed the true
  cost was not the scan itself (the scan already reads 2-bit packed data) but
  **re-deriving the subject's RC + packed strands for every query** (~28% of
  instructions: `revcomp_blastna` + `blast_compress_blastna_sequence` per
  query×subject task).  Fixed by the `PreparedSubject` cache on `SubjectDb`
  (`get_prepared`): RC + both packed strands are computed once per subject and
  shared (deterministic, bit-identical).  Measured longlib × chr20 single-thread:
  3:17.6 → 1:59.0 wall (−40%); the port now beats NCBI in all 16 swap combos.
  Cost: ~1.5× subject size cached per touched subject (chr20: +96 MB RSS).  The
  same cache feeds `search_phase2a`, which previously re-derived each subject per
  query chunk (~50× on chr22) in the standard orientation.
