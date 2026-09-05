//! Full search pipeline: lookup → ungapped filter → interval tree → gapped → output.
//!
//! Mirrors the NCBI rmblastn two-phase hot path:
//!   Phase 1 (BlastNaWordFinder in na_ungapped.c):
//!     NaLookupTable scan → exact extension → diagonal tracker → ungapped extension
//!     → collect all ungapped hits sorted by score.
//!
//!   Phase 2a (BLAST_GetGappedScore in blast_gapalign.c) — PRELIMINARY:
//!     For each ungapped hit (score descending):
//!       interval-tree containment check → preliminary gapped alignment (xdrop_gap=30)
//!       → score filter → add to interval tree.
//!
//!   Phase 2b (Blast_TracebackFromHSPList in blast_traceback.c) — TRACEBACK:
//!     For each accepted preliminary HSP:
//!       BlastGetStartForGappedAlignmentNucl (find best seed in preliminary region) →
//!       full gapped alignment with traceback (xdrop_gap_final=100) →
//!       score/complexity filter → output.
//!
//! Threading: parallelism is over subject sequences (see search_db_parallel in main.rs).

use crate::hits::score_compare_hsps;
use crate::search::diag_hash::BlastDiagHash;
use crate::search::itree::{BlastIntervalTree, ITreeHsp};
#[cfg(feature = "diagnostics")]
use std::sync::atomic::AtomicU64;

#[cfg(feature = "diagnostics")]
pub static COUNT_SEEDS: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "diagnostics")]
pub static COUNT_UNGAPPED_HITS: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "diagnostics")]
pub static COUNT_PRELIM_GAPPED: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "diagnostics")]
pub static COUNT_FINAL_GAPPED: AtomicU64 = AtomicU64::new(0);
#[cfg(feature = "diagnostics")]
pub static COUNT_FINAL_HITS: AtomicU64 = AtomicU64::new(0);

/// Diagnostic counter: times `improve_seed`'s minus branch produced a negative
/// `max_offset` and fell back to the unimproved seed.  Expected to stay at 0 — it is
/// only reachable at `score == 1` with the scan pinned to the array start.
#[cfg(feature = "diagnostics")]
pub static IMPROVE_SEED_NEGATIVE_OFFSET: AtomicU64 = AtomicU64::new(0);

use crate::blast::mblookup::{blast_mb_lookup_table_new, blast_mb_scan_subject, BlastMBLookupTable};
use crate::blast::nalookup::{blast_na_lookup_table_new, choose_lut_width, BlastNaLookupTable};
use crate::blast::nascan::blast_na_scan_subject;
use crate::blast::types::BlastOffsetPair;
#[cfg(test)]
use crate::blast::util::{blast_compress_blastna_sequence, seqblk_from_blastna};
use crate::seq::{prepare_subject_strands, PreparedSubject};
use crate::encoding::{revcomp_blastna, BLASTNA_COMPLEMENT};
use crate::filter::dust::{dust_mask, dust_mask_ncbi_compat, DUST_LEVEL, DUST_LINKER, DUST_WINDOW};
use crate::hits::{EditOp, EditScript, Hsp, Strand};
use crate::matrix::ScoreMatrix;
use crate::options::{SearchParams, SeedMode};
use crate::output::AlignResult;
use crate::search::complexity_adjust::apply_complexity_adjust;
use crate::search::gapped::{align_ex, extract_aligned, gapped_extend_score_only, AlignWorkspace, GapAlignResult};
use crate::search::ungapped::{extend_ungapped, UngappedResult};
use crate::stats::{blastna_to_iupac_aligned, compute_align_stats};

/// A post-deduplication seed in FWD-normalized coordinates.
///
/// `q_off` is chunk-relative (position within the current query chunk, not the full genome).
/// `s_off` is the subject left-end position on the forward strand (0-based, before ungapped extension).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeedRecord {
    pub q_off: u32,
    pub s_off: u32,
    pub strand: Strand,
}

/// Zero-cost abstraction for optional seed collection.
/// Monomorphized away (no-op path) in production code.
trait SeedOut {
    fn push_seed(&mut self, seed: SeedRecord);
}
struct DiscardSeeds;
impl SeedOut for DiscardSeeds {
    #[inline(always)] fn push_seed(&mut self, _: SeedRecord) {}
}
impl SeedOut for Vec<SeedRecord> {
    #[inline(always)] fn push_seed(&mut self, s: SeedRecord) { self.push(s); }
}

/// A passing ungapped alignment in FWD-normalized coordinates.
///
/// All coordinates are 0-based and relative to the query chunk (not the full genome).
/// `s_start`/`s_end` are on the forward subject strand regardless of strand.
/// Hits are collected in score-descending order (same order as they enter gapped alignment).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UngappedHit {
    pub q_start: u32,
    pub q_end:   u32,
    pub s_start: u32,
    pub s_end:   u32,
    pub score:   i32,
    pub strand:  Strand,
}

/// Zero-cost abstraction for optional ungapped-hit collection.
trait UngapOut {
    fn push_ungap(&mut self, hit: UngappedHit);
}
struct DiscardUngap;
impl UngapOut for DiscardUngap {
    #[inline(always)] fn push_ungap(&mut self, _: UngappedHit) {}
}
impl UngapOut for Vec<UngappedHit> {
    #[inline(always)] fn push_ungap(&mut self, h: UngappedHit) { self.push(h); }
}

/// Dispatch wrapper for BlastNaLookupTable (lut_word_length ≤ 8) or
/// BlastMBLookupTable (lut_word_length > 8, eMBLookupTable in NCBI).
enum NaLookup {
    Small(BlastNaLookupTable),
    Mega(BlastMBLookupTable),
}

impl NaLookup {
    #[inline] fn lut_word_length(&self) -> i32 {
        match self { NaLookup::Small(l) => l.lut_word_length, NaLookup::Mega(l) => l.lut_word_length }
    }
    #[inline] fn scan_step(&self) -> i32 {
        match self { NaLookup::Small(l) => l.scan_step, NaLookup::Mega(l) => l.scan_step }
    }
    #[inline] fn longest_chain(&self) -> i32 {
        match self { NaLookup::Small(l) => l.longest_chain, NaLookup::Mega(l) => l.longest_chain }
    }
    fn scan(
        &self, subject: &[u8], pairs: &mut [BlastOffsetPair],
        max_hits: i32, range: &mut [i32; 2],
    ) -> i32 {
        match self {
            NaLookup::Small(l) => blast_na_scan_subject(l, subject, pairs, max_hits, range),
            NaLookup::Mega(l)  => blast_mb_scan_subject(l, subject, pairs, max_hits, range),
        }
    }
}

/// Record of a passing ungapped extension, ready for the gapped phase.
struct UngappedRecord {
    /// Seed position in query (start of word_size match, after exact extension).
    q_seed: u32,
    /// Seed position in subject (start of word_size match, after exact extension).
    s_seed: u32,
    /// Ungapped alignment query start (0-indexed in forward-query space).
    q_start: u32,
    /// Ungapped alignment query end.
    q_end: u32,
    /// Ungapped alignment subject start (0-indexed, in the strand's subject space).
    s_start: u32,
    /// Ungapped alignment subject end.
    s_end: u32,
    /// Ungapped alignment score.
    score: i32,
    strand: Strand,
}

/// Preliminary gapped HSP — carries the seed and preliminary boundaries for the traceback phase.
#[derive(Clone)]
pub struct PrelimHsp {
    /// Seed used for the preliminary alignment (= ungapped hit seed, in strand subject space).
    pub q_seed: u32,
    pub s_seed: u32,
    /// Preliminary gapped alignment boundaries (in strand's subject space).
    pub q_start: u32,
    pub q_end: u32,
    pub s_start: u32,
    pub s_end: u32,
    pub score: i32,
    pub strand: Strand,
}

/// min_diag_separation value for rmblastn megablast (MB_HSP_CLOSE threshold).
/// Mirrors CBlastNucleotideOptionsHandle::SetHitSavingOptionsDefaults in NCBI.
const MIN_DIAG_SEP: i64 = 50;

/// Parameters for the opt-in prelim-cull fast mode (`--prelim-cull`).
///
/// NOT part of the faithful NCBI pipeline.  When enabled, preliminary HSPs
/// that are mask_level-style dominated by higher-scoring prelims are skipped
/// in Phase 2b (round 1), then re-tested against the survivors' FINAL
/// coordinates and scores and resurrected — with a second small Phase 2b
/// wave — where the dominance no longer holds (round 2).  The exact
/// output-level mask_level still runs afterwards.  Output differs from
/// faithful mode only in a small fraction of heavily-overlapped, mostly
/// low-scoring hits (measured on chr22 x longlib: balanced preset ~0.2% of
/// annotated bp, ~0.5% of hits, for ~68% of traceback DP avoided).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PrelimCullParams {
    /// Round-1 dominator score margin in percent (110 → the dominating prelim
    /// must score >= 1.10x the candidate).
    pub cull_margin_pct: u32,
    /// Round-1 coverage requirement in percent of the candidate's query span.
    pub cull_coverage: u32,
    /// Round-2 resurrection margin in percent (dominator FINAL score vs
    /// candidate PRELIM score; < 100 is conservative).
    pub resurrect_margin_pct: u32,
    /// Round-2 slack in bp added to each end of the candidate span before the
    /// dominance re-test (absorbs prelim-vs-final coordinate drift).
    pub resurrect_slack: u32,
}

/// Slop (bp) for deciding that a prelim endpoint abuts a query-chunk edge.
const CULL_BOUNDARY_EPS: u32 = 5;

/// True if `x` lies within [`CULL_BOUNDARY_EPS`] of any coordinate in the
/// sorted `boundaries` list.
fn abuts_boundary(x: u32, boundaries: &[u32]) -> bool {
    let lo = x.saturating_sub(CULL_BOUNDARY_EPS);
    let i = boundaries.partition_point(|&b| b < lo);
    i < boundaries.len() && boundaries[i] <= x.saturating_add(CULL_BOUNDARY_EPS)
}

/// Round 1 of the prelim cull: mark prelims whose query span is
/// `cull_coverage`%-covered by a single not-yet-culled prelim scoring at
/// least `cull_margin_pct`% of the candidate.  Mirrors the
/// [`apply_mask_level`] sweep (score-descending order, survivors-only
/// dominators, minus-strand +1 shift).  Returns per-subject keep masks.
///
/// `chunk_boundaries` holds the sorted interior query-chunk edge coordinates:
/// a prelim whose span abuts one was truncated by Phase 2a's chunking, so its
/// extent under-represents its final traceback — such prelims are never
/// culled (they still act as dominators).
pub fn cull_prelims(
    prelims_by_subject: &[Vec<PrelimHsp>],
    cull_coverage: u32,
    cull_margin_pct: u32,
    chunk_boundaries: &[u32],
) -> Vec<Vec<bool>> {
    struct Cand {
        subj: usize,
        idx: usize,
        qs: u32,
        qe: u32,
        score: i64,
    }
    let mut cands: Vec<Cand> = Vec::new();
    for (subj, v) in prelims_by_subject.iter().enumerate() {
        for (idx, p) in v.iter().enumerate() {
            let (qs, qe) = if p.strand == Strand::Minus {
                (p.q_start + 1, p.q_end + 1)
            } else {
                (p.q_start, p.q_end)
            };
            cands.push(Cand { subj, idx, qs, qe, score: p.score as i64 });
        }
    }
    // Deterministic order: score DESC, then query start ASC, then identity.
    cands.sort_unstable_by(|a, b| {
        b.score.cmp(&a.score)
            .then_with(|| a.qs.cmp(&b.qs))
            .then_with(|| a.subj.cmp(&b.subj))
            .then_with(|| a.idx.cmp(&b.idx))
    });

    let mut keep: Vec<Vec<bool>> =
        prelims_by_subject.iter().map(|v| vec![true; v.len()]).collect();
    let mut accepted: Vec<(u32, u32, i64)> = Vec::new(); // qs ASC
    let mut max_span: u32 = 0;
    for c in &cands {
        let span = c.qe.saturating_sub(c.qs) as i64;
        if span == 0 {
            continue; // keep — the final mask_level handles degenerate spans
        }
        if abuts_boundary(c.qs, chunk_boundaries) || abuts_boundary(c.qe, chunk_boundaries) {
            // Chunk-truncated span: exempt from culling, but let it dominate.
            let pos = accepted.partition_point(|a| a.0 < c.qs);
            accepted.insert(pos, (c.qs, c.qe, c.score));
            max_span = max_span.max(c.qe - c.qs);
            continue;
        }
        let upper = accepted.partition_point(|a| a.0 < c.qe);
        let min_start = c.qs.saturating_sub(max_span);
        let mut masked = false;
        let mut i = upper;
        while i > 0 {
            i -= 1;
            let (qs_j, qe_j, sc_j) = accepted[i];
            if qs_j < min_start { break; }
            if qe_j <= c.qs { continue; }
            if sc_j * 100 < c.score * cull_margin_pct as i64 { continue; }
            let ovlp = (qe_j.min(c.qe) as i64) - (qs_j.max(c.qs) as i64);
            if ovlp * 100 / span >= cull_coverage as i64 {
                masked = true;
                break;
            }
        }
        if masked {
            keep[c.subj][c.idx] = false;
        } else {
            let pos = accepted.partition_point(|a| a.0 < c.qs);
            accepted.insert(pos, (c.qs, c.qe, c.score));
            max_span = max_span.max(c.qe - c.qs);
        }
    }
    keep
}

/// Round 2 of the prelim cull: re-test culled prelims against the surviving
/// FINAL hits and return the indices (into `culled`) that must be
/// resurrected because no single survivor still dominates them.
///
/// `survivors` holds FWD-query `(qs, qe, score)` of the wave-1 results
/// (minus-strand +1-shifted), sorted by `qs` ascending.
pub fn resurrect_prelims(
    culled: &[(usize, PrelimHsp)],
    survivors: &[(u32, u32, i64)],
    surv_max_span: u32,
    mask_level: u32,
    resurrect_margin_pct: u32,
    resurrect_slack: u32,
) -> Vec<usize> {
    let mut out = Vec::new();
    for (ci, (_, p)) in culled.iter().enumerate() {
        let (qs0, qe0) = if p.strand == Strand::Minus {
            (p.q_start + 1, p.q_end + 1)
        } else {
            (p.q_start, p.q_end)
        };
        let qs = qs0.saturating_sub(resurrect_slack);
        let qe = qe0 + resurrect_slack;
        let span = (qe - qs) as i64;
        let score = p.score as i64;
        let upper = survivors.partition_point(|a| a.0 < qe);
        let min_start = qs.saturating_sub(surv_max_span);
        let mut dominated = false;
        let mut i = upper;
        while i > 0 {
            i -= 1;
            let (qs_j, qe_j, sc_j) = survivors[i];
            if qs_j < min_start { break; }
            if qe_j <= qs { continue; }
            if sc_j * 100 < score * resurrect_margin_pct as i64 { continue; }
            let ovlp = (qe_j.min(qe) as i64) - (qs_j.max(qs) as i64);
            if ovlp * 100 / span >= mask_level as i64 {
                dominated = true;
                break;
            }
        }
        if !dominated {
            out.push(ci);
        }
    }
    out
}

/// Find an improved gapped alignment seed within the preliminary alignment region.
/// Mirrors NCBI's BlastGetStartForGappedAlignmentNucl (blast_gapalign.c).
///
/// Searches the preliminary alignment region for the longest run of consecutive
/// identical base pairs and uses the middle as the new anchor.
/// Returns (new_q_seed, new_s_seed).
fn improve_seed(
    query: &[u8],
    subject: &[u8],
    q_seed: u32,
    s_seed: u32,
    prelim_q_start: u32,
    prelim_q_end: u32,
    prelim_s_start: u32,
    prelim_s_end: u32,
    is_minus: bool,
) -> (u32, u32) {
    const HSP_MAX_IDENT_RUN: i32 = 10;

    let q_seed = q_seed as usize;
    let s_seed = s_seed as usize;
    let prelim_q_start = prelim_q_start as usize;
    let prelim_q_end = prelim_q_end as usize;
    let prelim_s_start = prelim_s_start as usize;
    let prelim_s_end = prelim_s_end as usize;

    let hsp_max_ident_run = ((HSP_MAX_IDENT_RUN as f64) * 1.5) as i32;

    if !is_minus {
        // Plus strand: check if the current seed is already in a good 10-base identity run.
        let mut score = -1i32;
        let (mut qf, mut sf) = (q_seed, s_seed);
        while qf < prelim_q_end && qf < query.len() && sf < subject.len() && query[qf] == subject[sf] {
            score += 1;
            if score > HSP_MAX_IDENT_RUN { return (q_seed as u32, s_seed as u32); }
            qf += 1; sf += 1;
        }
        let (mut qb, mut sb) = (q_seed, s_seed);
        loop {
            if query[qb] != subject[sb] { break; }
            score += 1;
            if score > HSP_MAX_IDENT_RUN { return (q_seed as u32, s_seed as u32); }
            if qb == 0 || sb == 0 { break; }
            qb -= 1; sb -= 1;
        }
        // Plus strand (context=0): search increasing from BELOW seed.
        // Mirrors BlastGetStartForGappedAlignmentNucl context=0.
        let offset = q_seed.saturating_sub(prelim_q_start)
            .min(s_seed.saturating_sub(prelim_s_start));
        let search_q_start = q_seed - offset;
        let search_s_start = s_seed - offset;
        let q_len = (prelim_s_end.saturating_sub(search_s_start))
            .min(prelim_q_end.saturating_sub(search_q_start));

        let mut max_score = 0i32;
        let mut max_offset = search_q_start;
        let mut score = 0i32;
        let mut is_match = false;
        let mut prev_match = false;

        for (i, qi) in (search_q_start..search_q_start + q_len).enumerate() {
            let si = search_s_start + i;
            if qi >= query.len() || si >= subject.len() { break; }
            is_match = query[qi] == subject[si];
            if is_match != prev_match {
                prev_match = is_match;
                if is_match {
                    score = 1;
                } else if score > max_score {
                    max_score = score;
                    max_offset = qi - (score / 2) as usize;
                }
            } else if is_match {
                score += 1;
                if score > hsp_max_ident_run {
                    let new_q_seed = qi - (hsp_max_ident_run / 2) as usize;
                    let new_s_seed = search_s_start + (new_q_seed - search_q_start);
                    return (new_q_seed as u32, new_s_seed as u32);
                }
            }
        }
        if is_match && score > max_score {
            max_score = score;
            let qi = search_q_start + q_len;
            max_offset = qi - (score / 2) as usize;
        }
        if max_score > 0 {
            let new_s_seed = search_s_start + (max_offset - search_q_start);
            (max_offset as u32, new_s_seed as u32)
        } else {
            (q_seed as u32, s_seed as u32)
        }
    } else {
        // Minus strand (context=1) early-return check: NCBI's original code runs the same
        // check for both strands. In NCBI context=1 coords, "forward" = increasing RC-gen
        // (= decreasing FWD-gen in Rust) and "backward" = decreasing RC-gen (= increasing
        // FWD-gen in Rust). NCBI forward is bounded by hsp->query.end (= prelim RC-gen end
        // = prelim_q_start in FWD-gen); NCBI backward is bounded only by array start (= no
        // FWD-gen upper bound). So: Rust backward scan bounded by prelim_q_start; Rust forward
        // scan unbounded (only array bounds).
        let mut score = -1i32;
        let (mut qf, mut sf) = (q_seed, s_seed);
        while qf < query.len() && sf < subject.len() && query[qf] == subject[sf] {
            score += 1;
            if score > HSP_MAX_IDENT_RUN {
                return (q_seed as u32, s_seed as u32);
            }
            qf += 1; sf += 1;
        }
        let (mut qb, mut sb) = (q_seed, s_seed);
        loop {
            if query[qb] != subject[sb] { break; }
            score += 1;
            if score > HSP_MAX_IDENT_RUN {
                return (q_seed as u32, s_seed as u32);
            }
            if qb <= prelim_q_start || sb <= prelim_s_start { break; }
            qb -= 1; sb -= 1;
        }
        // Early return didn't fire: scan the prelim HSP region for the longest identity run.
        // Scan decreasing in FWD-genomic from prelim_q_end - 1, mirroring NCBI's context=1 direction.
        let offset = prelim_q_end.saturating_sub(1).saturating_sub(q_seed)
            .min(prelim_s_end.saturating_sub(1).saturating_sub(s_seed));
        let search_q_top = q_seed + offset;  // ≈ prelim_q_end - 1 (high FWD-genomic boundary)
        let search_s_top = s_seed + offset;  // corresponding rc-s boundary
        let q_len = (search_q_top.saturating_sub(prelim_q_start) + 1)
            .min(search_s_top.saturating_sub(prelim_s_start) + 1);

        let mut max_score = 0i32;
        // Signed, because the post-loop index below can legitimately land one position
        // BELOW the array start in this mirrored coordinate frame — see the comment at
        // the `is_match && score > max_score` block after the loop.
        let mut max_offset = q_seed as i64;
        let mut score = 0i32;
        let mut is_match = false;
        let mut prev_match = false;

        for i in 0..q_len {
            let qi = search_q_top - i;
            let si = search_s_top - i;
            if qi >= query.len() || si >= subject.len() { continue; }
            is_match = query[qi] == subject[si];
            if is_match != prev_match {
                prev_match = is_match;
                if is_match {
                    score = 1;
                } else if score > max_score {
                    // qi is the first non-match going downward; run was above qi.
                    max_score = score;
                    max_offset = qi as i64 + (score / 2) as i64;
                }
            } else if is_match {
                score += 1;
                if score > hsp_max_ident_run {
                    // qi is the lowest match in a long run; take midpoint upward.
                    let new_q_seed = qi + (hsp_max_ident_run / 2) as usize;
                    let new_s_seed = search_s_top - (search_q_top - new_q_seed);
                    return (new_q_seed as u32, new_s_seed as u32);
                }
            }
        }
        if is_match && score > max_score {
            max_score = score;
            // NCBI (blast_gapalign.c:3477) computes this in signed Int4: after its loop
            // `index == q_start + q_len`, the EXCLUSIVE end of the scan, and
            // `max_offset = index - score/2` steps back into range.  This branch mirrors
            // that scan in reversed (FWD-genomic) coordinates, so the exclusive end lands
            // one position BELOW `prelim_q_start` — i.e. -1 when the prelim HSP starts at
            // query offset 0, which is common when a TE library is the query.  The
            // `+ score/2` then brings it back, exactly as NCBI's `- score/2` does.
            // Computing the intermediate in usize underflowed: a debug panic, and in
            // release a double wrap that happened to land on the correct value.
            let qi = search_q_top as i64 - q_len as i64;
            max_offset = qi + (score / 2) as i64;
        }
        if max_score > 0 {
            if max_offset < 0 {
                // Only reachable at score == 1 with the scan pinned to the array start.
                // NCBI's frame keeps this non-negative, so there is no faithful value to
                // mirror; fall back to NCBI's no-improvement outcome rather than emit a
                // wrapped seed (which is what the usize arithmetic used to return).
                crate::diag_count!(IMPROVE_SEED_NEGATIVE_OFFSET);
                return (q_seed as u32, s_seed as u32);
            }
            let delta = search_q_top as i64 - max_offset;
            let new_s_seed = (search_s_top as i64 - delta).max(0);
            (max_offset as u32, new_s_seed as u32)
        } else {
            (q_seed as u32, s_seed as u32)
        }
    }
}

/// Pre-built lookup table for a query sequence.
///
/// Build once per query with [`build_query_lookup`], then pass to
/// [`search_with_query_lookup`] for each subject.  This avoids rebuilding
/// the (potentially very large) LUT for every subject sequence.
pub struct QueryLookup {
    lookup: NaLookup,
    /// Forward query length (without leading/trailing sentinels).
    fwd_n: u32,
    /// For Combined mode: offset of context-1 (RC query) in the combined buffer (= fwd_n + 1).
    /// For SeparateStrands: 0 (unused).
    c1_start: u32,
}

/// Build the lookup table from the query sequence.
///
/// Call this once per query before the parallel loop over subjects.
/// The result is `Send + Sync` and can be shared across rayon threads.
/// Returns `(QueryLookup, masked_query)` where `masked_query` is the DUST-masked
/// version of `query` (or a clone when dust=false).  Pass `masked_query` to
/// `search_phase2a` so that word extension uses the same masked sequence that was
/// indexed in the lookup table.
pub fn build_query_lookup(query: &[u8], params: &SearchParams) -> (QueryLookup, Vec<u8>) {
    let fwd_n = query.len().saturating_sub(2);
    let fwd_query_scan = if params.dust {
        let mut v = query.to_vec();
        dust_mask(&mut v, DUST_WINDOW, DUST_LEVEL, DUST_LINKER, 14);
        v
    } else {
        query.to_vec()
    };
    let query_bases = &fwd_query_scan[1..1 + fwd_n];

    // approx_entries: Combined counts both strand contexts (2×fwd_n), matching
    // NCBI's BlastChooseNaLookupTable call in BlastNaLookupTableNew.
    let approx_entries = match params.seed_mode {
        SeedMode::Combined       => 2 * fwd_n,
        SeedMode::SeparateStrands => fwd_n,
    };
    let lut_width = choose_lut_width(params.word_size, approx_entries);

    let (lookup, c1_start) = match params.seed_mode {
        SeedMode::Combined => {
            // Build combined buffer: [fwd_bases | sentinel(14) | rc_bases].
            // This mirrors NCBI's BlastNaLookupTableNew combined-context layout:
            //   context-0 at [0, fwd_n-1], separator at fwd_n, context-1 at [fwd_n+1, 2*fwd_n].
            let rc_bases = revcomp_blastna(query_bases);
            let mut combined: Vec<u8> = Vec::with_capacity(2 * fwd_n + 1);
            combined.extend_from_slice(query_bases);
            combined.push(14u8); // N sentinel separates the two contexts
            combined.extend_from_slice(&rc_bases);
            let locs: &[(i32, i32)] = &[(0, fwd_n as i32 - 1), (fwd_n as i32 + 1, 2 * fwd_n as i32)];
            let lut = if lut_width > 8 {
                NaLookup::Mega(blast_mb_lookup_table_new(
                    &combined, locs,
                    params.word_size as i32, lut_width as i32, approx_entries,
                ))
            } else {
                NaLookup::Small(blast_na_lookup_table_new(
                    &combined, locs,
                    params.word_size as i32, lut_width as i32,
                ))
            };
            (lut, fwd_n as u32 + 1)
        }
        SeedMode::SeparateStrands => {
            let lut = if lut_width > 8 {
                NaLookup::Mega(blast_mb_lookup_table_new(
                    query_bases, &[(0, fwd_n as i32 - 1)],
                    params.word_size as i32, lut_width as i32, approx_entries,
                ))
            } else {
                NaLookup::Small(blast_na_lookup_table_new(
                    query_bases, &[(0, fwd_n as i32 - 1)],
                    params.word_size as i32, lut_width as i32,
                ))
            };
            (lut, 0u32)
        }
    };

    (QueryLookup { lookup, fwd_n: fwd_n as u32, c1_start }, fwd_query_scan)
}

/// Build a lookup table from a *pre-masked* query (DUST already applied).
/// Unlike [`build_query_lookup`], DUST masking is NOT re-applied here.
/// Use this when the caller has already applied full-query DUST masking and
/// extracted a chunk from it, to avoid double-masking artifacts at chunk boundaries.
pub fn build_query_lookup_premask(masked_query: &[u8], params: &SearchParams) -> QueryLookup {
    let fwd_n = masked_query.len().saturating_sub(2);
    let query_bases = &masked_query[1..1 + fwd_n];

    let approx_entries = match params.seed_mode {
        SeedMode::Combined       => 2 * fwd_n,
        SeedMode::SeparateStrands => fwd_n,
    };
    let lut_width = choose_lut_width(params.word_size, approx_entries);

    let (lookup, c1_start) = match params.seed_mode {
        SeedMode::Combined => {
            let rc_bases = revcomp_blastna(query_bases);
            let mut combined: Vec<u8> = Vec::with_capacity(2 * fwd_n + 1);
            combined.extend_from_slice(query_bases);
            combined.push(14u8);
            combined.extend_from_slice(&rc_bases);
            let locs: &[(i32, i32)] = &[(0, fwd_n as i32 - 1), (fwd_n as i32 + 1, 2 * fwd_n as i32)];
            let lut = if lut_width > 8 {
                NaLookup::Mega(blast_mb_lookup_table_new(
                    &combined, locs,
                    params.word_size as i32, lut_width as i32, approx_entries,
                ))
            } else {
                NaLookup::Small(blast_na_lookup_table_new(
                    &combined, locs,
                    params.word_size as i32, lut_width as i32,
                ))
            };
            (lut, fwd_n as u32 + 1)
        }
        SeedMode::SeparateStrands => {
            let lut = if lut_width > 8 {
                NaLookup::Mega(blast_mb_lookup_table_new(
                    query_bases, &[(0, fwd_n as i32 - 1)],
                    params.word_size as i32, lut_width as i32, approx_entries,
                ))
            } else {
                NaLookup::Small(blast_na_lookup_table_new(
                    query_bases, &[(0, fwd_n as i32 - 1)],
                    params.word_size as i32, lut_width as i32,
                ))
            };
            (lut, 0u32)
        }
    };

    QueryLookup { lookup, fwd_n: fwd_n as u32, c1_start }
}

impl QueryLookup {
    pub fn lut_word_length(&self) -> u32 {
        self.lookup.lut_word_length() as u32
    }
}

/// Apply DUST masking to a query sequence (with leading/trailing sentinels).
/// When dust=false, returns a clone of the input unchanged.
/// The masked query is used for seeding (LUT construction, word extension, ungapped scoring)
/// but NOT for Phase 2a/2b gapped alignment — NCBI uses mask_at_hash=TRUE for rmblastn.
pub fn mask_query_for_alignment(query: &[u8], params: &SearchParams) -> Vec<u8> {
    if params.dust {
        let mut v = query.to_vec();
        // Multi-chunk path: NCBI applies RestrictToSeqInt to every DUST interval,
        // which stores the right end as GetToOpen() = e+1 (off-by-one). Apply compat.
        dust_mask_ncbi_compat(&mut v, DUST_WINDOW, DUST_LEVEL, DUST_LINKER, 14);
        v
    } else {
        query.to_vec()
    }
}

/// Search one query (chunk) against one subject using a pre-built lookup table.
///
/// `query` is the genome chunk (BLASTNA-encoded with leading/trailing sentinels).
/// All q coordinates returned in AlignResult are absolute (chunk_offset added).
/// `chunk_offset` is the start of this chunk in the full genome (0 for whole-query mode).
/// Obtain `ql` by calling [`build_query_lookup`] once for the chunk.
pub fn search_with_query_lookup(
    ql: &QueryLookup,
    query: &[u8],
    query_id: &str,
    subject_plus: &[u8],
    subject_id: &str,
    params: &SearchParams,
    matrix: &ScoreMatrix,
    chunk_offset: u32,
    n_mask: &[u8],
    prepared: Option<&PreparedSubject>,
) -> Vec<AlignResult> {
    let q_len = query.len().saturating_sub(2) as u32;
    let s_len = subject_plus.len().saturating_sub(2) as u32;

    if q_len < params.word_size as u32 || s_len < params.word_size as u32 {
        return Vec::new();
    }

    // Subject RC + NCBI2NA-packed strands for the scan: borrow from the shared
    // per-subject cache when supplied, else compute locally (identical bytes).
    let local_prep;
    let prep: &PreparedSubject = match prepared {
        Some(p) => p,
        None => {
            local_prep = prepare_subject_strands(subject_plus, n_mask);
            &local_prep
        }
    };
    let (subj_rc, packed_plus, packed_minus): (&[u8], &[u8], &[u8]) =
        (&prep.rc, &prep.packed_plus, &prep.packed_minus);

    // Generate masked query for DUST: used for exact extension, ungapped scoring, Phase 2a, Phase 2b.
    // When dust=false, masked_query == query (no allocation).
    let masked_query_buf;
    let masked_query: &[u8] = if params.dust {
        let mut v = query.to_vec();
        dust_mask(&mut v, DUST_WINDOW, DUST_LEVEL, DUST_LINKER, 14);
        masked_query_buf = v;
        &masked_query_buf
    } else {
        query
    };

    let mut ungapped = Vec::new();
    if ql.c1_start > 0 {
        // Combined LUT: single scan of FWD subject, interleaved plus/minus hits.
        collect_ungapped_combined(
            query, masked_query, q_len, subject_plus, &packed_plus[3..], subj_rc, s_len,
            ql.c1_start, ql.fwd_n, &ql.lookup, params, matrix, &mut DiscardSeeds, &mut ungapped,
        );
    } else {
        collect_ungapped(
            query, masked_query, q_len, subject_plus, &packed_plus[3..], s_len,
            Strand::Plus, &ql.lookup, params, matrix, &mut DiscardSeeds, &mut ungapped,
        );
        collect_ungapped(
            query, masked_query, q_len, subj_rc, &packed_minus[3..], s_len,
            Strand::Minus, &ql.lookup, params, matrix, &mut DiscardSeeds, &mut ungapped,
        );
    }
    let mut results = run_gapped_phase(
        query, query_id, subject_plus, subj_rc, subject_id,
        q_len, s_len, ungapped, params, matrix, &mut DiscardUngap,
        ql.lut_word_length(), n_mask, prep,
    );
    if chunk_offset > 0 {
        for r in &mut results {
            r.hsp.q_start += chunk_offset;
            r.hsp.q_end   += chunk_offset;
        }
    }
    results
}

/// Like [`search_with_query_lookup`] but also collects post-deduplication seeds.
///
/// Seeds in `seed_out` are in chunk-relative coordinates (not adjusted by `chunk_offset`).
/// Use this for regression testing the seeding and deduplication stages.
pub fn search_with_query_lookup_seeds(
    ql: &QueryLookup,
    query: &[u8],
    query_id: &str,
    subject_plus: &[u8],
    subject_id: &str,
    params: &SearchParams,
    matrix: &ScoreMatrix,
    chunk_offset: u32,
    seed_out: &mut Vec<SeedRecord>,
    n_mask: &[u8],
    prepared: Option<&PreparedSubject>,
) -> Vec<AlignResult> {
    let q_len = query.len().saturating_sub(2) as u32;
    let s_len = subject_plus.len().saturating_sub(2) as u32;

    if q_len < params.word_size as u32 || s_len < params.word_size as u32 {
        return Vec::new();
    }

    // Subject RC + NCBI2NA-packed strands for the scan: borrow from the shared
    // per-subject cache when supplied, else compute locally (identical bytes).
    let local_prep;
    let prep: &PreparedSubject = match prepared {
        Some(p) => p,
        None => {
            local_prep = prepare_subject_strands(subject_plus, n_mask);
            &local_prep
        }
    };
    let (subj_rc, packed_plus, packed_minus): (&[u8], &[u8], &[u8]) =
        (&prep.rc, &prep.packed_plus, &prep.packed_minus);

    let masked_query_buf_s;
    let masked_query_s: &[u8] = if params.dust {
        let mut v = query.to_vec();
        dust_mask(&mut v, DUST_WINDOW, DUST_LEVEL, DUST_LINKER, 14);
        masked_query_buf_s = v;
        &masked_query_buf_s
    } else {
        query
    };

    let mut ungapped = Vec::new();
    if ql.c1_start > 0 {
        collect_ungapped_combined(
            query, masked_query_s, q_len, subject_plus, &packed_plus[3..], subj_rc, s_len,
            ql.c1_start, ql.fwd_n, &ql.lookup, params, matrix, seed_out, &mut ungapped,
        );
    } else {
        collect_ungapped(
            query, masked_query_s, q_len, subject_plus, &packed_plus[3..], s_len,
            Strand::Plus, &ql.lookup, params, matrix, seed_out, &mut ungapped,
        );
        collect_ungapped(
            query, masked_query_s, q_len, subj_rc, &packed_minus[3..], s_len,
            Strand::Minus, &ql.lookup, params, matrix, seed_out, &mut ungapped,
        );
    }
    let mut results = run_gapped_phase(
        query, query_id, subject_plus, subj_rc, subject_id,
        q_len, s_len, ungapped, params, matrix, &mut DiscardUngap,
        ql.lut_word_length(), n_mask, prep,
    );
    if chunk_offset > 0 {
        for r in &mut results {
            r.hsp.q_start += chunk_offset;
            r.hsp.q_end   += chunk_offset;
        }
    }
    results
}

/// Like [`search_with_query_lookup`] but also collects post-sort ungapped hits.
///
/// Hits in `ungap_out` are in FWD-normalized chunk-relative coordinates and in
/// score-descending order (same order they enter gapped alignment).
pub fn search_with_query_lookup_ungapped(
    ql: &QueryLookup,
    query: &[u8],
    query_id: &str,
    subject_plus: &[u8],
    subject_id: &str,
    params: &SearchParams,
    matrix: &ScoreMatrix,
    chunk_offset: u32,
    ungap_out: &mut Vec<UngappedHit>,
    n_mask: &[u8],
    prepared: Option<&PreparedSubject>,
) -> Vec<AlignResult> {
    let q_len = query.len().saturating_sub(2) as u32;
    let s_len = subject_plus.len().saturating_sub(2) as u32;

    if q_len < params.word_size as u32 || s_len < params.word_size as u32 {
        return Vec::new();
    }

    // Subject RC + NCBI2NA-packed strands for the scan: borrow from the shared
    // per-subject cache when supplied, else compute locally (identical bytes).
    let local_prep;
    let prep: &PreparedSubject = match prepared {
        Some(p) => p,
        None => {
            local_prep = prepare_subject_strands(subject_plus, n_mask);
            &local_prep
        }
    };
    let (subj_rc, packed_plus, packed_minus): (&[u8], &[u8], &[u8]) =
        (&prep.rc, &prep.packed_plus, &prep.packed_minus);

    let masked_query_buf_u;
    let masked_query_u: &[u8] = if params.dust {
        let mut v = query.to_vec();
        dust_mask(&mut v, DUST_WINDOW, DUST_LEVEL, DUST_LINKER, 14);
        masked_query_buf_u = v;
        &masked_query_buf_u
    } else {
        query
    };

    let mut ungapped = Vec::new();
    if ql.c1_start > 0 {
        collect_ungapped_combined(
            query, masked_query_u, q_len, subject_plus, &packed_plus[3..], subj_rc, s_len,
            ql.c1_start, ql.fwd_n, &ql.lookup, params, matrix, &mut DiscardSeeds, &mut ungapped,
        );
    } else {
        collect_ungapped(
            query, masked_query_u, q_len, subject_plus, &packed_plus[3..], s_len,
            Strand::Plus, &ql.lookup, params, matrix, &mut DiscardSeeds, &mut ungapped,
        );
        collect_ungapped(
            query, masked_query_u, q_len, subj_rc, &packed_minus[3..], s_len,
            Strand::Minus, &ql.lookup, params, matrix, &mut DiscardSeeds, &mut ungapped,
        );
    }
    let mut results = run_gapped_phase(
        query, query_id, subject_plus, subj_rc, subject_id,
        q_len, s_len, ungapped, params, matrix, ungap_out,
        ql.lut_word_length(), n_mask, prep,
    );
    if chunk_offset > 0 {
        for r in &mut results {
            r.hsp.q_start += chunk_offset;
            r.hsp.q_end   += chunk_offset;
        }
    }
    results
}

/// Search one query against one subject sequence (plus and minus strands).
///
/// Builds the lookup table on every call.  Use [`build_query_lookup`] +
/// [`search_with_query_lookup`] when searching multiple subjects against the
/// same query to avoid rebuilding the LUT for each subject.
///
/// `query` and `subject_plus` are BLASTNA-encoded with leading/trailing sentinels.
pub fn search_query_vs_subject(
    query: &[u8],
    query_id: &str,
    subject_plus: &[u8],
    subject_id: &str,
    params: &SearchParams,
    matrix: &ScoreMatrix,
) -> Vec<AlignResult> {
    let (ql, _) = build_query_lookup(query, params);
    search_with_query_lookup(&ql, query, query_id, subject_plus, subject_id, params, matrix, 0, &[], None)
}

/// Filter HSPs by masklevel: an HSP is dropped if any single higher-scoring HSP
/// covers more than `mask_level`% of its query span.
///
/// `subject_names` controls sort order and tiebreaking:
/// - Non-empty (cross-subject mode): mirrors NCBI's `Blast_HSPResultsApplyMasklevel` —
///   all HSPs from all subjects are sorted together by (score DESC, oid DESC) where oid
///   is the index of the subject in `subject_names`.  Must be called once after all
///   subjects' results are collected for a query.
/// - Empty (per-subject / unit-test mode): uses the per-subject ScoreCompareHSPs sort
///   (score DESC, s_start ASC, s_end DESC, q_off ASC, q_end DESC).
pub fn apply_mask_level(results: &mut Vec<AlignResult>, mask_level: u32, subject_names: &[String]) {
    if mask_level >= 100 {
        return;
    }
    if subject_names.is_empty() {
        // Per-subject sort: ScoreCompareHSPs order.
        // For minus strand, NCBI query.offset = qlen - FWD_qe (ASC → FWD_qe DESC)
        // and query.end = qlen - FWD_qs (DESC → FWD_qs ASC).
        let qa_len = results.first().map(|r| r.hsp.q_len).unwrap_or(0);
        results.sort_unstable_by(|a, b| {
            let q_off_a = match a.hsp.strand { Strand::Plus => a.hsp.q_start, Strand::Minus => qa_len - a.hsp.q_end };
            let q_off_b = match b.hsp.strand { Strand::Plus => b.hsp.q_start, Strand::Minus => qa_len - b.hsp.q_end };
            let q_end_a = match a.hsp.strand { Strand::Plus => a.hsp.q_end, Strand::Minus => qa_len - a.hsp.q_start };
            let q_end_b = match b.hsp.strand { Strand::Plus => b.hsp.q_end, Strand::Minus => qa_len - b.hsp.q_start };
            b.hsp.score.cmp(&a.hsp.score)
                .then_with(|| a.hsp.s_start.cmp(&b.hsp.s_start))
                .then_with(|| b.hsp.s_end.cmp(&a.hsp.s_end))
                .then_with(|| q_off_a.cmp(&q_off_b))
                .then_with(|| q_end_b.cmp(&q_end_a))
        });
    } else {
        // Cross-subject sort: mirrors s_SortHspWrapRawScore in blast_hits.c.
        // Primary: score DESC, oid DESC.
        // Tiebreakers: ScoreCompareHSPs (blast_hits.c), which is the per-subject
        // pre-sort order.  NCBI uses a stable sort for s_SortHspWrapRawScore, so
        // the ScoreCompareHSPs pre-order is preserved for equal-score, equal-oid
        // hits.  Adding these tiebreakers replicates that stable-sort behavior.
        let oid_map: std::collections::HashMap<&str, usize> = subject_names.iter()
            .enumerate()
            .map(|(i, n)| (n.as_str(), i))
            .collect();
        results.sort_unstable_by(|a, b| {
            let oid_a = oid_map.get(a.subject_id.as_str()).copied().unwrap_or(0);
            let oid_b = oid_map.get(b.subject_id.as_str()).copied().unwrap_or(0);
            b.hsp.score.cmp(&a.hsp.score)
                .then_with(|| oid_b.cmp(&oid_a))
                // ScoreCompareHSPs tiebreakers: subject.offset ASC, subject.end DESC
                .then_with(|| a.hsp.s_start.cmp(&b.hsp.s_start))
                .then_with(|| b.hsp.s_end.cmp(&a.hsp.s_end))
                // query.offset ASC: plus=q_start, minus=q_len-q_end
                .then_with(|| {
                    let q_off_a = if a.hsp.strand == Strand::Plus { a.hsp.q_start } else { a.hsp.q_len - a.hsp.q_end };
                    let q_off_b = if b.hsp.strand == Strand::Plus { b.hsp.q_start } else { b.hsp.q_len - b.hsp.q_end };
                    q_off_a.cmp(&q_off_b)
                })
                // query.end DESC: plus=q_end, minus=q_len-q_start
                .then_with(|| {
                    let q_end_a = if a.hsp.strand == Strand::Plus { a.hsp.q_end } else { a.hsp.q_len - a.hsp.q_start };
                    let q_end_b = if b.hsp.strand == Strand::Plus { b.hsp.q_end } else { b.hsp.q_len - b.hsp.q_start };
                    q_end_b.cmp(&q_end_a)
                })
        });
    }

    let n = results.len();
    if n == 0 { return; }

    // (q_start, q_end, score) — sorted by q_start for binary search.
    let mut accepted: Vec<(u32, u32, i32)> = Vec::new();
    let mut max_accepted_span: u32 = 0;
    let mut keep = vec![true; n];

    for i in 0..n {
        // Replicate NCBI's minus-strand coordinate shift in blast_itree.c.
        // When eQueryOnlyStrandIndifferent is used, minus-strand HSPs are
        // converted using context_start = genome_len+1, so their FWD
        // coordinates are +1 relative to the true 0-based positions.
        // This makes cross-strand overlap exactly match NCBI at boundary cases.
        let (qs_i, qe_i) = if results[i].hsp.strand == Strand::Minus {
            (results[i].hsp.q_start + 1, results[i].hsp.q_end + 1)
        } else {
            (results[i].hsp.q_start, results[i].hsp.q_end)
        };
        let sc_i = results[i].hsp.score;
        let q_span_i = (qe_i - qs_i) as i64;
        if q_span_i == 0 {
            keep[i] = false;
            continue;
        }

        let upper = accepted.partition_point(|a| a.0 < qe_i);
        let min_start = qs_i.saturating_sub(max_accepted_span);

        let mut masked = false;
        let mut idx = upper;
        while idx > 0 {
            idx -= 1;
            let (qs_j, qe_j, sc_j) = accepted[idx];
            if qs_j < min_start {
                break;
            }
            if qe_j <= qs_i {
                continue;
            }
            // Mirror NCBI s_HSPQueryRangeIsMasklevelContained: skip only if accepted
            // entry has strictly lower score (in_score > tree_hsp->score → return 0).
            if sc_j < sc_i {
                continue;
            }
            let ovlp = (qe_j.min(qe_i) as i64) - (qs_j.max(qs_i) as i64);
            // NCBI: (Int4)(100*(double)ovlp/span) >= masklevel
            if ovlp * 100 / q_span_i >= mask_level as i64 {
                masked = true;
                break;
            }
        }

        if masked {
            keep[i] = false;
        } else {
            let pos = accepted.partition_point(|a| a.0 < qs_i);
            accepted.insert(pos, (qs_i, qe_i, sc_i));
            let span = (qe_i - qs_i) as u32;
            if span > max_accepted_span {
                max_accepted_span = span;
            }
        }
    }

    let mut j = 0;
    results.retain(|_| { let k = keep[j]; j += 1; k });
}

/// Puts one query's results in the order NCBI hands to its formatters.
///
/// After traceback NCBI sorts the hit list with `s_EvalueCompareHSPLists`
/// (best e-value ASC, top score DESC, oid DESC) and each subject's HSPs with
/// `s_EvalueCompareHSPs` (e-value ASC, then `ScoreCompareHSPs`), and both the
/// tabular and the pairwise formatter print in that order: every HSP of one
/// subject, then the next subject.  Within one query the e-value is a monotone
/// function of the score (and a constant sentinel for table matrices), so the
/// e-value keys never reorder anything the score keys do not, and this
/// function sorts on the score keys only.  `oid` is the index in `subject_names`; unknown names get
/// oid 0, as in `apply_mask_level`.
pub fn sort_hit_list_order(results: &mut [AlignResult], subject_names: &[String]) {
    let oid_map: std::collections::HashMap<&str, usize> = subject_names.iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    // (top score, oid) per subject, keyed by owned name so the map does not
    // borrow `results` while we sort it.
    let mut rank: std::collections::HashMap<String, (i32, usize)> = std::collections::HashMap::new();
    for r in results.iter() {
        let oid = oid_map.get(r.subject_id.as_str()).copied().unwrap_or(0);
        let e = rank.entry(r.subject_id.clone()).or_insert((i32::MIN, oid));
        e.0 = e.0.max(r.hsp.score);
    }
    results.sort_by(|a, b| {
        let (best_a, oid_a) = rank[a.subject_id.as_str()];
        let (best_b, oid_b) = rank[b.subject_id.as_str()];
        best_b.cmp(&best_a)
            .then_with(|| oid_b.cmp(&oid_a))
            .then_with(|| score_compare_hsps(&a.hsp, &b.hsp))
    });
}

/// Phase 1: scan and collect all passing ungapped hits for one strand.
///
/// Mirrors NCBI's BlastNaWordFinder loop: seed scan → exact extension →
/// diagonal tracker → ungapped extension → save.
///
/// `subject`        — BLASTNA (1 byte/base) with sentinels; used for ungapped extension.
/// `packed_subject` — NCBI2NA packed (4 bases/byte) starting at base-0 byte
///                    (i.e., `seq_blk.compressed_nuc_seq()` = `packed[3..]`).
fn collect_ungapped<S: SeedOut>(
    query: &[u8],
    ext_query: &[u8],
    q_len: u32,
    subject: &[u8],
    packed_subject: &[u8],
    s_len: u32,
    strand: Strand,
    lookup: &NaLookup,
    params: &SearchParams,
    matrix: &ScoreMatrix,
    seed_out: &mut S,
    ungapped: &mut Vec<UngappedRecord>,
) {
    let _ = (q_len, matrix);

    let s_real = subject.len().saturating_sub(2);
    let lw   = lookup.lut_word_length() as usize;
    let step = lookup.scan_step() as usize;

    // scan_start: for minus strand match NCBI's context-1 start offset so
    // that covered subject positions align with the plus-strand FWD scan.
    let scan_start = match strand {
        Strand::Plus  => 0usize,
        Strand::Minus => {
            if step <= 1 || s_real < lw { 0 } else { (s_real - lw) % step }
        }
    };

    let end_range = s_real.saturating_sub(lw) as i32;
    if scan_start as i32 > end_range {
        return;
    }

    let ext_needed = params.word_size as i64 - lookup.lut_word_length() as i64;
    let q_real_len = query.len().saturating_sub(2) as i64;
    let s_real_len = subject.len().saturating_sub(2) as i64;
    let lut_w      = lookup.lut_word_length() as i64;
    // KA-derived cutoff when the matrix+gap combo has ALP params (set per-search
    // in search_db_parallel); otherwise the historical fixed min_raw_gapped_score/2.
    let ungapped_cutoff = params
        .ungapped_cutoff
        .unwrap_or(params.min_raw_gapped_score / 2);

    // Batch must exceed longest_chain so max_hits = batch_size - longest_chain > 0.
    // A batch smaller than longest_chain causes the scan to break with 0 hits at the
    // first hit position and never advance scan_range — an infinite loop.
    let batch_size = (4096usize).max(lookup.longest_chain() as usize + 4096);
    let mut batch  = vec![BlastOffsetPair::default(); batch_size];

    // -----------------------------------------------------------------------
    // Plus strand: scan ascending FWD_s → dedup in ascending FWD_s (correct).
    //
    // Minus strand: scan ascending RC_s = descending FWD_s.  NCBI processes
    // minus-strand hits in ascending FWD_s order (leftmost-FWD seed wins on
    // each diagonal).  To reproduce that without a global O(n log n) sort:
    //   • Pass 1 (scan loop): collect all exactly-extended seeds into a Vec.
    //   • Pass 2: iterate the Vec in REVERSE (= descending RC_s = ascending FWD_s)
    //     and apply the same FWD-ascending dedup + ungapped extension.
    // This is O(n) vs the old O(n log n) sort while giving identical results.
    // -----------------------------------------------------------------------
    let mut minus_seeds: Vec<(u32, u32)> = if strand == Strand::Minus { Vec::new() } else { Vec::new() };

    let mut diag_hash = BlastDiagHash::new();

    let mut scan_range = [scan_start as i32, end_range];
    loop {
        let n = lookup.scan(packed_subject, &mut batch, batch_size as i32, &mut scan_range);

        for i in 0..n as usize {
            let q0 = batch[i].q_off as i64;
            let s0 = batch[i].s_off as i64;

            // Exact extension: extend the lut_width-mer to a full word_size match.
            let (q_off, s_off) = if ext_needed > 0 {
                // NCBI's combined-LUT context-1 scans the forward subject and extends
                // "left" in the combined buffer first.  For context-1 (RC query), "left"
                // in the combined buffer corresponds to going RIGHT in forward-query space.
                // Rust scans RC(subject) with the forward LUT, so "left" in Rust's scan
                // means going LEFT in forward-query space — the opposite direction.
                // Mirror NCBI: for minus strand try right (= NCBI ext_left, anchor-shifting)
                // first; for plus strand try left first.
                // NCBI: mask_at_hash=TRUE — masked query for exact word extension.
                let (ext_left, ext_right) = if strand == Strand::Minus {
                    // Minus strand: right extension (higher FWD-q, lower FWD-s) mirrors NCBI's
                    // ext_left (left in ctx1 = right in FWD-q).  Try it first so that when
                    // bases are available on both sides we land on the same seed as NCBI.
                    let mut ext_right = 0i64;
                    while ext_right < ext_needed {
                        let qi = q0 + lut_w + ext_right;
                        let si = s0 + lut_w + ext_right;
                        if qi >= q_real_len || si >= s_real_len { break; }
                        let qb = ext_query[1 + qi as usize];
                        let sb = subject[1 + si as usize];
                        if qb >= 4 || qb != sb { break; }
                        ext_right += 1;
                    }
                    let mut ext_left = 0i64;
                    while ext_left + ext_right < ext_needed {
                        let qi = q0 - ext_left - 1;
                        let si = s0 - ext_left - 1;
                        if qi < 0 || si < 0 { break; }
                        let qb = ext_query[1 + qi as usize];
                        let sb = subject[1 + si as usize];
                        if qb >= 4 || qb != sb { break; }
                        ext_left += 1;
                    }
                    (ext_left, ext_right)
                } else {
                    // Plus strand: left extension (lower FWD-q, lower FWD-s) mirrors NCBI ext_left.
                    let mut ext_left = 0i64;
                    while ext_left < ext_needed {
                        let qi = q0 - ext_left - 1;
                        let si = s0 - ext_left - 1;
                        if qi < 0 || si < 0 { break; }
                        let qb = ext_query[1 + qi as usize];
                        let sb = subject[1 + si as usize];
                        if qb >= 4 || qb != sb { break; }
                        ext_left += 1;
                    }
                    let mut ext_right = 0i64;
                    while ext_left + ext_right < ext_needed {
                        let qi = q0 + lut_w + ext_right;
                        let si = s0 + lut_w + ext_right;
                        if qi >= q_real_len || si >= s_real_len { break; }
                        let qb = ext_query[1 + qi as usize];
                        let sb = subject[1 + si as usize];
                        if qb >= 4 || qb != sb { break; }
                        ext_right += 1;
                    }
                    (ext_left, ext_right)
                };
                if ext_left + ext_right < ext_needed { continue; }
                // Anchor position: plus strand shifts left by ext_left; minus strand
                // shifts right by ext_right (= NCBI ext_left, the anchor-shifting direction).
                // The stored (q_off, s_off) is the anchor used for diagonal tracking and
                // ungapped extension; FWD-normalized seed left end = q_off - ext_needed (minus)
                // or q_off (plus).
                if strand == Strand::Minus {
                    ((q0 + ext_right) as u32, (s0 + ext_right) as u32)
                } else {
                    ((q0 - ext_left) as u32, (s0 - ext_left) as u32)
                }
            } else {
                (batch[i].q_off, batch[i].s_off)
            };

            match strand {
                Strand::Minus => {
                    // Defer dedup/extension: collect seeds for reverse-order pass below.
                    minus_seeds.push((q_off, s_off));
                }
                Strand::Plus => {
                    // Plus: ascending FWD_s — dedup and extend inline.
                    let diag = (s_off as i32).wrapping_sub(q_off as i32);
                    let last_se = diag_hash.get(diag);
                    if s_off < last_se { continue; }

                    seed_out.push_seed(SeedRecord { q_off, s_off, strand: Strand::Plus });
                    crate::diag_count!(COUNT_SEEDS);
                    let raw_p = extend_ungapped(query, subject, q_off, s_off, matrix, params.xdrop_ungap, i32::MIN, true).unwrap();
                    let s_end_store = if raw_p.score >= ungapped_cutoff { raw_p.s_end } else { s_off + params.word_size as u32 };
                    diag_hash.insert(diag, s_end_store, s_off, 1 - params.word_size as i32);
                    if raw_p.score < ungapped_cutoff { continue; }
                    let ung = raw_p;
                    crate::diag_count!(COUNT_UNGAPPED_HITS);
                    ungapped.push(UngappedRecord {
                        q_seed: q_off, s_seed: s_off,
                        q_start: ung.q_start, q_end: ung.q_end,
                        s_start: ung.s_start, s_end: ung.s_end,
                        score: ung.score, strand,
                    });
                }
            }
        }

        if scan_range[0] > end_range { break; }
    }

    // Minus-strand pass 2: process collected seeds in reverse (descending RC_s =
    // ascending FWD_s), mirroring the old sorted-by-(diag,-s_off) dedup order.
    // diag_hash tracks rightmost FWD_s_end covered on each diagonal (FWD-space dedup).
    if strand == Strand::Minus {
        // Extend in NCBI's native (RC_query, FWD_subject) frame — see the detailed note
        // in collect_ungapped_combined.  `subject` here is the RC subject (subj_rc); its
        // revcomp is the forward subject.  `query` is forward; its revcomp is RC_query.
        let q_len_local = query.len().saturating_sub(2) as u32;
        let query_rc_full = revcomp_blastna(query);
        let subject_fwd = revcomp_blastna(subject);
        for &(q_off, s_off) in minus_seeds.iter().rev() {
            let diag = (s_off as i32).wrapping_sub(q_off as i32);
            // FWD-space coords: s_off is RC position; FWD right end = s_len - s_off.
            let fwd_s_start = s_len.saturating_sub(s_off + params.word_size as u32);
            let fwd_s_end   = s_len - s_off;
            let last_fwd_se = diag_hash.get(diag);
            if fwd_s_start < last_fwd_se { continue; }

            let s_seed = s_len.saturating_sub(s_off + params.word_size as u32);
            seed_out.push_seed(SeedRecord { q_off, s_off: s_seed, strand: Strand::Minus });
            crate::diag_count!(COUNT_SEEDS);
            let qn = q_len_local - 1 - q_off;
            let sn = s_len - 1 - s_off;
            let raw_ncbi = extend_ungapped(&query_rc_full, &subject_fwd, qn, sn, matrix, params.xdrop_ungap, i32::MIN, true).unwrap();
            let raw_m = UngappedResult {
                score:   raw_ncbi.score,
                q_start: q_len_local - raw_ncbi.q_end,
                q_end:   q_len_local - raw_ncbi.q_start,
                s_start: s_len - raw_ncbi.s_end,
                s_end:   s_len - raw_ncbi.s_start,
            };
            let fwd_se_store = if raw_ncbi.score >= ungapped_cutoff { raw_ncbi.s_end } else { fwd_s_end };
            diag_hash.insert(diag, fwd_se_store, fwd_s_start, 1 - params.word_size as i32);
            if raw_m.score < ungapped_cutoff { continue; }
            let ung = raw_m;
            crate::diag_count!(COUNT_UNGAPPED_HITS);
            ungapped.push(UngappedRecord {
                q_seed: q_off, s_seed: s_off,
                q_start: ung.q_start, q_end: ung.q_end,
                s_start: ung.s_start, s_end: ung.s_end,
                score: ung.score, strand,
            });
        }
    }
}

/// Combined-LUT Phase 1: single FWD-subject scan, interleaved plus/minus hits.
///
/// Mirrors NCBI's BlastNaWordFinder with a combined (two-context) lookup table:
///   context-0 (q_off < c1_start): plus-strand seed, FWD_query vs FWD_subject.
///   context-1 (q_off >= c1_start): minus-strand seed, RC_query vs FWD_subject.
///     Converted to (FWD_query, RC_subject) space for exact-extension and ungapped
///     extension so the same extend_ungapped() call is reused.
///
/// Because both contexts are returned in ascending FWD_s order, both can be deduped
/// inline with separate diagonal maps — no collection+reversal pass needed for minus.
#[allow(clippy::too_many_arguments)]
fn collect_ungapped_combined<S: SeedOut>(
    query: &[u8],
    ext_query: &[u8],
    q_len: u32,
    subject_plus: &[u8],
    packed_fwd_subject: &[u8],
    subj_rc: &[u8],
    s_len: u32,
    c1_start: u32,
    fwd_n: u32,
    lookup: &NaLookup,
    params: &SearchParams,
    matrix: &ScoreMatrix,
    seed_out: &mut S,
    ungapped: &mut Vec<UngappedRecord>,
) {
    let _ = (q_len, matrix);

    let s_real = subject_plus.len().saturating_sub(2);
    let lw    = lookup.lut_word_length() as i64;
    let step  = lookup.scan_step() as usize;

    let end_range = s_real.saturating_sub(lw as usize) as i32;
    if end_range < 0 { return; }

    let ext_needed    = params.word_size as i64 - lw;
    let q_real_len    = query.len().saturating_sub(2) as i64;
    let s_real_len    = s_len as i64; // same for FWD and RC subject
    // RC of the FULL query (sentinel-wrapped), so minus-strand ungapped extension can be
    // performed in NCBI's native (RC_query, FWD_subject) frame.  Extending in the mirror
    // (FWD_query, RC_subject) frame is score-equivalent but NOT decision-equivalent: the
    // adaptive-xdrop right extension scores bases in a different order, so X_current
    // tightens at different points and the committed endpoint can differ (bug #36).
    let query_rc_full = revcomp_blastna(query);
    // KA-derived cutoff when the matrix+gap combo has ALP params (set per-search
    // in search_db_parallel); otherwise the historical fixed min_raw_gapped_score/2.
    let ungapped_cutoff = params
        .ungapped_cutoff
        .unwrap_or(params.min_raw_gapped_score / 2);
    let c1 = c1_start as i64;
    let fn_ = fwd_n as i64;

    let batch_size = (4096usize).max(lookup.longest_chain() as usize + 4096);
    let mut batch  = vec![BlastOffsetPair::default(); batch_size];

    // Single shared hash table for both strands, matching NCBI's ewp->hash_table which
    // covers both context-0 (plus) and context-1 (minus) seeds in one 512-bucket table.
    let mut diag_hash = BlastDiagHash::new();

    let mut scan_range = [0i32, end_range];
    loop {
        let n = lookup.scan(packed_fwd_subject, &mut batch, batch_size as i32, &mut scan_range);

        for i in 0..n as usize {
            let q0_combined = batch[i].q_off as i64;
            let s_fwd       = batch[i].s_off as i64;

            if q0_combined < c1 {
                // ── Context-0: plus strand ────────────────────────────────────
                // NCBI: mask_at_hash=TRUE — masked query for exact word extension.
                // (Seeding uses masked query; extension also uses masked query;
                // ungapped scoring uses unmasked query — same as NCBI's behavior.)
                let (q_off, s_off) = if ext_needed > 0 {
                    let mut ext_left = 0i64;
                    while ext_left < ext_needed {
                        let qi = q0_combined - ext_left - 1;
                        let si = s_fwd - ext_left - 1;
                        if qi < 0 || si < 0 { break; }
                        let qb = ext_query[1 + qi as usize];
                        let sb = subject_plus[1 + si as usize];
                        if qb >= 4 || sb >= 4 || qb != sb { break; }
                        ext_left += 1;
                    }
                    let mut ext_right = 0i64;
                    while ext_left + ext_right < ext_needed {
                        let qi = q0_combined + lw + ext_right;
                        let si = s_fwd + lw + ext_right;
                        if qi >= q_real_len || si >= s_real_len { break; }
                        let qb = ext_query[1 + qi as usize];
                        let sb = subject_plus[1 + si as usize];
                        if qb >= 4 || sb >= 4 || qb != sb { break; }
                        ext_right += 1;
                    }
                    if ext_left + ext_right < ext_needed { continue; }
                    ((q0_combined - ext_left) as u32, (s_fwd - ext_left) as u32)
                } else {
                    (batch[i].q_off, batch[i].s_off)
                };

                // NCBI sign convention: diag = s_off - q_off (i32, wrapping).
                // ext_left cancels: (s_fwd - ext_left) - (q0_combined - ext_left) = s_fwd - q0_combined.
                let diag = (s_off as i32).wrapping_sub(q_off as i32);
                let last_se = diag_hash.get(diag);
                // Use extended position s_off for dedup check, matching NCBI's s_BlastNaExtend
                // which passes the left-shifted s_offset to s_BlastnDiagHashExtendInitialHit.
                if s_off < last_se { continue; }

                if crate::diag_enabled!("BLAST_DUMP_SEEDS") {
                    eprintln!("SEED strand=+ q={} s={}", q_off, s_off);
                }
                seed_out.push_seed(SeedRecord { q_off, s_off, strand: Strand::Plus });
                crate::diag_count!(COUNT_SEEDS);
                // NCBI: mask_at_hash=TRUE — unmasked query for ungapped extension scoring.
                let raw_plus = extend_ungapped(query, subject_plus, q_off, s_off, matrix, params.xdrop_ungap, i32::MIN, true).unwrap();
                let s_end_store = if raw_plus.score >= ungapped_cutoff { raw_plus.s_end } else { s_off + params.word_size as u32 };
                diag_hash.insert(diag, s_end_store, s_off, 1 - params.word_size as i32);
                if raw_plus.score < ungapped_cutoff { continue; }
                let ung = raw_plus;
                crate::diag_count!(COUNT_UNGAPPED_HITS);
                ungapped.push(UngappedRecord {
                    q_seed: q_off, s_seed: s_off,
                    q_start: ung.q_start, q_end: ung.q_end,
                    s_start: ung.s_start, s_end: ung.s_end,
                    score: ung.score, strand: Strand::Plus,
                });
            } else {
                // ── Context-1: minus strand ───────────────────────────────────
                // RC_query[rc_q_off .. rc_q_off+lw] matched FWD_subject[s_fwd .. s_fwd+lw].
                // Convert to (FWD_query, RC_subject) space so extend_ungapped() can be reused.
                //   q0 = fwd_n - rc_q_off - lw  (left end in FWD_query of the lw match)
                //   s0 = s_len - s_fwd - lw      (left end in RC_subject)
                let rc_q_off = q0_combined - c1;
                let q0 = fn_ - rc_q_off - lw;
                let s0 = s_real_len - s_fwd - lw;
                if q0 < 0 || s0 < 0 { continue; }

                // Exact extension in (FWD_q, RC_s) space: try RIGHT first, then LEFT.
                // NCBI: mask_at_hash=TRUE — masked query for exact word extension.
                let (q_off, s_off) = if ext_needed > 0 {
                    let mut ext_right = 0i64;
                    while ext_right < ext_needed {
                        let qi = q0 + lw + ext_right;
                        let si = s0 + lw + ext_right;
                        if qi >= q_real_len || si >= s_real_len { break; }
                        let qb = ext_query[1 + qi as usize];
                        let sb = subj_rc[1 + si as usize];
                        if qb >= 4 || sb >= 4 || qb != sb { break; }
                        ext_right += 1;
                    }
                    let mut ext_left = 0i64;
                    while ext_left + ext_right < ext_needed {
                        let qi = q0 - ext_left - 1;
                        let si = s0 - ext_left - 1;
                        if qi < 0 || si < 0 { break; }
                        let qb = ext_query[1 + qi as usize];
                        let sb = subj_rc[1 + si as usize];
                        if qb >= 4 || sb >= 4 || qb != sb { break; }
                        ext_left += 1;
                    }
                    if ext_left + ext_right < ext_needed { continue; }
                    ((q0 + ext_right) as u32, (s0 + ext_right) as u32)
                } else {
                    (q0 as u32, s0 as u32)
                };

                // Diagonal dedup in FWD-subject space (ascending s_fwd order).
                // NCBI sign convention: diag = s_fwd - q_combined (combined-buffer / FWD-subject
                // space, same as NCBI's s_BlastnDiagHashExtendInitialHit line 872: diag = s_off - q_off).
                let diag = (s_fwd as i32).wrapping_sub(q0_combined as i32);
                // FWD_s position after extension: ext_right in RC_s == ext_left in FWD_s,
                // so s_fwd_eff = s_fwd - ext_right = s_len - s_off - lw.
                let s_fwd_eff = s_len - s_off - lw as u32;
                // fwd_s_end: NCBI stores s_off + word_length (not lut_word_length) for failed seeds.
                let fwd_s_end = s_fwd_eff + params.word_size as u32;
                let last_fwd_se = diag_hash.get(diag);
                if s_fwd_eff < last_fwd_se { continue; }

                let q_seed = q_off.saturating_sub(ext_needed.max(0) as u32);
                let s_seed = s_len.saturating_sub(s_off + lw as u32);
                if crate::diag_enabled!("BLAST_DUMP_SEEDS") {
                    eprintln!("SEED strand=- q={} s={}", q_seed, s_seed);
                }
                seed_out.push_seed(SeedRecord { q_off: q_seed, s_off: s_seed, strand: Strand::Minus });
                crate::diag_count!(COUNT_SEEDS);
                // Ungapped extension in NCBI's native (RC_query, FWD_subject) frame.
                // Map the (FWD_q, RC_s) anchor to (RC_q, FWD_s): qn = q_len-1-q_off,
                // sn = s_len-1-s_off.  Using the plus-strand code (left=fixed, right=adaptive)
                // is what makes this faithful: NCBI's adaptive RIGHT extension begins AT the
                // seed base, so the matching word commits positive score and tightens
                // X_current before the extension crosses any low-scoring region.  The old
                // (FWD_q, RC_s) path ran the adaptive extension starting just PAST the word,
                // so X_current never tightened from the word and the extension over-ran into
                // adjacent regions, diverging from NCBI and dedup-suppressing real seeds (#36).
                let qn = q_len - 1 - q_off;
                let sn = s_len - 1 - s_off;
                let raw_ncbi = extend_ungapped(&query_rc_full, subject_plus, qn, sn, matrix, params.xdrop_ungap, i32::MIN, true).unwrap();
                // Map (RC_q, FWD_s) result back to (FWD_q, RC_s) for UngappedRecord storage.
                let raw_minus = UngappedResult {
                    score:   raw_ncbi.score,
                    q_start: q_len - raw_ncbi.q_end,
                    q_end:   q_len - raw_ncbi.q_start,
                    s_start: s_len - raw_ncbi.s_end,
                    s_end:   s_len - raw_ncbi.s_start,
                };
                // fwd_se_store: FWD-subject end of the ungapped extension (= raw_ncbi.s_end).
                let fwd_se_store = if raw_ncbi.score >= ungapped_cutoff { raw_ncbi.s_end } else { fwd_s_end };
                diag_hash.insert(diag, fwd_se_store, s_fwd_eff, 1 - params.word_size as i32);
                if raw_minus.score < ungapped_cutoff { continue; }
                let ung = raw_minus;
                crate::diag_count!(COUNT_UNGAPPED_HITS);
                ungapped.push(UngappedRecord {
                    q_seed: q_off, s_seed: s_off,
                    q_start: ung.q_start, q_end: ung.q_end,
                    s_start: ung.s_start, s_end: ung.s_end,
                    score: ung.score, strand: Strand::Minus,
                });
            }
        }

        if scan_range[0] > end_range { break; }
    }
    let _ = step; // scan_step is encoded in the packed subject scan; kept for parity with collect_ungapped
}

/// Phase 2a inner loop: for each pre-sorted ungapped hit, run preliminary gapped alignment
/// and populate the interval-tree containment filter.  Returns preliminary HSPs (already
/// purged via purge_prelim_common_endpoints) in the caller's coordinate space.
///
/// `ungapped` must already be sorted (same order as run_gapped_phase lines 1192-1218).
/// `qa_len` / `sa_len` are the non-sentinel lengths of query / subject.
fn run_phase2a_inner(
    query: &[u8],
    query_rc: &[u8],
    subject_plus: &[u8],
    qa_len: u32,
    sa_len: u32,
    ungapped: &[UngappedRecord],
    params: &SearchParams,
    matrix: &ScoreMatrix,
    lut_word_length: u32,
    _chunk_offset: u32,
) -> Vec<PrelimHsp> {
    let mut accepted_plus  = BlastIntervalTree::new(0, qa_len + 1, 0, sa_len + 1);
    let mut accepted_minus = BlastIntervalTree::new(0, qa_len + 1, 0, sa_len + 1);
    let mut prelim_hsps: Vec<PrelimHsp> = Vec::new();
    let mut ws = AlignWorkspace::new();

    for ung in ungapped {
        let (cand_s_start, cand_s_end) = match ung.strand {
            Strand::Plus  => (ung.s_start, ung.s_end),
            Strand::Minus => (sa_len - ung.s_end, sa_len - ung.s_start),
        };
        let plus = ung.strand == Strand::Plus;
        let contained = {
            let tree = if plus { &accepted_plus } else { &accepted_minus };
            tree.contains(ung.q_start, ung.q_end, cand_s_start, cand_s_end,
                          ung.score, plus, MIN_DIAG_SEP)
        };
        if contained {
            continue;
        }

        let (q_seed_prelim, s_seed_prelim, q_seed_tb, s_seed_tb) = match ung.strand {
            Strand::Plus => {
                let (tb_q, tb_s) = if ung.s_end >= ung.s_seed + 8 {
                    (ung.q_seed + 3, ung.s_seed + 3)
                } else {
                    (ung.q_seed, ung.s_seed)
                };
                let offset_adj = 4 - (tb_s % 4);
                let mut pre_s = tb_s + offset_adj;
                let mut pre_q = tb_q + offset_adj;
                if pre_s >= sa_len || pre_q >= qa_len {
                    pre_s -= 4;
                    pre_q -= 4;
                }
                (pre_q, pre_s, tb_q, tb_s)
            }
            Strand::Minus => {
                let lw_u32 = lut_word_length as u32;
                // NCBI gapped_start = rc_q_off + 3 (DynProg path adds 3).
                // In FWD-genome: q_seed = q0 + ext_right, so q_seed + lw - 4 = rc_q_off + 3 - ext_right.
                // Using q_seed (not q0_raw) correctly accounts for the extension shift.
                //
                // The +3 shift applies iff the ungapped alignment extends ≥8 bases to the RIGHT
                // of the saved seed in FORWARD-subject space (NCBI blast_gapalign.c:4109
                // `s_end >= init_hsp->s_off + 8`).  Mapping the (FWD_query, RC_subject) ungapped
                // record to FWD-subject: saved seed = sa_len - ung.s_seed - lw, ungapped end =
                // sa_len - ung.s_start.  The old condition `ung.s_end >= ung.s_seed + 8` tested
                // the RC-subject end against the RC-subject anchor (wrong frame/direction) — see
                // the matching fix and rationale in run_gapped_phase (bug #38).
                let s_off_saved_fwd = sa_len - ung.s_seed - lw_u32;
                let s_end_fwd       = sa_len - ung.s_start;
                let (tb_q, tb_s) = if s_end_fwd >= s_off_saved_fwd + 8 {
                    (ung.q_seed + lw_u32 - 4, ung.s_seed + lw_u32 - 4)
                } else {
                    (ung.q_seed, ung.s_seed)
                };
                let fwd_te_seed = sa_len - 1 - tb_s;
                let offset_adj = 4 - (fwd_te_seed % 4);
                let pre_s = tb_s.saturating_sub(offset_adj);
                let pre_q = tb_q.saturating_sub(offset_adj);
                (pre_q, pre_s, tb_q, tb_s)
            }
        };

        // NCBI: mask_at_hash=TRUE for blastn/rmblastn — DUST only affects lookup table,
        // not gapped alignment. Use unmasked query for Phase 2a.
        let (pq_seq, ps_seq, pq_seed, ps_seed) = if ung.strand == Strand::Minus {
            (query_rc, subject_plus, qa_len - 1 - q_seed_prelim, sa_len - 1 - s_seed_prelim)
        } else {
            (query, subject_plus, q_seed_prelim, s_seed_prelim)
        };
        let phase2a_xdrop = params.xdrop_gap.min(ung.score);
        crate::diag_count!(COUNT_PRELIM_GAPPED);
        let prelim = gapped_extend_score_only(
            pq_seq, ps_seq, pq_seed, ps_seed,
            params.gap_open, params.gap_extend,
            phase2a_xdrop,
            matrix, &mut ws.dp_score,
            false,
        );

        let (prelim_score, q_start, q_end, s_start, s_end) = match prelim {
            None => { continue; }
            Some((sc, q0, q1, s0, s1)) => {
                if ung.strand == Strand::Minus {
                    (sc, qa_len - q1, qa_len - q0, sa_len - s1, sa_len - s0)
                } else {
                    (sc, q0, q1, s0, s1)
                }
            }
        };

        if prelim_score < params.min_raw_gapped_score {
            continue;
        }

        let (norm_s_start, norm_s_end) = match ung.strand {
            Strand::Plus => (s_start, s_end),
            Strand::Minus => (sa_len - s_end, sa_len - s_start),
        };

        {
            let tree = if ung.strand == Strand::Plus { &mut accepted_plus } else { &mut accepted_minus };
            let _ = tree.add(ITreeHsp { q_start, q_end, s_start: norm_s_start, s_end: norm_s_end,
                                        score: prelim_score }, plus);
        }

        prelim_hsps.push(PrelimHsp {
            q_seed: q_seed_tb,
            s_seed: s_seed_tb,
            q_start,
            q_end,
            s_start,
            s_end,
            score: prelim_score,
            strand: ung.strand,
        });
    }

    purge_prelim_common_endpoints(&mut prelim_hsps, sa_len);
    prelim_hsps
}

/// Phase 2b traceback + Phase 2c containment + endpoint purge + mask_level filter.
///
/// Mirrors Blast_TracebackFromHSPList.  `prelim_hsps` must be in global coordinates
/// (q_start, q_end, q_seed offset to the full genome; s coords in strand-native space).
/// `query` is the FULL genome sequence (with sentinels); `qa_len` is derived internally.
pub fn run_phase2b(
    query: &[u8],
    query_rc: &[u8],
    query_id: &str,
    subject_plus: &[u8],
    subject_id: &str,
    mut prelim_hsps: Vec<PrelimHsp>,
    params: &SearchParams,
    matrix: &ScoreMatrix,
    n_mask: &[u8],
) -> Vec<AlignResult> {
    let qa_len = (query.len() - 2) as u32;
    let s_len  = (subject_plus.len() - 2) as u32;
    let sa_len = s_len;

    // query_rc (reverse-complement of the FULL query) is computed once by the caller
    // and shared by reference across all per-subject Phase 2b tasks — recomputing it
    // here per task held one ~full-query copy per concurrent thread (see project_context
    // memory note on the chr22×longlib footprint).
    // CRandom bases at N/IUPAC positions are used only for seeding/ungapped extension.
    // For gapped alignment and improve_seed, NCBI scores ambiguous positions using the
    // original ambiguity code (14=N, or 4-13 for S/M/R/Y/K/W etc.).
    let subject_align_vec: Vec<u8> = if n_mask.is_empty() {
        subject_plus.to_vec()
    } else {
        let mut s = subject_plus.to_vec();
        for (i, &code) in n_mask.iter().enumerate() {
            if code != 0 { s[i + 1] = code; }
        }
        s
    };
    let subject_plus: &[u8] = &subject_align_vec;
    let subj_rc  = revcomp_blastna(subject_plus);
    let mut ws   = AlignWorkspace::new();

    // ScoreCompareHSPs sort: score DESC, FWD s_start ASC, FWD s_end DESC,
    // query.offset ASC, query.end DESC.  NCBI coordinate transform per strand:
    //   plus:  query.offset = FWD_qs,        query.end = FWD_qe
    //   minus: query.offset = qlen - FWD_qe, query.end = qlen - FWD_qs
    // Using NCBI-space values gives a valid total order for both pure and mixed-strand lists.
    prelim_hsps.sort_unstable_by(|a, b| {
        let fwd_sa = match a.strand { Strand::Plus => a.s_start, Strand::Minus => s_len - a.s_end };
        let fwd_sb = match b.strand { Strand::Plus => b.s_start, Strand::Minus => s_len - b.s_end };
        let fwd_ea = match a.strand { Strand::Plus => a.s_end, Strand::Minus => s_len - a.s_start };
        let fwd_eb = match b.strand { Strand::Plus => b.s_end, Strand::Minus => s_len - b.s_start };
        let q_off_a = match a.strand { Strand::Plus => a.q_start, Strand::Minus => qa_len - a.q_end };
        let q_off_b = match b.strand { Strand::Plus => b.q_start, Strand::Minus => qa_len - b.q_end };
        let q_end_a = match a.strand { Strand::Plus => a.q_end, Strand::Minus => qa_len - a.q_start };
        let q_end_b = match b.strand { Strand::Plus => b.q_end, Strand::Minus => qa_len - b.q_start };
        b.score.cmp(&a.score)
            .then_with(|| fwd_sa.cmp(&fwd_sb))
            .then_with(|| fwd_eb.cmp(&fwd_ea))
            .then_with(|| q_off_a.cmp(&q_off_b))
            .then_with(|| q_end_b.cmp(&q_end_a))
    });

    let mut tb_accepted_plus  = BlastIntervalTree::new(0, qa_len + 1, 0, s_len + 1);
    let mut tb_accepted_minus = BlastIntervalTree::new_rc(0, qa_len + 1, 0, s_len + 1);
    let mut results: Vec<AlignResult> = Vec::new();

    for phsp in &prelim_hsps {
        let (ps_start, ps_end) = match phsp.strand {
            Strand::Plus  => (phsp.s_start, phsp.s_end),
            Strand::Minus => (s_len - phsp.s_end, s_len - phsp.s_start),
        };
        let plus2 = phsp.strand == Strand::Plus;
        let prelim_contained = {
            let tree = if plus2 { &tb_accepted_plus } else { &tb_accepted_minus };
            let (cq_start, cq_end) = if plus2 {
                (phsp.q_start, phsp.q_end)
            } else {
                (qa_len - phsp.q_end, qa_len - phsp.q_start)
            };
            tree.contains(cq_start, cq_end, ps_start, ps_end,
                          phsp.score, plus2, MIN_DIAG_SEP)
        };
        if prelim_contained {
            continue;
        }

        let sa_for_seed = if phsp.strand == Strand::Minus {
            &subj_rc[1..subj_rc.len() - 1]
        } else {
            &subject_plus[1..subject_plus.len() - 1]
        };
        // NCBI: mask_at_hash=TRUE — unmasked query for improve_seed and Phase 2b traceback.
        let qa_for_seed = &query[1..query.len() - 1];

        let (new_q_seed, new_s_seed) = improve_seed(
            qa_for_seed, sa_for_seed,
            phsp.q_seed, phsp.s_seed,
            phsp.q_start, phsp.q_end,
            phsp.s_start, phsp.s_end,
            phsp.strand == Strand::Minus,
        );
        if crate::diag_enabled!("RMBLAST_DUMP_IMPROVE") {
            eprintln!("RUST_IMPROVE strand={:?} prelim q=[{},{}] s=[{},{}] seed=({},{}) -> ({},{})",
                phsp.strand, phsp.q_start, phsp.q_end, phsp.s_start, phsp.s_end,
                phsp.q_seed, phsp.s_seed, new_q_seed, new_s_seed);
        }

        let (tq_seq, ts_seq, tq_seed, ts_seed) = if phsp.strand == Strand::Minus {
            let tq = qa_len - 1 - new_q_seed;
            let ts = sa_len - 1 - new_s_seed;
            (&query_rc[..], subject_plus, tq, ts)
        } else {
            (query, subject_plus, new_q_seed, new_s_seed)
        };

        crate::diag_count!(COUNT_FINAL_GAPPED);
        let gapped = gapped_extend_bidirectional(
            tq_seq, ts_seq, tq_seed, ts_seed,
            params.gap_open, params.gap_extend,
            params.xdrop_gap_final,
            matrix, &mut ws,
        );
        if crate::diag_enabled!("RMBLAST_DUMP_IMPROVE") {
            if let Some((sc, q0, q1, s0, s1, _)) = &gapped {
                eprintln!("RUST_TB seed=({},{}) -> score={} q=[{},{}] s=[{},{}]",
                    tq_seed, ts_seed, sc, q0, q1, s0, s1);
            }
        }
        let (mut score, q_start, q_end, s_start, s_end, mut edit_script, q_bases, s_bases) = match gapped {
            None => { continue; }
            Some((sc, q0, q1, s0, s1, es)) => {
                if phsp.strand == Strand::Minus {
                    // Extract from (query_rc, subject_plus) using the unreversed edit script.
                    // revcomp_blastna then flips both to (query_fwd, subj_rc) orientation.
                    // This preserves NCBI's exact gap placement within homopolymer runs —
                    // reversing the script after extraction shifts gaps in runs by one position.
                    let (qb, sb) = extract_aligned(
                        &query_rc[1..query_rc.len() - 1],
                        &subject_plus[1..subject_plus.len() - 1],
                        q0 as usize,
                        s0 as usize,
                        &es,
                        n_mask,
                    );
                    (sc, qa_len - q1, qa_len - q0, sa_len - s1, sa_len - s0, es,
                     revcomp_blastna(&qb), revcomp_blastna(&sb))
                } else {
                    let (qb, sb) = extract_aligned(
                        &query[1..query.len() - 1],
                        &subject_plus[1..subject_plus.len() - 1],
                        q0 as usize,
                        s0 as usize,
                        &es,
                        n_mask,
                    );
                    (sc, q0, q1, s0, s1, es, qb, sb)
                }
            }
        };
        if phsp.strand == Strand::Minus {
            edit_script.reverse();
        }

        if params.complexity_adjust {
            match apply_complexity_adjust(
                score, query, q_start, &edit_script, matrix,
                params.min_raw_gapped_score,
            ) {
                Some(adj) => score = adj,
                None => { continue; }
            }
        }
        // Note: NCBI does NOT re-apply min_raw_gapped_score after Phase 2b traceback.
        // The threshold is only applied in Phase 2a. Hits that pass Phase 2a but whose
        // final Phase 2b score drops below the threshold are still output.

        let (report_q_start, report_q_end, report_s_start, report_s_end) = match phsp.strand {
            Strand::Plus  => (q_start, q_end, s_start, s_end),
            Strand::Minus => (q_start, q_end, s_len - s_end, s_len - s_start),
        };

        crate::diag_count!(COUNT_FINAL_HITS);

        {
            let tree = if phsp.strand == Strand::Plus { &mut tb_accepted_plus } else { &mut tb_accepted_minus };
            let (tree_qs, tree_qe) = if phsp.strand == Strand::Plus {
                (report_q_start, report_q_end)
            } else {
                (qa_len - report_q_end, qa_len - report_q_start)
            };
            tree.add(ITreeHsp { q_start: tree_qs, q_end: tree_qe,
                                s_start: report_s_start, s_end: report_s_end, score }, true);
        }

        let q_iupac = blastna_to_iupac_aligned(&q_bases);
        let s_iupac = blastna_to_iupac_aligned(&s_bases);
        let stats = compute_align_stats(&q_iupac, &s_iupac, false);


        results.push(AlignResult {
            hsp: Hsp {
                score,
                q_start: report_q_start,
                q_end: report_q_end,
                q_len: qa_len,
                s_start: report_s_start,
                s_end: report_s_end,
                s_len,
                strand: phsp.strand,
                edit_script,
                q_seq: q_bases,
                s_seq: s_bases,
            },
            query_id: query_id.to_string(),
            subject_id: subject_id.to_string(),
            stats,
        });
    }

    // ── NCBI post-processing: purge common endpoints with cut + reevaluate ────
    // Mirrors blast_traceback.c: Blast_HSPListPurgeHSPsWithCommonEndpoints(FALSE)
    // + Blast_HSPReevaluateWithAmbiguitiesGapped loop + PurgeHSPs(TRUE).
    {
        // NCBI blast_traceback.c:658/688 — purge(FALSE) cut pass, re-trace the cut
        // remainders (which sit at arr[extra_start..]) in array order, PurgeNull,
        // then purge(TRUE) delete pass.
        let mut arr: Vec<Option<AlignResult>> =
            std::mem::take(&mut results).into_iter().map(Some).collect();
        let extra_start = purge_hsps_with_common_endpoints(&mut arr, false);

        // Reevaluate each cut remainder and keep survivors.
        let q_core   = &query[1..query.len() - 1];
        let src_core = &subject_plus[1..subject_plus.len() - 1];
        let rc_core  = &subj_rc[1..subj_rc.len() - 1];
        // Re-extraction for minus strand uses rc_core (reversed), so n_mask must be
        // reversed and each IUPAC code complemented to match the RC strand orientation.
        let n_mask_rc_gapped: Vec<u8> = if n_mask.is_empty() { Vec::new() } else {
            n_mask.iter().rev().map(|&c| if c == 0 { 0 } else { BLASTNA_COMPLEMENT[c as usize] }).collect()
        };

        for slot in arr[extra_start..].iter_mut() {
            let mut r = match slot.take() { Some(r) => r, None => continue };
            let is_minus = r.hsp.strand == Strand::Minus;
            // Re-extract aligned sequences in the same orientation as the original call:
            //   Plus:  query[1..] from q_start, subject_plus[1..] from s_start
            //   Minus: query[1..] from q_start (= qa_len - DP.q1), subj_rc[1..] from sa_len - s_end
            let a_start = r.hsp.q_start as usize;
            let b_start = if is_minus {
                (sa_len - r.hsp.s_end) as usize
            } else {
                r.hsp.s_start as usize
            };
            let b_core = if is_minus { rc_core } else { src_core };
            let b_mask_r = if is_minus { &n_mask_rc_gapped[..] } else { n_mask };
            let (q_seq, s_seq) = extract_aligned(q_core, b_core, a_start, b_start, &r.hsp.edit_script, b_mask_r);

            let delete = reevaluate_gapped(
                &mut r.hsp,
                &q_seq,
                &s_seq,
                matrix,
                params.gap_open,
                params.gap_extend,
                // NCBI Blast_HSPReevaluateWithAmbiguitiesGapped uses cutoff_score
                // (= cutoffs[context].cutoff_score = min_raw_gapped_score for rmblastn)
                // both for the run-restart logic and the final keep/delete decision
                // (s_UpdateReevaluatedHSP: keep iff score >= cutoff_score).  Cut
                // secondaries that re-score below this are deleted, not kept.
                params.min_raw_gapped_score,
                is_minus,
                q_core,
                a_start,
                b_core,
                b_start,
            );
            if delete { continue; }

            // Rebuild q_seq/s_seq/stats for the trimmed region.
            let a_start2 = r.hsp.q_start as usize;
            let b_start2 = if is_minus {
                (sa_len - r.hsp.s_end) as usize
            } else {
                r.hsp.s_start as usize
            };
            let (q2, s2) = extract_aligned(q_core, b_core, a_start2, b_start2, &r.hsp.edit_script, b_mask_r);
            let q_iupac2 = blastna_to_iupac_aligned(&q2);
            let s_iupac2 = blastna_to_iupac_aligned(&s2);
            r.stats = compute_align_stats(&q_iupac2, &s_iupac2, false);
            r.hsp.q_seq = q2;
            r.hsp.s_seq = s2;

            *slot = Some(r);
        }

        // PurgeNull: compact survivors + surviving remainders (preserves order).
        results = arr.into_iter().flatten().collect();

        // Second pass (purge=TRUE): delete-all common-endpoint duplicates.
        let mut arr2: Vec<Option<AlignResult>> =
            std::mem::take(&mut results).into_iter().map(Some).collect();
        purge_hsps_with_common_endpoints(&mut arr2, true);
        results = arr2.into_iter().flatten().collect();
    }

    // Phase 2c: final containment filter.  Sort by ScoreCompareHSPs (NCBI coordinate space):
    // score DESC, s_start ASC, s_end DESC, query.offset ASC, query.end DESC.
    // For minus strand: query.offset = qlen - FWD_qe (ASC → FWD_qe DESC),
    //                   query.end   = qlen - FWD_qs (DESC → FWD_qs ASC).
    results.sort_by(|a, b| {
        let q_off_a = match a.hsp.strand { Strand::Plus => a.hsp.q_start, Strand::Minus => qa_len - a.hsp.q_end };
        let q_off_b = match b.hsp.strand { Strand::Plus => b.hsp.q_start, Strand::Minus => qa_len - b.hsp.q_end };
        let q_end_a = match a.hsp.strand { Strand::Plus => a.hsp.q_end, Strand::Minus => qa_len - a.hsp.q_start };
        let q_end_b = match b.hsp.strand { Strand::Plus => b.hsp.q_end, Strand::Minus => qa_len - b.hsp.q_start };
        b.hsp.score.cmp(&a.hsp.score)
            .then_with(|| a.hsp.s_start.cmp(&b.hsp.s_start))
            .then_with(|| b.hsp.s_end.cmp(&a.hsp.s_end))
            .then_with(|| q_off_a.cmp(&q_off_b))
            .then_with(|| q_end_b.cmp(&q_end_a))
    });
    let mut final_plus  = BlastIntervalTree::new(0, qa_len + 1, 0, s_len + 1);
    let mut final_minus = BlastIntervalTree::new(0, qa_len + 1, 0, s_len + 1);
    results.retain(|r| {
        let plus3 = r.hsp.strand == Strand::Plus;
        let contained = {
            let tree = if plus3 { &final_plus } else { &final_minus };
            tree.contains(r.hsp.q_start, r.hsp.q_end, r.hsp.s_start, r.hsp.s_end,
                          r.hsp.score, plus3, MIN_DIAG_SEP)
        };
        if !contained {
            let tree = if plus3 { &mut final_plus } else { &mut final_minus };
            tree.add_simple(ITreeHsp { q_start: r.hsp.q_start, q_end: r.hsp.q_end,
                                       s_start: r.hsp.s_start, s_end: r.hsp.s_end,
                                       score: r.hsp.score });
        }
        !contained
    });

    results
}

/// Phase 1 (ungapped) + Phase 2a (preliminary gapped) for one chunk × subject pair.
///
/// Returns preliminary HSPs in **global** coordinates: q_start, q_end, q_seed are
/// offset by `chunk_offset`; s coords are in strand-native subject space (unchanged).
/// Call [`merge_chunk_prelims`] to combine adjacent-chunk results, then [`run_phase2b`]
/// with the full query for the final traceback.
/// Phase 2a: seeding (word extension + ungapped extension) + preliminary gapped alignment.
///
/// NCBI applies DUST masking as soft masking: the masked query is used for seeding
/// (LUT + word extension check + ungapped extension scoring) AND for preliminary gapped
/// Phase 2a: preliminary gapped alignment for multi-chunk searches.
///
/// `query`        — original unmasked chunk; used for Phase 2a gapped alignment (mask_at_hash=TRUE).
/// `masked_query` — DUST-masked chunk; used for seeding (LUT lookup, word extension, ungapped scoring).
///                  When dust=false, pass the same slice as query.
pub fn search_phase2a(
    ql: &QueryLookup,
    query: &[u8],
    masked_query: &[u8],
    subject_plus: &[u8],
    _subject_id: &str,
    params: &SearchParams,
    matrix: &ScoreMatrix,
    chunk_offset: u32,
    prepared: Option<&PreparedSubject>,
) -> Vec<PrelimHsp> {
    let q_len = query.len().saturating_sub(2) as u32;
    let s_len = subject_plus.len().saturating_sub(2) as u32;

    if q_len < params.word_size as u32 || s_len < params.word_size as u32 {
        return Vec::new();
    }

    // Subject RC + NCBI2NA-packed strands for the scan: borrow from the shared
    // per-subject cache when supplied, else compute locally (identical bytes).
    let local_prep;
    let prep: &PreparedSubject = match prepared {
        Some(p) => p,
        None => {
            // Phase 2a never touches the align views, so no n_mask is needed here.
            local_prep = prepare_subject_strands(subject_plus, &[]);
            &local_prep
        }
    };
    let (subj_rc, packed_plus, packed_minus): (&[u8], &[u8], &[u8]) =
        (&prep.rc, &prep.packed_plus, &prep.packed_minus);

    let mut ungapped = Vec::new();
    if ql.c1_start > 0 {
        collect_ungapped_combined(
            query, masked_query, q_len, subject_plus, &packed_plus[3..], subj_rc, s_len,
            ql.c1_start, ql.fwd_n, &ql.lookup, params, matrix, &mut DiscardSeeds, &mut ungapped,
        );
    } else {
        collect_ungapped(
            query, masked_query, q_len, subject_plus, &packed_plus[3..], s_len,
            Strand::Plus, &ql.lookup, params, matrix, &mut DiscardSeeds, &mut ungapped,
        );
        collect_ungapped(
            query, masked_query, q_len, subj_rc, &packed_minus[3..], s_len,
            Strand::Minus, &ql.lookup, params, matrix, &mut DiscardSeeds, &mut ungapped,
        );
    }

    // Same sort as run_gapped_phase (Blast_InitHitListSortByScore).
    ungapped.sort_unstable_by(|a, b| {
        let a_ss = match a.strand { Strand::Plus => a.s_start, Strand::Minus => s_len - a.s_end };
        let b_ss = match b.strand { Strand::Plus => b.s_start, Strand::Minus => s_len - b.s_end };
        let a_qk = match a.strand {
            Strand::Plus  => a.q_start as u64,
            Strand::Minus => (q_len as u64 + 1) + (q_len as u64).saturating_sub(a.q_end as u64),
        };
        let b_qk = match b.strand {
            Strand::Plus  => b.q_start as u64,
            Strand::Minus => (q_len as u64 + 1) + (q_len as u64).saturating_sub(b.q_end as u64),
        };
        let a_sk = match a.strand { Strand::Plus => a.s_seed, Strand::Minus => s_len - a.s_seed };
        let b_sk = match b.strand { Strand::Plus => b.s_seed, Strand::Minus => s_len - b.s_seed };
        b.score.cmp(&a.score)
            .then_with(|| a_ss.cmp(&b_ss))
            .then_with(|| {
                let la = a.s_end - a.s_start;
                let lb = b.s_end - b.s_start;
                lb.cmp(&la)
            })
            .then_with(|| a_qk.cmp(&b_qk))
            .then_with(|| a_sk.cmp(&b_sk))
    });

    let query_rc = revcomp_blastna(query);
    let qa_len   = q_len;
    let sa_len   = s_len;

    let mut prelims = run_phase2a_inner(
        query, &query_rc, subject_plus, qa_len, sa_len, &ungapped,
        params, matrix, ql.lut_word_length(), chunk_offset,
    );

    if chunk_offset > 0 {
        for p in &mut prelims {
            p.q_start += chunk_offset;
            p.q_end   += chunk_offset;
            p.q_seed  += chunk_offset;
        }
    }

    prelims
}

/// Phase 2: preliminary gapped + traceback.
///
/// Phase 2a (BLAST_GetGappedScore): sorts ungapped hits by score descending, runs
/// preliminary gapped alignment with xdrop_gap (30) to populate the interval tree.
///
/// Phase 2b (Blast_TracebackFromHSPList): for each accepted preliminary HSP, finds
/// a better seed via BlastGetStartForGappedAlignmentNucl, then runs full gapped
/// alignment with traceback using xdrop_gap_final (100).
fn run_gapped_phase<U: UngapOut>(
    query: &[u8],
    query_id: &str,
    subject_plus: &[u8],
    _subj_rc: &[u8],
    subject_id: &str,
    q_len: u32,
    s_len: u32,
    mut ungapped: Vec<UngappedRecord>,
    params: &SearchParams,
    matrix: &ScoreMatrix,
    ungap_out: &mut U,
    lut_word_length: u32,
    n_mask: &[u8],
    prep: &PreparedSubject,
) -> Vec<AlignResult> {
    // N/IUPAC handling: NCBI uses CRandom (2-bit scan) for Phase 2a (score-only filter)
    // and the original ambiguity code for Phase 2b (full traceback + improve_seed + reevaluate).
    // Phase 2a uses the original subject_plus/subj_rc parameters (CRandom at ambig positions).
    // Phase 2b uses the align views below (original BLASTNA code at ambig positions),
    // borrowed from the per-subject cache (built with this same n_mask).
    let subject_align: &[u8] = &prep.align;
    let subj_rc_align: &[u8] = &prep.align_rc;

    // Mirrors Blast_InitHitListSortByScore / score_compare_match in blast_extend.c.
    // NCBI sorts by combined-query q_start, which for plus strand = FWD q_start (< q_len)
    // and for minus strand = c1_offset + (q_len - FWD_q_end) (>> q_len).  This naturally
    // places all plus-strand hits before all minus-strand hits in the q_start tiebreaker.
    // We replicate that by mapping: plus → q_start, minus → (q_len+1) + (q_len - q_end).
    // NCBI also stores s_start in FWD-subject space (from the FWD scan), so we FWD-normalize.
    ungapped.sort_unstable_by(|a, b| {
        let a_ss = match a.strand { Strand::Plus => a.s_start, Strand::Minus => s_len - a.s_end };
        let b_ss = match b.strand { Strand::Plus => b.s_start, Strand::Minus => s_len - b.s_end };
        let a_qk = match a.strand { Strand::Plus => a.q_start as u64,
                                    Strand::Minus => (q_len as u64 + 1) + (q_len as u64).saturating_sub(a.q_end as u64) };
        let b_qk = match b.strand { Strand::Plus => b.q_start as u64,
                                    Strand::Minus => (q_len as u64 + 1) + (q_len as u64).saturating_sub(b.q_end as u64) };
        // s_seed ASC: word-hit position in FWD-subject coords — matches NCBI's effective
        // insertion order (left-to-right scan) for hits with otherwise equal sort keys.
        let a_sk = match a.strand { Strand::Plus => a.s_seed, Strand::Minus => s_len - a.s_seed };
        let b_sk = match b.strand { Strand::Plus => b.s_seed, Strand::Minus => s_len - b.s_seed };
        b.score.cmp(&a.score)
            .then_with(|| a_ss.cmp(&b_ss))
            .then_with(|| {
                let la = a.s_end - a.s_start;
                let lb = b.s_end - b.s_start;
                lb.cmp(&la)
            })
            .then_with(|| a_qk.cmp(&b_qk))
            .then_with(|| a_sk.cmp(&b_sk))
    });

    // Emit ungapped hits (sorted, FWD-normalized) for collection and optional debug dump.
    for ung in &ungapped {
        let (ss, se) = match ung.strand {
            Strand::Plus  => (ung.s_start, ung.s_end),
            Strand::Minus => (s_len - ung.s_end, s_len - ung.s_start),
        };
        ungap_out.push_ungap(UngappedHit {
            q_start: ung.q_start, q_end: ung.q_end,
            s_start: ss, s_end: se,
            score: ung.score, strand: ung.strand,
        });
    }

    // Interval trees for Phase 2a containment (one per strand, mirrors NCBI BlastIntervalTree).
    // Tree range = [0, query_non_sentinel_length + 1] matching NCBI's Blast_IntervalTreeInit.
    let q_tree_end = (query.len() - 2) as u32 + 1;
    let mut accepted_plus  = BlastIntervalTree::new(0, q_tree_end, 0, s_len + 1);
    let mut accepted_minus = BlastIntervalTree::new(0, q_tree_end, 0, s_len + 1);
    let mut prelim_hsps: Vec<PrelimHsp> = Vec::new();
    let mut ws = AlignWorkspace::new();

    // For minus-strand gapped alignment, use NCBI orientation (RC-query, FWD-subject) so the
    // seed always lands in the LEFT (reverse) extension, matching BLAST_GappedAlignmentWithTraceback.
    // Results are converted back to FWD-genomic + RC-subject before storage.
    // NCBI: mask_at_hash=TRUE — masked query only affects lookup table, not gapped alignment.
    let query_rc        = revcomp_blastna(query);
    let qa_len = (query.len() - 2) as u32;   // non-sentinel bases in query
    let sa_len = (subject_plus.len() - 2) as u32; // non-sentinel bases in subject

    // Phase 2a: preliminary gapped alignment with xdrop_gap (mirrors BLAST_GetGappedScore).
    for ung in &ungapped {
        // Normalize candidate s coords to plus-strand for the interval tree query.
        let (cand_s_start, cand_s_end) = match ung.strand {
            Strand::Plus  => (ung.s_start, ung.s_end),
            Strand::Minus => (s_len - ung.s_end, s_len - ung.s_start),
        };
        let plus = ung.strand == Strand::Plus;
        let contained = {
            let tree = if plus { &accepted_plus } else { &accepted_minus };
            tree.contains(ung.q_start, ung.q_end, cand_s_start, cand_s_end,
                          ung.score, plus, MIN_DIAG_SEP)
        };
        if contained {
            continue;
        }

        // NCBI shifts the seed +3 before Phase 2a if the ungapped extension extends ≥8
        // bases past the word-start (blast_gapalign.c lines 4041-4044).  It then rounds
        // the Phase 2a pivot to the next 4-base subject boundary inside
        // s_BlastDynProgNtGappedAlignment (offset_adjustment = 4 - (s_off % 4)).
        // This rounding changes the prelim alignment boundaries, which in turn affects the
        // search range used by BlastGetStartForGappedAlignmentNucl (improve_seed) in Phase
        // 2b.  We must replicate the rounding so our improve_seed search covers the same
        // region, even though 1-byte-per-base alignment doesn't require it algorithmically.
        //
        // Phase 2a uses (q_seed_prelim, s_seed_prelim) — the rounded pivot.
        // Phase 2b uses (q_seed_tb, s_seed_tb) — the unrounded seed (stored as gapped_start).
        //
        // For plus strand: both +3 shift and 4-base rounding apply to Phase 2a.
        //   offset_adj = 4 - ((s_seed + 3) % 4)  (always 1..=4)
        //   prelim seed = (q_seed + 3 + offset_adj, s_seed + 3 + offset_adj)
        //   tb seed     = (q_seed + 3, s_seed + 3)
        //
        // For minus strand: NCBI context-1 FWD-TE shifts correspond to RC-TE shifts in
        //   the opposite direction.
        //   NCBI Phase 2b seed (context-1): (q_off+3, s_off+3) → Rust FWD: (q_seed+wm4, s_seed+wm4).
        //   NCBI Phase 2a rounding in FWD-TE: offset_adj = 4 - ((sa_len-1-s_seed_tb) % 4).
        //   In RC-TE this subtracts offset_adj: prelim_s = s_seed_tb - offset_adj.
        let (q_seed_prelim, s_seed_prelim, q_seed_tb, s_seed_tb) = match ung.strand {
            Strand::Plus => {
                let (tb_q, tb_s) = if ung.s_end >= ung.s_seed + 8 {
                    (ung.q_seed + 3, ung.s_seed + 3)
                } else {
                    (ung.q_seed, ung.s_seed)
                };
                // Round Phase 2a pivot to next 4-base subject boundary (offset_adjustment).
                let offset_adj = 4 - (tb_s % 4);
                let mut pre_s = tb_s + offset_adj;
                let mut pre_q = tb_q + offset_adj;
                // Bound check mirrors NCBI: if pivot overshoots, pull back by 4.
                if pre_s >= sa_len || pre_q >= qa_len {
                    pre_s -= 4;
                    pre_q -= 4;
                }
                (pre_q, pre_s, tb_q, tb_s)
            }
            Strand::Minus => {
                let lw_u32 = lut_word_length as u32;
                // NCBI's +3 gapped_start shift (blast_gapalign.c:4109) applies iff the ungapped
                // alignment extends ≥8 bases to the RIGHT of the saved seed in FORWARD-subject
                // space: `s_end >= init_hsp->s_off + 8`, where init_hsp->s_off is the saved seed
                // (FWD-subject) and s_end is the ungapped end (FWD-subject).  Rust stores the
                // ungapped record in (FWD_query, RC_subject) coords, so map to FWD-subject:
                //   saved seed FWD_s = sa_len - ung.s_seed - lw   (== NCBI init_hsp->s_off)
                //   ungapped end FWD_s = sa_len - ung.s_start
                // The previous condition `ung.s_end >= ung.s_seed + 8` tested the RC-subject end
                // against the RC-subject anchor — the wrong frame AND the wrong direction.  It
                // happened to agree with NCBI for most hits but gave the wrong answer when the
                // ungapped extension to the right (FWD_s) of the seed landed within 8 bases of
                // the boundary, dropping the +3 shift → traceback gapped_start 3bp off → 1bp
                // gap-placement shift in homopolymer runs whenever improve_seed early-returns
                // (bug #38).  This is the minus-strand analog of the plus branch's natural
                // FWD-subject condition above.
                let s_off_saved_fwd = sa_len - ung.s_seed - lw_u32; // saved seed, FWD-subject
                let s_end_fwd       = sa_len - ung.s_start;         // ungapped end, FWD-subject
                let extends = s_end_fwd >= s_off_saved_fwd + 8;
                let (tb_q, tb_s) = if extends {
                    (ung.q_seed + lw_u32 - 4, ung.s_seed + lw_u32 - 4)
                } else {
                    (ung.q_seed, ung.s_seed)
                };
                // Phase 2a seed: mirror NCBI's s_BlastDynProgNtGappedAlignment exactly, in
                // NCBI's own (RC_query, FWD_subject) frame.  NCBI takes the ungapped seed
                // offsets (with the +3 gapped_start shift above), then
                // `offset_adjustment = 4 - (s_off % 4)` shifts the start UP to the next 4-base
                // subject boundary (so the packed-nucleotide left extension starts on a byte
                // boundary), with a guard that subtracts 4 when the shifted start would run past
                // the END of either sequence.
                //
                // Computing the pivot in this FWD frame (rather than the mirror RC_subject
                // frame) is essential: NCBI's guard trims from the HIGH end, so the pivot can
                // never fall below 0.  The previous RC-frame formula applied the guard in the
                // wrong direction and underflowed for seeds at the subject boundary
                // (ung.s_seed≈0), producing a negative seed → out-of-bounds panic (mirslib,
                // word_size 6).  For all non-boundary seeds this yields the identical pivot.
                let shift = if extends { 3 } else { 0 };
                let q_off_n = (qa_len - ung.q_seed - lw_u32) + shift; // RC_q word-left (+3)
                let s_off_n = (sa_len - ung.s_seed - lw_u32) + shift; // FWD_s word-left (+3)
                let oa = 4 - (s_off_n % 4); // offset_adjustment ∈ [1,4]
                let mut q_len_piv = q_off_n + oa;
                let mut s_len_piv = s_off_n + oa;
                if q_len_piv > qa_len || s_len_piv > sa_len {
                    q_len_piv -= 4;
                    s_len_piv -= 4;
                }
                // Since bug #42, gapped_extend_score_only uses the pivot as NCBI's
                // q_length/s_length directly — the EXCLUSIVE left bound / FIRST base of the
                // right extension (left = [0, pivot)).  So the prelim pivot in NCBI's
                // (RC_q, FWD_s) frame must equal q_len_piv / s_len_piv.  The caller maps
                // back via `qa_len-1-pre_q` / `sa_len-1-pre_s`, so store the FWD-frame
                // complement that yields pivot == q_len_piv: `qa_len-1-q_len_piv`.  (The
                // plus branch already passes q_length directly; before #42 this used
                // `qa_len-q_len_piv` to give pivot q_length-1, matching the old pivot-in-left
                // split — that one-base offset became a latent minus-only bug at #42, only
                // visible at the prelim-filter boundary, e.g. cc L2d_3end 181-vs-179.)
                let pre_q = qa_len - 1 - q_len_piv;
                let pre_s = sa_len - 1 - s_len_piv;
                (pre_q, pre_s, tb_q, tb_s)
            }
        };

        // Preliminary gapped alignment anchored at seed, with smaller xdrop_gap.
        // Minus strand: use NCBI orientation (RC-query, FWD-TE) so seed lands in LEFT.
        // NCBI: mask_at_hash=TRUE — DUST only affects lookup table, not gapped alignment.
        let (pq_seq, ps_seq, pq_seed, ps_seed) = if ung.strand == Strand::Minus {
            (&query_rc[..], subject_plus, qa_len - 1 - q_seed_prelim, sa_len - 1 - s_seed_prelim)
        } else {
            (query, subject_plus, q_seed_prelim, s_seed_prelim)
        };
        // NCBI caps Phase 2a xdrop at min(xdrop_gap, ungapped_score):
        // blast_gapalign.c s_BlastDynProgNtGappedAlignment lines 2970-2972.
        let phase2a_xdrop = params.xdrop_gap.min(ung.score);
        crate::diag_count!(COUNT_PRELIM_GAPPED);
        let prelim = gapped_extend_score_only(
            pq_seq, ps_seq, pq_seed, ps_seed,
            params.gap_open, params.gap_extend,
            phase2a_xdrop,
            matrix, &mut ws.dp_score,
            false,
        );

        let (prelim_score, q_start, q_end, s_start, s_end) = match prelim {
            None => {
                continue;
            }
            Some((sc, q0, q1, s0, s1)) => {
                if ung.strand == Strand::Minus {
                    // RC-query + FWD-TE → FWD-genomic + RC-TE
                    (sc, qa_len - q1, qa_len - q0, sa_len - s1, sa_len - s0)
                } else {
                    (sc, q0, q1, s0, s1)
                }
            }
        };

        if prelim_score < params.min_raw_gapped_score {
            continue;
        }

        // Normalize subject coordinates to plus-strand for the interval tree.
        let (norm_s_start, norm_s_end) = match ung.strand {
            Strand::Plus => (s_start, s_end),
            Strand::Minus => (s_len - s_end, s_len - s_start),
        };

        // Add to interval tree using PRELIMINARY (compact) boundaries.
        {
            let tree = if ung.strand == Strand::Plus { &mut accepted_plus } else { &mut accepted_minus };
            tree.add(ITreeHsp { q_start, q_end, s_start: norm_s_start, s_end: norm_s_end, score: prelim_score }, plus);
        }

        prelim_hsps.push(PrelimHsp {
            q_seed: q_seed_tb,
            s_seed: s_seed_tb,
            q_start,
            q_end,
            s_start,  // in strand's subject space
            s_end,
            score: prelim_score,
            strand: ung.strand,
        });
    }

    // Mirrors Blast_HSPListPurgeHSPsWithCommonEndpoints (called in blast_engine.c between
    // Phase 2a and Phase 2b). Removes prelim HSPs that share the same start or end
    // coordinates, keeping only the best-scoring one in each group.
    purge_prelim_common_endpoints(&mut prelim_hsps, s_len);

    // Phase 2b: traceback with xdrop_gap_final (mirrors Blast_TracebackFromHSPList).
    // NCBI sorts the preliminary HSP list by ScoreCompareHSPs: score DESC, then
    // subject.offset ASC (FWD s_start), subject.end DESC, query.offset ASC, query.end DESC.
    // NCBI coordinate transform per strand:
    //   plus:  query.offset = FWD_qs,        query.end = FWD_qe
    //   minus: query.offset = qlen - FWD_qe, query.end = qlen - FWD_qs
    prelim_hsps.sort_unstable_by(|a, b| {
        let fwd_sa = match a.strand { Strand::Plus => a.s_start, Strand::Minus => s_len - a.s_end };
        let fwd_sb = match b.strand { Strand::Plus => b.s_start, Strand::Minus => s_len - b.s_end };
        let fwd_ea = match a.strand { Strand::Plus => a.s_end, Strand::Minus => s_len - a.s_start };
        let fwd_eb = match b.strand { Strand::Plus => b.s_end, Strand::Minus => s_len - b.s_start };
        let q_off_a = match a.strand { Strand::Plus => a.q_start, Strand::Minus => qa_len - a.q_end };
        let q_off_b = match b.strand { Strand::Plus => b.q_start, Strand::Minus => qa_len - b.q_end };
        let q_end_a = match a.strand { Strand::Plus => a.q_end, Strand::Minus => qa_len - a.q_start };
        let q_end_b = match b.strand { Strand::Plus => b.q_end, Strand::Minus => qa_len - b.q_start };
        b.score.cmp(&a.score)
            .then_with(|| fwd_sa.cmp(&fwd_sb))
            .then_with(|| fwd_eb.cmp(&fwd_ea))
            .then_with(|| q_off_a.cmp(&q_off_b))
            .then_with(|| q_end_b.cmp(&q_end_a))
    });

    // Interval trees for Phase 2b containment (mirrors NCBI BlastIntervalTree in
    // Blast_TracebackFromHSPList — same structure as Phase 2a, separate per strand).
    let mut tb_accepted_plus  = BlastIntervalTree::new(0, qa_len + 1, 0, s_len + 1);
    let mut tb_accepted_minus = BlastIntervalTree::new_rc(0, qa_len + 1, 0, s_len + 1);
    let mut results: Vec<AlignResult> = Vec::new();

    for phsp in &prelim_hsps {
        // Containment check: is the preliminary HSP enclosed by an already-accepted
        // traceback HSP? Mirrors BlastIntervalTreeContainsHSP in Blast_TracebackFromHSPList.
        let (ps_start, ps_end) = match phsp.strand {
            Strand::Plus  => (phsp.s_start, phsp.s_end),
            Strand::Minus => (s_len - phsp.s_end, s_len - phsp.s_start),
        };
        let plus2 = phsp.strand == Strand::Plus;
        let prelim_contained = {
            let tree = if plus2 { &tb_accepted_plus } else { &tb_accepted_minus };
            let (cq_start, cq_end) = if plus2 {
                (phsp.q_start, phsp.q_end)
            } else {
                (qa_len - phsp.q_end, qa_len - phsp.q_start)
            };
            tree.contains(cq_start, cq_end, ps_start, ps_end,
                          phsp.score, plus2, MIN_DIAG_SEP)
        };
        if prelim_contained {
            continue;
        }

        // improve_seed still operates in FWD-genomic + RC-TE space (PrelimHsp coords).
        // NCBI: mask_at_hash=TRUE — unmasked query for improve_seed and Phase 2b traceback.
        let sa_for_seed = if phsp.strand == Strand::Minus {
            &subj_rc_align[1..subj_rc_align.len() - 1]
        } else {
            &subject_align[1..subject_align.len() - 1]
        };
        let qa_for_seed = &query[1..query.len() - 1];

        let (new_q_seed, new_s_seed) = improve_seed(
            qa_for_seed, sa_for_seed,
            phsp.q_seed, phsp.s_seed,
            phsp.q_start, phsp.q_end,
            phsp.s_start, phsp.s_end,
            phsp.strand == Strand::Minus,
        );
        if crate::diag_enabled!("RMBLAST_DUMP_IMPROVE") {
            eprintln!("RUST_IMPROVE strand={:?} prelim q=[{},{}] s=[{},{}] seed=({},{}) -> ({},{})",
                phsp.strand, phsp.q_start, phsp.q_end, phsp.s_start, phsp.s_end,
                phsp.q_seed, phsp.s_seed, new_q_seed, new_s_seed);
        }

        // Full traceback: minus strand uses NCBI orientation (RC-query, FWD-TE) so seed is in LEFT.
        let (tq_seq, ts_seq, tq_seed, ts_seed) = if phsp.strand == Strand::Minus {
            let tq = qa_len - 1 - new_q_seed;
            let ts = sa_len - 1 - new_s_seed;
            (&query_rc[..], &subject_align[..], tq, ts)
        } else {
            (query, &subject_align[..], new_q_seed, new_s_seed)
        };

        crate::diag_count!(COUNT_FINAL_GAPPED);
        let gapped = gapped_extend_bidirectional(
            tq_seq, ts_seq, tq_seed, ts_seed,
            params.gap_open, params.gap_extend,
            params.xdrop_gap_final,
            matrix, &mut ws,
        );
        if crate::diag_enabled!("RMBLAST_DUMP_IMPROVE") {
            if let Some((sc, q0, q1, s0, s1, _)) = &gapped {
                eprintln!("RUST_TB seed=({},{}) -> score={} q=[{},{}] s=[{},{}]",
                    tq_seed, ts_seed, sc, q0, q1, s0, s1);
            }
        }
        let (mut score, q_start, q_end, s_start, s_end, mut edit_script, q_bases, s_bases) = match gapped {
            None => continue,
            Some((sc, q0, q1, s0, s1, es)) => {
                if phsp.strand == Strand::Minus {
                    let (qb, sb) = extract_aligned(
                        &query_rc[1..query_rc.len() - 1],
                        &subject_align[1..subject_align.len() - 1],
                        q0 as usize,
                        s0 as usize,
                        &es,
                        n_mask,
                    );
                    let (fqs, fqe, fss, fse) = (qa_len - q1, qa_len - q0, sa_len - s1, sa_len - s0);
                    (sc, fqs, fqe, fss, fse, es,
                     revcomp_blastna(&qb), revcomp_blastna(&sb))
                } else {
                    let (qb, sb) = extract_aligned(
                        &query[1..query.len() - 1],
                        &subject_align[1..subject_align.len() - 1],
                        q0 as usize,
                        s0 as usize,
                        &es,
                        n_mask,
                    );
                    (sc, q0, q1, s0, s1, es, qb, sb)
                }
            }
        };
        if phsp.strand == Strand::Minus {
            edit_script.reverse();
        }

        // Score filter: complexity_adjust only. NCBI does not re-apply min_raw_gapped_score
        // after Phase 2b traceback — hits that passed Phase 2a are kept regardless of final score.
        if params.complexity_adjust {
            match apply_complexity_adjust(
                score, query, q_start, &edit_script, matrix,
                params.min_raw_gapped_score,
            ) {
                Some(adj) => score = adj,
                None => {
                    continue;
                }
            }
        }

        // Normalize subject coordinates to plus-strand for output.
        let (report_q_start, report_q_end, report_s_start, report_s_end) = match phsp.strand {
            Strand::Plus => (q_start, q_end, s_start, s_end),
            Strand::Minus => (q_start, q_end, s_len - s_end, s_len - s_start),
        };

        crate::diag_count!(COUNT_FINAL_HITS);

        // Add traceback result to the interval tree so subsequent prelim HSPs
        // can be checked against it.
        {
            let tree = if phsp.strand == Strand::Plus { &mut tb_accepted_plus } else { &mut tb_accepted_minus };
            let (tree_qs, tree_qe) = if phsp.strand == Strand::Plus {
                (report_q_start, report_q_end)
            } else {
                (qa_len - report_q_end, qa_len - report_q_start)
            };
            tree.add(ITreeHsp { q_start: tree_qs, q_end: tree_qe,
                                s_start: report_s_start, s_end: report_s_end, score }, true);
        }

        let q_iupac = blastna_to_iupac_aligned(&q_bases);
        let s_iupac = blastna_to_iupac_aligned(&s_bases);
        // subject is ancestral (repeat consensus), query is derived (genome) — matches NCBI
        let stats = compute_align_stats(&q_iupac, &s_iupac, false);

        results.push(AlignResult {
            hsp: Hsp {
                score,
                q_start: report_q_start,
                q_end: report_q_end,
                q_len,
                s_start: report_s_start,
                s_end: report_s_end,
                s_len,
                strand: phsp.strand,
                edit_script,
                q_seq: q_bases,
                s_seq: s_bases,
            },
            query_id: query_id.to_string(),
            subject_id: subject_id.to_string(),
            stats,
        });
    }

    // NCBI PurgeHSPsWithCommonEndpoints (blast_traceback.c lines 638-687):
    // First pass (purge=FALSE): cut secondaries and reevaluate via
    // Blast_HSPReevaluateWithAmbiguitiesGapped.  Second pass (purge=TRUE): delete.
    {
        // NCBI blast_traceback.c:658/688 — purge(FALSE) cut pass, re-trace the cut
        // remainders (arr[extra_start..]) in array order, PurgeNull, purge(TRUE).
        let mut arr: Vec<Option<AlignResult>> =
            std::mem::take(&mut results).into_iter().map(Some).collect();
        let extra_start = purge_hsps_with_common_endpoints(&mut arr, false);

        let q_core   = &query[1..query.len() - 1];
        let src_core = &subject_align[1..subject_align.len() - 1];
        let rc_core  = &subj_rc_align[1..subj_rc_align.len() - 1];
        let n_mask_rc_2b: &[u8] = &prep.n_mask_rc;

        for slot in arr[extra_start..].iter_mut() {
            let mut r = match slot.take() { Some(r) => r, None => continue };
            let is_minus = r.hsp.strand == Strand::Minus;
            let a_start = r.hsp.q_start as usize;
            let b_start = if is_minus {
                (sa_len - r.hsp.s_end) as usize
            } else {
                r.hsp.s_start as usize
            };
            let b_core = if is_minus { rc_core } else { src_core };
            let b_mask_2b = if is_minus { &n_mask_rc_2b[..] } else { n_mask };
            let (q_seq, s_seq) = extract_aligned(q_core, b_core, a_start, b_start, &r.hsp.edit_script, b_mask_2b);

            let delete = reevaluate_gapped(
                &mut r.hsp,
                &q_seq,
                &s_seq,
                matrix,
                params.gap_open,
                params.gap_extend,
                // NCBI Blast_HSPReevaluateWithAmbiguitiesGapped uses cutoff_score
                // (= cutoffs[context].cutoff_score = min_raw_gapped_score for rmblastn)
                // both for the run-restart logic and the final keep/delete decision
                // (s_UpdateReevaluatedHSP: keep iff score >= cutoff_score).  Cut
                // secondaries that re-score below this are deleted, not kept.
                params.min_raw_gapped_score,
                is_minus,
                q_core,
                a_start,
                b_core,
                b_start,
            );
            if delete { continue; }

            // Rebuild q_seq/s_seq/stats for the trimmed region.
            let a_start2 = r.hsp.q_start as usize;
            let b_start2 = if is_minus {
                (sa_len - r.hsp.s_end) as usize
            } else {
                r.hsp.s_start as usize
            };
            let (q2, s2) = extract_aligned(q_core, b_core, a_start2, b_start2, &r.hsp.edit_script, b_mask_2b);
            let q_iupac2 = blastna_to_iupac_aligned(&q2);
            let s_iupac2 = blastna_to_iupac_aligned(&s2);
            r.stats = compute_align_stats(&q_iupac2, &s_iupac2, false);
            r.hsp.q_seq = q2;
            r.hsp.s_seq = s2;

            *slot = Some(r);
        }

        results = arr.into_iter().flatten().collect();

        let mut arr2: Vec<Option<AlignResult>> =
            std::mem::take(&mut results).into_iter().map(Some).collect();
        purge_hsps_with_common_endpoints(&mut arr2, true);
        results = arr2.into_iter().flatten().collect();
    }

    // NCBI Phase 2c: third containment check (blast_traceback.c lines 678-692).
    // Sort traceback results by score descending, then for each result drop it if
    // its BOTH endpoints fall inside a higher-scoring (or equal) already-accepted
    // result on the same strand.  Mirrors BlastIntervalTreeContainsHSP, which
    // requires BOTH spatial containment AND diagonal proximity (MB_HSP_CLOSE):
    // at least one endpoint must be within MIN_DIAG_SEP diagonals of the accepted
    // hit.  Sub-alignments on different diagonals are NOT removed.
    // Sort traceback results by ScoreCompareHSPs (NCBI coordinate space):
    // score DESC, s_start ASC, s_end DESC, query.offset ASC, query.end DESC.
    // For minus strand: query.offset = qlen - FWD_qe (ASC → FWD_qe DESC),
    //                   query.end   = qlen - FWD_qs (DESC → FWD_qs ASC).
    results.sort_by(|a, b| {
        let q_off_a = match a.hsp.strand { Strand::Plus => a.hsp.q_start, Strand::Minus => qa_len - a.hsp.q_end };
        let q_off_b = match b.hsp.strand { Strand::Plus => b.hsp.q_start, Strand::Minus => qa_len - b.hsp.q_end };
        let q_end_a = match a.hsp.strand { Strand::Plus => a.hsp.q_end, Strand::Minus => qa_len - a.hsp.q_start };
        let q_end_b = match b.hsp.strand { Strand::Plus => b.hsp.q_end, Strand::Minus => qa_len - b.hsp.q_start };
        b.hsp.score.cmp(&a.hsp.score)
            .then_with(|| a.hsp.s_start.cmp(&b.hsp.s_start))
            .then_with(|| b.hsp.s_end.cmp(&a.hsp.s_end))
            .then_with(|| q_off_a.cmp(&q_off_b))
            .then_with(|| q_end_b.cmp(&q_end_a))
    });
    // Phase 2c containment: interval trees per strand (mirrors Phase 2a structure).
    let mut final_plus  = BlastIntervalTree::new(0, qa_len + 1, 0, s_len + 1);
    let mut final_minus = BlastIntervalTree::new(0, qa_len + 1, 0, s_len + 1);
    results.retain(|r| {
        let plus3 = r.hsp.strand == Strand::Plus;
        let contained = {
            let tree = if plus3 { &final_plus } else { &final_minus };
            tree.contains(r.hsp.q_start, r.hsp.q_end, r.hsp.s_start, r.hsp.s_end,
                          r.hsp.score, plus3, MIN_DIAG_SEP)
        };
        if !contained {
            let tree = if plus3 { &mut final_plus } else { &mut final_minus };
            tree.add_simple(ITreeHsp { q_start: r.hsp.q_start, q_end: r.hsp.q_end,
                                s_start: r.hsp.s_start, s_end: r.hsp.s_end,
                                score: r.hsp.score });
        }
        !contained
    });

    results
}

/// Remove prelim HSPs that share a start or end point with a higher-scoring prelim HSP.
/// Mirrors the call to Blast_HSPListPurgeHSPsWithCommonEndpoints in blast_engine.c:545,
/// which runs after Phase 2a and before Phase 2b traceback.
///
/// NCBI compares (context, query.offset, subject.offset) for start-dedup and
/// (context, query.end, subject.end) for end-dedup. In Rust PrelimHsp coordinates
/// (all in strand-specific space; minus-strand s coords are RC-TE):
///   Plus  strand: q_off = q_start, s_off = s_start; q_end = q_end, s_end = s_end
///   Minus strand: q_off = q_end,   s_off = sa_len - s_end (FWD-TE start)
///                 q_end = q_start, s_end = sa_len - s_start (FWD-TE end)
fn purge_prelim_common_endpoints(prelim_hsps: &mut Vec<PrelimHsp>, sa_len: u32) {
    let n = prelim_hsps.len();
    if n < 2 { return; }

    let mut keep = vec![true; n];

    // FWD-TE start/end helpers for minus strand (convert RC-TE to FWD-TE).
    let fwd_s_start = |h: &PrelimHsp| match h.strand {
        Strand::Plus  => h.s_start,
        Strand::Minus => sa_len - h.s_end,
    };
    let fwd_s_end = |h: &PrelimHsp| match h.strand {
        Strand::Plus  => h.s_end,
        Strand::Minus => sa_len - h.s_start,
    };
    // NCBI's query.offset (ascending = left of alignment) per strand:
    let q_off = |h: &PrelimHsp| match h.strand {
        Strand::Plus  => h.q_start,
        Strand::Minus => h.q_end,
    };
    // NCBI's query.end (ascending = right of alignment) per strand:
    let q_end_key = |h: &PrelimHsp| match h.strand {
        Strand::Plus  => h.q_end,
        Strand::Minus => h.q_start,
    };
    let strand_ctx = |h: &PrelimHsp| if h.strand == Strand::Plus { 0u8 } else { 1u8 };

    // Pass 1: deduplicate by common start (query.offset, subject.offset).
    // Sort: context ASC, q_off ASC, s_off ASC, score DESC, q_end DESC, s_end DESC.
    // When all keys tie, preserve insertion order (index ASC) — matches NCBI's qsort
    // practical behavior where equal elements are not reordered.
    {
        let mut idx: Vec<usize> = (0..n).collect();
        idx.sort_unstable_by(|&a, &b| {
            let ha = &prelim_hsps[a];
            let hb = &prelim_hsps[b];
            // NCBI s_QueryOffsetCompareHSPs: after (ctx, q_off, s_off, score DESC),
            // tiebreaker is query.end DESC.  For plus strand, query.end = FWD_qe → DESC.
            // For minus strand, NCBI query.end = qlen - FWD_qs → DESC means SMALLER
            // FWD_qs (more leftward start) wins, so we sort ASC on q_end_key (= FWD_qs).
            // NCBI s_QueryOffsetCompareHSPs: after (ctx, q_off, s_off, score DESC),
            // tiebreaker is query.end DESC.  For plus strand, query.end = FWD_qe → DESC.
            // For minus strand, NCBI query.end = qlen - FWD_qs → DESC means SMALLER
            // FWD_qs (more leftward start) wins, so we sort ASC on q_end_key (= FWD_qs).
            let qek_cmp = match ha.strand {
                Strand::Plus  => q_end_key(hb).cmp(&q_end_key(ha)), // DESC on FWD_qe
                Strand::Minus => q_end_key(ha).cmp(&q_end_key(hb)), // ASC  on FWD_qs
            };
            strand_ctx(ha).cmp(&strand_ctx(hb))
                .then_with(|| q_off(ha).cmp(&q_off(hb)))
                .then_with(|| fwd_s_start(ha).cmp(&fwd_s_start(hb)))
                .then_with(|| hb.score.cmp(&ha.score))
                .then_with(|| qek_cmp)
                .then_with(|| fwd_s_end(hb).cmp(&fwd_s_end(ha)))
                .then_with(|| a.cmp(&b))
        });
        let mut i = 0;
        while i < n {
            if !keep[idx[i]] { i += 1; continue; }
            let h0 = &prelim_hsps[idx[i]];
            let (ctx0, qo0, so0) = (strand_ctx(h0), q_off(h0), fwd_s_start(h0));
            let mut j = i + 1;
            while j < n {
                let hj = &prelim_hsps[idx[j]];
                if strand_ctx(hj) != ctx0 || q_off(hj) != qo0 || fwd_s_start(hj) != so0 { break; }
                keep[idx[j]] = false;
                j += 1;
            }
            i = j;
        }
    }

    // Pass 2: deduplicate by common end (query.end, subject.end).
    // Sort: context ASC, q_end ASC, s_end ASC, score DESC, q_off DESC, s_off DESC.
    // Same insertion-order tiebreaker as Pass 1.
    {
        let mut idx: Vec<usize> = (0..n).collect();
        idx.sort_unstable_by(|&a, &b| {
            let ha = &prelim_hsps[a];
            let hb = &prelim_hsps[b];
            // NCBI s_QueryEndCompareHSPs tiebreaker after (ctx, query.end, subject.end):
            // score DESC, then query.offset DESC, then subject.offset DESC.
            // For minus strand NCBI query.offset = qlen - FWD_qe, so query.offset DESC means
            // SMALLER FWD_qe first → ASC on q_off (= FWD_qe).  Plus strand: q_off = FWD_qs → DESC.
            // (Mirrors the strand-aware qek_cmp used in Pass 1.)
            let qoff_cmp = match ha.strand {
                Strand::Plus  => q_off(hb).cmp(&q_off(ha)), // DESC on FWD_qs
                Strand::Minus => q_off(ha).cmp(&q_off(hb)), // ASC  on FWD_qe
            };
            strand_ctx(ha).cmp(&strand_ctx(hb))
                .then_with(|| q_end_key(ha).cmp(&q_end_key(hb)))
                .then_with(|| fwd_s_end(ha).cmp(&fwd_s_end(hb)))
                .then_with(|| hb.score.cmp(&ha.score))
                .then_with(|| qoff_cmp)
                .then_with(|| fwd_s_start(hb).cmp(&fwd_s_start(ha)))
                .then_with(|| a.cmp(&b))
        });
        let mut i = 0;
        while i < n {
            if !keep[idx[i]] { i += 1; continue; }
            let h0 = &prelim_hsps[idx[i]];
            let (ctx0, qe0, se0) = (strand_ctx(h0), q_end_key(h0), fwd_s_end(h0));
            let mut j = i + 1;
            while j < n {
                let hj = &prelim_hsps[idx[j]];
                if strand_ctx(hj) != ctx0 || q_end_key(hj) != qe0 || fwd_s_end(hj) != se0 { break; }
                keep[idx[j]] = false;
                j += 1;
            }
            i = j;
        }
    }

    let mut j = 0;
    for i in 0..n {
        if keep[i] { prelim_hsps.swap(i, j); j += 1; }
    }
    prelim_hsps.truncate(j);
}

// ── Common-endpoint cut helpers ───────────────────────────────────────────────

/// Walk from the FRONT of `ops`, find the first point where qid >= need_q AND sid >= need_s,
/// and KEEP the prefix (remove the suffix).  Mirrors NCBI `s_CutOffGapEditScript` cut_begin=FALSE.
/// Returns (ok, actual_qid_consumed, actual_sid_consumed).
fn cut_ops_keep_front(ops: &mut Vec<(EditOp, u32)>, need_q: u32, need_s: u32) -> (bool, u32, u32) {
    if need_q == 0 && need_s == 0 { return (true, 0, 0); }
    let mut gq = 0u32;
    let mut gs = 0u32;
    let mut cut_i = 0usize;
    let mut partial = 0u32;
    let mut found = false;
    'scan: for i in 0..ops.len() {
        let (op, cnt) = ops[i];
        match op {
            EditOp::Sub => {
                for j in 0..cnt {
                    gq += 1; gs += 1;
                    if gq >= need_q && gs >= need_s {
                        cut_i = i; partial = j + 1; found = true; break 'scan;
                    }
                }
            }
            EditOp::GapInQuery => {
                gs += cnt;
                if gq >= need_q && gs >= need_s { cut_i = i; partial = cnt; found = true; break; }
            }
            EditOp::GapInSubject => {
                gq += cnt;
                if gq >= need_q && gs >= need_s { cut_i = i; partial = cnt; found = true; break; }
            }
        }
    }
    if !found { return (false, gq, gs); }
    let (op, cnt) = ops[cut_i];
    if partial < cnt {
        ops[cut_i] = (op, partial);
    }
    ops.truncate(cut_i + 1);
    (true, gq, gs)
}

/// Remove the leading prefix of `ops` that consumes exactly `need_q` query-bases
/// AND `need_s` subject-bases (NCBI `s_CutOffGapEditScript` cut_begin=TRUE for plus strand).
/// Returns (ok, actual_gq_consumed, actual_gs_consumed).
/// actual_gs_consumed can exceed need_s when GapInQuery ops span the cut point.
fn cut_ops_front(ops: &mut Vec<(EditOp, u32)>, need_q: u32, need_s: u32) -> (bool, u32, u32) {
    if need_q == 0 && need_s == 0 { return (true, 0, 0); }
    let mut gq = 0u32;
    let mut gs = 0u32;
    let mut cut_i = 0usize;
    let mut partial = 0u32;
    let mut found = false;
    'scan: for i in 0..ops.len() {
        let (op, cnt) = ops[i];
        match op {
            EditOp::Sub => {
                for j in 0..cnt {
                    gq += 1; gs += 1;
                    if gq >= need_q && gs >= need_s {
                        cut_i = i; partial = j + 1; found = true; break 'scan;
                    }
                }
            }
            EditOp::GapInQuery => {
                gs += cnt;
                if gq >= need_q && gs >= need_s { cut_i = i; partial = cnt; found = true; break; }
            }
            EditOp::GapInSubject => {
                gq += cnt;
                if gq >= need_q && gs >= need_s { cut_i = i; partial = cnt; found = true; break; }
            }
        }
    }
    if !found { return (false, gq, gs); }
    let (op, cnt) = ops[cut_i];
    let rem = cnt - partial;
    ops.drain(0..=cut_i);
    if rem > 0 { ops.insert(0, (op, rem)); }
    (true, gq, gs)
}

/// Remove the trailing suffix of `ops` that consumes `need_q` query-bases AND
/// `need_s` subject-bases from the back (NCBI cut_begin=TRUE for minus strand,
/// or cut_begin=FALSE for plus strand common-end purge).
/// Returns (ok, actual_gq_consumed, actual_gs_consumed).
/// actual_gs_consumed can exceed need_s when GapInQuery ops span the cut point.
fn cut_ops_back(ops: &mut Vec<(EditOp, u32)>, need_q: u32, need_s: u32) -> (bool, u32, u32) {
    if need_q == 0 && need_s == 0 { return (true, 0, 0); }
    let n = ops.len();
    let mut gq = 0u32;
    let mut gs = 0u32;
    let mut cut_i = n;
    let mut partial = 0u32;
    let mut found = false;
    'scan: for i in (0..n).rev() {
        let (op, cnt) = ops[i];
        match op {
            EditOp::Sub => {
                for j in 0..cnt {
                    gq += 1; gs += 1;
                    if gq >= need_q && gs >= need_s {
                        cut_i = i; partial = j + 1; found = true; break 'scan;
                    }
                }
            }
            EditOp::GapInQuery => {
                gs += cnt;
                if gq >= need_q && gs >= need_s { cut_i = i; partial = cnt; found = true; break; }
            }
            EditOp::GapInSubject => {
                gq += cnt;
                if gq >= need_q && gs >= need_s { cut_i = i; partial = cnt; found = true; break; }
            }
        }
    }
    if !found { return (false, gq, gs); }
    let (op, cnt) = ops[cut_i];
    let keep = cnt - partial;
    ops.truncate(cut_i);
    if keep > 0 { ops.push((op, keep)); }
    (true, gq, gs)
}

/// Walk `ops` from the BACK, consuming `need_q` query-bases AND `need_s` subject-bases,
/// then REMOVE the front (overlapping) portion and KEEP the back (non-overlapping) portion.
///
/// Mirrors NCBI's `s_CutOffGapEditScript` with cut_begin=FALSE applied from the
/// non-overlapping side of a minus-strand common-END secondary (the extension beyond
/// the primary).  Gap ops are processed whole-run at a time, same as NCBI.
///
/// Returns (ok, actual_gq_consumed, actual_gs_consumed).
fn keep_ops_back(ops: &mut Vec<(EditOp, u32)>, need_q: u32, need_s: u32) -> (bool, u32, u32) {
    if need_q == 0 && need_s == 0 { return (true, 0, 0); }
    let mut gq = 0u32;
    let mut gs = 0u32;
    let mut cut_i = 0usize;
    let mut partial = 0u32;
    let mut found = false;
    let n = ops.len();
    'scan: for i in (0..n).rev() {
        let (op, cnt) = ops[i];
        match op {
            EditOp::Sub => {
                for j in 0..cnt {
                    gq += 1; gs += 1;
                    if gq >= need_q && gs >= need_s {
                        cut_i = i; partial = j + 1; found = true; break 'scan;
                    }
                }
            }
            EditOp::GapInQuery => {
                gs += cnt;
                if gq >= need_q && gs >= need_s { cut_i = i; partial = cnt; found = true; break; }
            }
            EditOp::GapInSubject => {
                gq += cnt;
                if gq >= need_q && gs >= need_s { cut_i = i; partial = cnt; found = true; break; }
            }
        }
    }
    if !found { return (false, gq, gs); }
    let (op, cnt) = ops[cut_i];
    ops.drain(0..cut_i);
    if partial < cnt {
        // Partial cut only valid for Sub ops (same constraint as NCBI's ASSERT).
        ops[0] = (op, partial);
    }
    (true, gq, gs)
}

/// Mirrors `Blast_HSPReevaluateWithAmbiguitiesGapped`.
///
/// `q_seq` / `s_seq` are the pre-extracted aligned columns (BLASTNA, gap=15) for
/// this HSP's edit script (obtained from `extract_aligned`).
///
/// Works in "script-space" coordinates: positions advance from (0,0) along the ops.
/// For plus strand: script-q maps 1:1 to original-query; script-s maps 1:1 to original-subject.
/// For minus strand: script-q maps to original-query (same direction); script-s maps
/// to subj_rc (sa_len − hsp.s_end .. sa_len − hsp.s_start), so `sa_len` is needed
/// to convert the s-delta back to forward-subject coordinates.
///
/// Returns `true` if the HSP should be deleted (score < min_score after re-evaluation).
fn reevaluate_gapped(
    hsp: &mut Hsp,
    q_seq: &[u8],
    s_seq: &[u8],
    matrix: &ScoreMatrix,
    gap_open: i32,
    gap_extend: i32,
    min_score: i32,
    is_minus: bool,
    raw_q: &[u8],
    q_raw_start: usize,
    raw_s: &[u8],
    s_raw_start: usize,
) -> bool {
    if hsp.edit_script.ops.is_empty() { return true; }

    let mut sum = 0i32;
    let mut score = 0i32;
    let mut seq_i = 0usize; // aligned column index into q_seq / s_seq

    // Non-gap positions consumed so far (track deltas from script start)
    let mut q_pos = 0u32; // non-gap query bases consumed (= script-q distance from start)
    let mut s_pos = 0u32; // non-gap subject bases consumed (= script-s distance from start)

    // Best and current segment tracking
    let mut best_q_start = 0u32;
    let mut best_s_start = 0u32;
    let mut best_q_end   = 0u32;
    let mut best_s_end   = 0u32;
    let mut cur_q_start  = 0u32;
    let mut cur_s_start  = 0u32;

    let mut best_start_i   = 0usize;
    let mut best_end_i     = 0usize;
    let mut best_end_num   = 0u32;
    let mut cur_start_i    = 0usize;
    let mut cur_start_skip = 0u32; // units to skip at ops[cur_start_i]
    let mut best_start_skip = 0u32;

    let ops_snap: Vec<(EditOp, u32)> = hsp.edit_script.ops.clone();
    let esp_len = ops_snap.len();

    // NCBI Blast_HSPReevaluateWithAmbiguitiesGapped scans the alignment in
    // SUBJECT-FORWARD order.  For plus strand that equals query-forward (our native
    // script order).  For minus strand, subject-forward = query-BACKWARD, so the
    // best-subsegment Kadane scan must run over a REVERSED view of the alignment;
    // otherwise a zero-net tail at the segment boundary is kept/dropped on the wrong
    // side (e.g. L1PBa1_5end @42117052).  We scan the reversed columns/ops here, then
    // convert the resulting segment markers back to forward indexing so the existing
    // trim/extend/coord code (below) is unchanged.  The aligned arrays are tiny (one
    // HSP), so reversing them per call is cheap — the huge raw genome is never reversed.
    let scan_reverse = is_minus;
    let rev_ops;
    let rev_q;
    let rev_s;
    let (scan_ops, scan_q, scan_s): (&[(EditOp, u32)], &[u8], &[u8]) = if scan_reverse {
        rev_ops = ops_snap.iter().rev().cloned().collect::<Vec<_>>();
        rev_q = q_seq.iter().rev().cloned().collect::<Vec<u8>>();
        rev_s = s_seq.iter().rev().cloned().collect::<Vec<u8>>();
        (&rev_ops, &rev_q, &rev_s)
    } else {
        (&ops_snap, q_seq, s_seq)
    };

    for (ei, &(op, cnt)) in scan_ops.iter().enumerate() {
        match op {
            EditOp::Sub => {
                for k in 0..cnt {
                    let qb0 = scan_q.get(seq_i).copied().unwrap_or(14) & 0xf;
                    let sb0 = scan_s.get(seq_i).copied().unwrap_or(14) & 0xf;
                    // Minus strand: the cut/reevaluate path aligns FWD-query (q_core) vs
                    // RC-subject (rc_core), i.e. matrix.scores[q_fwd][comp(s_fwd)].  The
                    // preliminary DP (and NCBI) score minus in the (RC_query, FWD_subject)
                    // frame = matrix.scores[comp(q_fwd)][s_fwd].  For NON-complement-symmetric
                    // matrices (the RepeatMasker p##g family, e.g. 20p37g[G][A]=-8 vs [C][T]=-7)
                    // these two frames differ, so reevaluated/cut minus HSPs were scored a few
                    // points off NCBI.  Complement both bases for minus to score in NCBI's frame
                    // (no-op for ACGT-symmetric matrices like comparison.matrix). Same class as
                    // bug #30 (ungapped extension frame).
                    let (qb, sb) = if is_minus {
                        (BLASTNA_COMPLEMENT[qb0 as usize] & 0xf, BLASTNA_COMPLEMENT[sb0 as usize] & 0xf)
                    } else {
                        (qb0, sb0)
                    };
                    sum += matrix.scores[qb as usize][sb as usize];
                    seq_i += 1; q_pos += 1; s_pos += 1;
                    let op_idx = k + 1; // units consumed in this op so far

                    if sum < 0 {
                        // New run starts after this base
                        if op_idx < cnt {
                            cur_start_i = ei;
                            cur_start_skip = op_idx;
                        } else {
                            cur_start_i = ei + 1;
                            cur_start_skip = 0;
                        }
                        sum = 0;
                        cur_q_start = q_pos;
                        cur_s_start = s_pos;
                        if score < min_score {
                            best_q_start = q_pos; best_s_start = s_pos; score = 0;
                            best_start_i = cur_start_i; best_start_skip = cur_start_skip;
                            best_end_i = cur_start_i; best_end_num = 0;
                        }
                    } else if sum > score {
                        score = sum;
                        best_q_start = cur_q_start; best_s_start = cur_s_start;
                        best_q_end = q_pos; best_s_end = s_pos;
                        best_start_i = cur_start_i; best_start_skip = cur_start_skip;
                        best_end_i = ei; best_end_num = op_idx;
                    }
                }
            }
            EditOp::GapInQuery => {
                sum -= gap_open + gap_extend * cnt as i32;
                for _ in 0..cnt { seq_i += 1; s_pos += 1; }
                if sum < 0 {
                    cur_start_i = ei + 1; cur_start_skip = 0;
                    sum = 0; cur_q_start = q_pos; cur_s_start = s_pos;
                    if score < min_score {
                        best_q_start = q_pos; best_s_start = s_pos; score = 0;
                        best_start_i = cur_start_i; best_start_skip = 0;
                        best_end_i = cur_start_i; best_end_num = 0;
                    }
                } else if sum > score {
                    score = sum;
                    best_q_start = cur_q_start; best_s_start = cur_s_start;
                    best_q_end = q_pos; best_s_end = s_pos;
                    best_start_i = cur_start_i; best_start_skip = cur_start_skip;
                    best_end_i = ei; best_end_num = cnt;
                }
            }
            EditOp::GapInSubject => {
                sum -= gap_open + gap_extend * cnt as i32;
                for _ in 0..cnt { seq_i += 1; q_pos += 1; }
                if sum < 0 {
                    cur_start_i = ei + 1; cur_start_skip = 0;
                    sum = 0; cur_q_start = q_pos; cur_s_start = s_pos;
                    if score < min_score {
                        best_q_start = q_pos; best_s_start = s_pos; score = 0;
                        best_start_i = cur_start_i; best_start_skip = 0;
                        best_end_i = cur_start_i; best_end_num = 0;
                    }
                } else if sum > score {
                    score = sum;
                    best_q_start = cur_q_start; best_s_start = cur_s_start;
                    best_q_end = q_pos; best_s_end = s_pos;
                    best_start_i = cur_start_i; best_start_skip = cur_start_skip;
                    best_end_i = ei; best_end_num = cnt;
                }
            }
        }
    }

    if score < min_score { return true; }

    // Convert the reversed-frame best-segment markers back to forward indexing so the
    // trim/extend/coord code below (which operates on the forward `ops_snap` /
    // `hsp.edit_script.ops` and forward raw sequences) works unchanged.
    if scan_reverse {
        // Total query / subject bases in the alignment.
        let mut qtot = 0u32;
        let mut stot = 0u32;
        for &(op, c) in &ops_snap {
            match op {
                EditOp::Sub => { qtot += c; stot += c; }
                EditOp::GapInSubject => { qtot += c; }
                EditOp::GapInQuery => { stot += c; }
            }
        }
        // Reversed op index r maps to forward op index (esp_len-1-r).  The reversed
        // segment [best_start_i, best_end_i] becomes forward [f_start_i, f_end_i].
        // best_*_i point at Sub ops (gaps never improve the running score).
        let f_start_i = esp_len - 1 - best_end_i;
        let f_end_i   = esp_len - 1 - best_start_i;
        // Front-skip on the reversed start op = back-skip on the forward END op.
        let f_end_num = ops_snap[f_end_i].1 - best_start_skip;
        // Units kept from front of reversed end op = kept from back of forward START op.
        let f_start_skip = ops_snap[f_start_i].1 - best_end_num;
        let f_q_start = qtot - best_q_end;
        let f_q_end   = qtot - best_q_start;
        let f_s_start = stot - best_s_end;
        let f_s_end   = stot - best_s_start;
        best_start_i = f_start_i;
        best_end_i = f_end_i;
        best_start_skip = f_start_skip;
        best_end_num = f_end_num;
        best_q_start = f_q_start;
        best_q_end = f_q_end;
        best_s_start = f_s_start;
        best_s_end = f_s_end;
    }

    // Trim edit script to [best_start_i, best_end_i] with adjustments
    {
        let ops = &mut hsp.edit_script.ops;
        // Trim tail
        if best_end_i < esp_len.saturating_sub(1) || best_end_num < ops_snap[best_end_i].1 {
            ops.truncate(best_end_i + 1);
            if best_end_num > 0 { ops[best_end_i].1 = best_end_num; } else { ops.truncate(best_end_i); }
        }
        // Trim head
        if best_start_i > 0 {
            ops.drain(0..best_start_i);
        }
        // Trim within first op
        if best_start_skip > 0 && !ops.is_empty() {
            let (op0, cnt0) = ops[0];
            if best_start_skip < cnt0 {
                ops[0] = (op0, cnt0 - best_start_skip);
            } else {
                ops.remove(0);
            }
        }
    }

    // Post-processing: extend the trimmed best segment with exact raw-sequence matches,
    // mirroring NCBI's Blast_HSPReevaluateWithAmbiguitiesGapped extension step.
    // Use i64 so the left extension can go before script position 0 without wrapping.
    let mut ext_bqs = best_q_start as i64;
    let mut ext_bss = best_s_start as i64;
    let mut ext_bqe = best_q_end as i64;
    let mut ext_bse = best_s_end as i64;
    {
        // Left extension: walk backward through raw sequences from ext_bqs/ext_bss.
        loop {
            let qi = q_raw_start as i64 + ext_bqs;
            let si = s_raw_start as i64 + ext_bss;
            if qi <= 0 || si <= 0 { break; }
            let qb = raw_q[(qi - 1) as usize];
            let sb = raw_s[(si - 1) as usize];
            if qb < 4 && qb == sb {
                // score_params->reward == 0 for rmblastn (matrix-only scoring),
                // so extension updates coordinates/script but does not change score.
                ext_bqs -= 1;
                ext_bss -= 1;
                let ops = &mut hsp.edit_script.ops;
                if ops.first().map(|(op, _)| *op) == Some(EditOp::Sub) {
                    ops[0].1 += 1;
                } else {
                    ops.insert(0, (EditOp::Sub, 1));
                }
            } else {
                break;
            }
        }

        // Right extension: walk forward through raw sequences from ext_bqe/ext_bse.
        loop {
            let qi = q_raw_start as i64 + ext_bqe;
            let si = s_raw_start as i64 + ext_bse;
            if qi >= raw_q.len() as i64 || si >= raw_s.len() as i64 { break; }
            let qb = raw_q[qi as usize];
            let sb = raw_s[si as usize];
            if qb < 4 && qb == sb {
                // score_params->reward == 0 for rmblastn (matrix-only scoring),
                // so extension updates coordinates/script but does not change score.
                ext_bqe += 1;
                ext_bse += 1;
                let ops = &mut hsp.edit_script.ops;
                if ops.last().map(|(op, _)| *op) == Some(EditOp::Sub) {
                    ops.last_mut().unwrap().1 += 1;
                } else {
                    ops.push((EditOp::Sub, 1));
                }
            } else {
                break;
            }
        }
    }

    // Update hsp coordinates using i64 arithmetic (ext_bqs/ext_bss can be negative).
    hsp.score = score;

    if is_minus {
        // For minus strand: q_start increases (same as plus), s is mirrored.
        // Script s_origin = sa_len - old_hsp.s_end.  After trimming ext_bss from front
        // and (total_s - ext_bse) from back:
        //   new s_end (hsp fwd) = old_hsp.s_end - ext_bss
        //   new s_start (hsp fwd) = old_hsp.s_end - ext_bse
        let old_q_start = hsp.q_start as i64;
        let old_s_end   = hsp.s_end as i64;
        hsp.q_start = (old_q_start + ext_bqs).max(0) as u32;
        hsp.q_end   = (old_q_start + ext_bqe).max(0) as u32;
        hsp.s_end   = (old_s_end - ext_bss).max(0) as u32;
        hsp.s_start = (old_s_end - ext_bse).max(0) as u32;
    } else {
        let old_q_start = hsp.q_start as i64;
        let old_s_start = hsp.s_start as i64;
        let new_q_start = (old_q_start + ext_bqs).max(0) as u32;
        let new_s_start = (old_s_start + ext_bss).max(0) as u32;
        hsp.q_start = new_q_start;
        hsp.q_end   = new_q_start + (ext_bqe - ext_bqs) as u32;
        hsp.s_start = new_s_start;
        hsp.s_end   = new_s_start + (ext_bse - ext_bss) as u32;
    }

    false
}

// ── Faithful port of Blast_HSPListPurgeHSPsWithCommonEndpoints ──────────────
//
// NCBI (blast_hits.c:2576) operates IN PLACE on `BlastHSP** hsp_array` with a
// `purge` flag: FALSE = cut secondaries (keep the trimmed remainder), TRUE =
// delete them.  It runs two sub-passes — common START (query.offset,
// subject.offset) then common END (query.end, subject.end) — and in each pass
// MOVES every removed HSP (cut or freed) to the end of the shrinking array.
// The caller (blast_traceback.c:658/688) calls it twice: FALSE (cut) → the
// returned count is `extra_start`, and the cut remainders sit at
// hsp_array[extra_start..hspcnt] → re-traced (reevaluated) in array order →
// PurgeNull → TRUE (delete).
//
// Earlier this was approximated by pulling cut pieces into a separate vec and
// `.reverse()`-ing it to coax the keep-first dedup into NCBI's selection (#32B,
// behaviourally correct but structurally back-worked).  The functions below are
// the faithful structure: `None` models a freed/NULL slot, and the move-to-end
// ordering is reproduced exactly so the re-trace order emerges naturally.

/// NCBI `s_QueryOffsetCompareHSPs` comparator (blast_hits.c).  Strand ASC,
/// query.offset ASC, subject.offset ASC, score DESC, query.end DESC,
/// subject.end DESC.  A stable sort supplies NCBI's quasi-stable qsort tie
/// behaviour for the small per-subject arrays, so no explicit index key is
/// needed.
fn cmp_query_offset(a: &AlignResult, b: &AlignResult) -> std::cmp::Ordering {
    let ha = &a.hsp; let hb = &b.hsp;
    let sa = ha.strand == Strand::Minus;
    let sb = hb.strand == Strand::Minus;
    if sa != sb { return sa.cmp(&sb); }
    let kqa = match ha.strand { Strand::Plus => ha.q_start, Strand::Minus => ha.q_end };
    let kqb = match hb.strand { Strand::Plus => hb.q_start, Strand::Minus => hb.q_end };
    kqa.cmp(&kqb)
        .then_with(|| ha.s_start.cmp(&hb.s_start))   // subject.offset ASC
        .then_with(|| hb.score.cmp(&ha.score))       // score DESC
        .then_with(|| match ha.strand {              // query.end DESC
            Strand::Plus  => hb.q_end.cmp(&ha.q_end),
            Strand::Minus => ha.q_start.cmp(&hb.q_start),
        })
        .then_with(|| hb.s_end.cmp(&ha.s_end))        // subject.end DESC
}

/// NCBI `s_QueryEndCompareHSPs` comparator.  Strand ASC, query.end ASC,
/// subject.end ASC, score DESC, query.offset DESC, subject.offset DESC.
fn cmp_query_end(a: &AlignResult, b: &AlignResult) -> std::cmp::Ordering {
    let ha = &a.hsp; let hb = &b.hsp;
    let sa = ha.strand == Strand::Minus;
    let sb = hb.strand == Strand::Minus;
    if sa != sb { return sa.cmp(&sb); }
    let kqa = match ha.strand { Strand::Plus => ha.q_end, Strand::Minus => ha.q_start };
    let kqb = match hb.strand { Strand::Plus => hb.q_end, Strand::Minus => hb.q_start };
    kqa.cmp(&kqb)
        .then_with(|| ha.s_end.cmp(&hb.s_end))         // subject.end ASC
        .then_with(|| hb.score.cmp(&ha.score))         // score DESC
        .then_with(|| match ha.strand {                // query.offset DESC
            Strand::Plus  => hb.q_start.cmp(&ha.q_start),
            Strand::Minus => ha.q_end.cmp(&hb.q_end),
        })
        .then_with(|| hb.s_start.cmp(&ha.s_start))      // subject.offset DESC
}

/// Whether two HSPs share a common START (NCBI: context + query.offset +
/// subject.offset + subject.frame).
fn same_group_start(a: &AlignResult, b: &AlignResult) -> bool {
    if a.hsp.strand != b.hsp.strand { return false; }
    let ka = match a.hsp.strand { Strand::Plus => a.hsp.q_start, Strand::Minus => a.hsp.q_end };
    let kb = match b.hsp.strand { Strand::Plus => b.hsp.q_start, Strand::Minus => b.hsp.q_end };
    ka == kb && a.hsp.s_start == b.hsp.s_start
}

/// Whether two HSPs share a common END (NCBI: context + query.end +
/// subject.end + subject.frame).
fn same_group_end(a: &AlignResult, b: &AlignResult) -> bool {
    if a.hsp.strand != b.hsp.strand { return false; }
    let ka = match a.hsp.strand { Strand::Plus => a.hsp.q_end, Strand::Minus => a.hsp.q_start };
    let kb = match b.hsp.strand { Strand::Plus => b.hsp.q_end, Strand::Minus => b.hsp.q_start };
    ka == kb && a.hsp.s_end == b.hsp.s_end
}

/// Cut a common-START secondary against `primary` (NCBI: `!purge && query.end >
/// primary.query.end` → `s_CutOffGapEditScript(..., cut_begin=TRUE)`).  Returns
/// `true` if the secondary should be KEPT (trimmed remainder), `false` if it
/// must be freed (it does not extend past the primary, or the cut emptied it).
fn cut_start_secondary(primary: &AlignResult, sec: &mut AlignResult) -> bool {
    let strand_p = primary.hsp.strand;
    let secondary_extends = match strand_p {
        Strand::Plus  => sec.hsp.q_end > primary.hsp.q_end,
        Strand::Minus => sec.hsp.q_start < primary.hsp.q_start,
    };
    if !secondary_extends { return false; }
    match strand_p {
        Strand::Plus => {
            let dq = primary.hsp.q_end.saturating_sub(sec.hsp.q_start);
            let ds = primary.hsp.s_end.saturating_sub(sec.hsp.s_start);
            let old_q_start = sec.hsp.q_start;
            let old_s_start = sec.hsp.s_start;
            let (ok, gq, gs) = cut_ops_front(&mut sec.hsp.edit_script.ops, dq, ds);
            if ok { sec.hsp.q_start = old_q_start + gq; sec.hsp.s_start = old_s_start + gs; }
        }
        Strand::Minus => {
            let dq = sec.hsp.q_end.saturating_sub(primary.hsp.q_start);
            let ds = primary.hsp.s_end.saturating_sub(sec.hsp.s_start);
            let old_s_start = sec.hsp.s_start;
            let (ok, _gq, gs) = cut_ops_back(&mut sec.hsp.edit_script.ops, dq, ds);
            if ok { sec.hsp.q_end = primary.hsp.q_start; sec.hsp.s_start = old_s_start + gs; }
        }
    }
    !sec.hsp.edit_script.ops.is_empty()
}

/// Cut a common-END secondary against `primary` (NCBI: `!purge && query.offset <
/// primary.query.offset` → `s_CutOffGapEditScript(..., cut_begin=FALSE)`).
fn cut_end_secondary(primary: &AlignResult, sec: &mut AlignResult) -> bool {
    let strand_p = primary.hsp.strand;
    let secondary_extends = match strand_p {
        Strand::Plus  => sec.hsp.q_start < primary.hsp.q_start,
        Strand::Minus => sec.hsp.q_end   > primary.hsp.q_end,
    };
    if !secondary_extends { return false; }
    match strand_p {
        Strand::Plus => {
            let need_q = primary.hsp.q_start.saturating_sub(sec.hsp.q_start);
            let need_s = primary.hsp.s_start.saturating_sub(sec.hsp.s_start);
            let (ok, qid, sid) = cut_ops_keep_front(&mut sec.hsp.edit_script.ops, need_q, need_s);
            if ok { sec.hsp.q_end = sec.hsp.q_start + qid; sec.hsp.s_end = sec.hsp.s_start + sid; }
        }
        Strand::Minus => {
            let need_q = sec.hsp.q_end.saturating_sub(primary.hsp.q_end);
            let need_s = primary.hsp.s_start.saturating_sub(sec.hsp.s_start);
            let (ok, gq, gs) = keep_ops_back(&mut sec.hsp.edit_script.ops, need_q, need_s);
            if ok {
                sec.hsp.q_start = sec.hsp.q_end.saturating_sub(gq);
                sec.hsp.s_end = sec.hsp.s_start + gs;
            }
        }
    }
    !sec.hsp.edit_script.ops.is_empty()
}

/// Faithful port of `Blast_HSPListPurgeHSPsWithCommonEndpoints` (blast_hits.c:2576).
///
/// `arr` is NCBI's `hsp_array`; `None` models a freed/NULL slot.  Survivors end
/// up in `arr[0..count)`; removed HSPs (cut remainders as `Some`, freed as
/// `None`) are moved to `arr[count..]` exactly as NCBI's shift-to-end does, so a
/// `purge=false` caller can re-trace the cut remainders in array order.  Returns
/// the survivor count.
fn purge_hsps_with_common_endpoints(arr: &mut [Option<AlignResult>], purge: bool) -> usize {
    let mut hsp_count = arr.len();
    if hsp_count < 2 { return hsp_count; }

    // Pass 1: common START.  qsort s_QueryOffsetCompareHSPs over arr[0..hsp_count).
    arr[..hsp_count].sort_by(|a, b| cmp_query_offset(a.as_ref().unwrap(), b.as_ref().unwrap()));
    let mut i = 0;
    while i < hsp_count {
        // NCBI's inner `while` never increments j: it repeatedly collapses the
        // element at i+1 as later elements shift down into it.
        while i + 1 < hsp_count
            && same_group_start(arr[i].as_ref().unwrap(), arr[i + 1].as_ref().unwrap())
        {
            hsp_count -= 1;
            let mut sec = arr[i + 1].take().unwrap();
            let keep = if !purge { cut_start_secondary(arr[i].as_ref().unwrap(), &mut sec) } else { false };
            let moved = if keep { Some(sec) } else { None };
            for k in (i + 1)..hsp_count {
                arr[k] = arr[k + 1].take();
            }
            arr[hsp_count] = moved;
        }
        i += 1;
    }

    // Pass 2: common END.  qsort s_QueryEndCompareHSPs over the survivors.
    arr[..hsp_count].sort_by(|a, b| cmp_query_end(a.as_ref().unwrap(), b.as_ref().unwrap()));
    let mut i = 0;
    while i < hsp_count {
        while i + 1 < hsp_count
            && same_group_end(arr[i].as_ref().unwrap(), arr[i + 1].as_ref().unwrap())
        {
            hsp_count -= 1;
            let mut sec = arr[i + 1].take().unwrap();
            let keep = if !purge { cut_end_secondary(arr[i].as_ref().unwrap(), &mut sec) } else { false };
            let moved = if keep { Some(sec) } else { None };
            for k in (i + 1)..hsp_count {
                arr[k] = arr[k + 1].take();
            }
            arr[hsp_count] = moved;
        }
        i += 1;
    }

    hsp_count
}

/// Remove HSPs that share a start or end point with a higher-scoring HSP.
/// Mirrors Blast_HSPListPurgeHSPsWithCommonEndpoints (blast_hits.c:2486).
///
/// NCBI sorts by (context-query.offset, subject.offset) and removes duplicate-start
/// HSPs, then sorts by (context-query.end, subject.end) and removes duplicate-end HSPs.
/// For blastn with forward-normalised Rust coordinates the context-query transforms are:
///   Plus  strand: context-start = (q_start, s_start); context-end = (q_end, s_end)
///   Minus strand: context-start = (q_end,   s_start); context-end = (q_start, s_end)
///
/// Input `results` must be sorted by score descending so the first member of each
/// duplicate group is the highest-scoring one (which we keep).
#[cfg(test)] // retained only as a reference implementation exercised by the unit tests
fn purge_common_endpoints(results: &mut Vec<AlignResult>, _s_len: u32) {
    let n = results.len();
    if n < 2 { return; }

    let mut keep = vec![true; n];

    // ── pass 1: common start point ──────────────────────────────────────────
    {
        // Sort indices by (strand, start_key, score-desc) so that within each
        // duplicate group the first (kept) element has the highest score.
        // start_key: (q_start, s_start) for plus; (q_end, s_start) for minus.
        let mut idx: Vec<usize> = (0..n).collect();
        idx.sort_unstable_by(|&a, &b| {
            let ha = &results[a].hsp;
            let hb = &results[b].hsp;
            let strand_a = ha.strand == Strand::Minus;
            let strand_b = hb.strand == Strand::Minus;
            if strand_a != strand_b { return strand_a.cmp(&strand_b); }
            let key_qa = match ha.strand { Strand::Plus => ha.q_start, Strand::Minus => ha.q_end };
            let key_qb = match hb.strand { Strand::Plus => hb.q_start, Strand::Minus => hb.q_end };
            if key_qa != key_qb { return key_qa.cmp(&key_qb); }
            if ha.s_start != hb.s_start { return ha.s_start.cmp(&hb.s_start); }
            // tiebreakers within group: score DESC, then end offsets DESC
            if ha.score != hb.score { return hb.score.cmp(&ha.score); }
            let end_qa = match ha.strand { Strand::Plus => ha.q_end, Strand::Minus => ha.q_start };
            let end_qb = match hb.strand { Strand::Plus => hb.q_end, Strand::Minus => hb.q_start };
            if end_qa != end_qb { return end_qb.cmp(&end_qa); }
            hb.s_end.cmp(&ha.s_end)
        });

        let mut i = 0;
        while i < n {
            if !keep[idx[i]] { i += 1; continue; }
            let h0 = &results[idx[i]].hsp;
            let key0_q = match h0.strand {
                Strand::Plus  => h0.q_start,
                Strand::Minus => h0.q_end,
            };
            let key0_s  = h0.s_start;
            let strand0 = h0.strand;
            let mut j = i + 1;
            while j < n {
                let hj = &results[idx[j]].hsp;
                if hj.strand != strand0 { break; }
                let key_q = match hj.strand {
                    Strand::Plus  => hj.q_start,
                    Strand::Minus => hj.q_end,
                };
                if key_q != key0_q || hj.s_start != key0_s { break; }
                keep[idx[j]] = false;
                j += 1;
            }
            i = j;
        }
    }

    // ── pass 2: common end point ─────────────────────────────────────────────
    {
        // end_key: (q_end, s_end) for plus; (q_start, s_end) for minus.
        // Mirrors s_QueryEndCompareHSPs: (context, query.end ASC, subject.end ASC,
        //   score DESC, query.offset DESC, subject.offset DESC).
        let mut idx: Vec<usize> = (0..n).collect();
        idx.sort_unstable_by(|&a, &b| {
            let ha = &results[a].hsp;
            let hb = &results[b].hsp;
            let strand_a = ha.strand == Strand::Minus;
            let strand_b = hb.strand == Strand::Minus;
            if strand_a != strand_b { return strand_a.cmp(&strand_b); }
            let key_qa = match ha.strand { Strand::Plus => ha.q_end, Strand::Minus => ha.q_start };
            let key_qb = match hb.strand { Strand::Plus => hb.q_end, Strand::Minus => hb.q_start };
            if key_qa != key_qb { return key_qa.cmp(&key_qb); }
            if ha.s_end != hb.s_end { return ha.s_end.cmp(&hb.s_end); }
            // tiebreakers within group: score DESC, then start offsets DESC
            if ha.score != hb.score { return hb.score.cmp(&ha.score); }
            let off_qa = match ha.strand { Strand::Plus => ha.q_start, Strand::Minus => ha.q_end };
            let off_qb = match hb.strand { Strand::Plus => hb.q_start, Strand::Minus => hb.q_end };
            if off_qa != off_qb { return off_qb.cmp(&off_qa); }
            hb.s_start.cmp(&ha.s_start)
        });

        let mut i = 0;
        while i < n {
            if !keep[idx[i]] { i += 1; continue; }
            let h0 = &results[idx[i]].hsp;
            let key0_q = match h0.strand {
                Strand::Plus  => h0.q_end,
                Strand::Minus => h0.q_start,
            };
            let key0_s  = h0.s_end;
            let strand0 = h0.strand;
            let mut j = i + 1;
            while j < n {
                let hj = &results[idx[j]].hsp;
                if hj.strand != strand0 { break; }
                let key_q = match hj.strand {
                    Strand::Plus  => hj.q_end,
                    Strand::Minus => hj.q_start,
                };
                if key_q != key0_q || hj.s_end != key0_s { break; }
                keep[idx[j]] = false;
                j += 1;
            }
            i = j;
        }
    }

    let mut j = 0;
    results.retain(|_| { let k = keep[j]; j += 1; k });
}

/// Bidirectional gapped alignment anchored at (q_seed, s_seed).
///
/// Matches BLAST_GappedAlignmentWithTraceback convention exactly:
///   LEFT extension (reverse DP, `reverse=true`):  includes the seed.
///   RIGHT extension (forward DP, `reverse=false`): excludes the seed.
///
/// Callers are responsible for presenting sequences and seeds in the correct
/// orientation so that LEFT is always the extension that includes the seed,
/// exactly as in NCBI's implementation.
///
/// Returns (score, q_start, q_end, s_start, s_end, combined_edit_script),
/// all positions 0-indexed into the real (non-sentinel) bases.
fn gapped_extend_bidirectional(
    query: &[u8],
    subject: &[u8],
    q_seed: u32,
    s_seed: u32,
    gap_open: i32,
    gap_extend: i32,
    xdrop: i32,
    matrix: &ScoreMatrix,
    ws: &mut AlignWorkspace,
) -> Option<(i32, u32, u32, u32, u32, EditScript)> {
    let qa = &query[1..query.len() - 1];
    let sa = &subject[1..subject.len() - 1];

    let left_m = q_seed as usize + 1;
    let left_n = s_seed as usize + 1;
    let left = align_ex(
        &qa[..left_m], &sa[..left_n],
        left_m, left_n,
        gap_open, gap_extend, xdrop, matrix,
        true, ws,
    );
    // NCBI keeps an empty left half (BLAST_GappedAlignmentWithTraceback: sl=0,
    // alignment starts at seed+1).  This happens when improve_seed snaps the
    // seed into an N-x-N identity run whose columns all score negatively —
    // dropping the HSP here loses hits NCBI reports (bug #47).
    let q_start = q_seed + 1 - left.a_len as u32;
    let s_start = s_seed + 1 - left.b_len as u32;

    let rq = q_seed as usize + 1;
    let rs = s_seed as usize + 1;
    let right_m = qa.len().saturating_sub(rq);
    let right_n = sa.len().saturating_sub(rs);
    let right = if right_m > 0 && right_n > 0 {
        align_ex(
            &qa[rq..], &sa[rs..],
            right_m, right_n,
            gap_open, gap_extend, xdrop, matrix,
            false, ws,
        )
    } else {
        GapAlignResult { score: 0, a_len: 0, b_len: 0, edit_script: EditScript::new() }
    };

    let mut total_score = left.score + right.score;

    let mut q_start = q_start;
    let mut s_start = s_start;
    let mut q_end = q_seed + 1 + right.a_len as u32;
    let mut s_end = s_seed + 1 + right.b_len as u32;

    let mut edit_script = left.edit_script;
    edit_script.reverse(); // flip left-extension to forward direction
    for (op, count) in right.edit_script.ops {
        edit_script.push(op, count);
    }

    // NCBI BLAST_GappedAlignmentWithTraceback (blast_gapalign.c): "rarely ... it
    // is possible to compute an optimal alignment with a leading or trailing
    // gap.  Prune these unneeded gaps here and update the score and alignment
    // boundaries."  This fires when the forced seed cell sits in an N-run and
    // the adjacent path enters/leaves it through a gap op.
    while let Some(&(op, n)) = edit_script.ops.first() {
        if op == EditOp::Sub { break; }
        total_score += gap_open + n as i32 * gap_extend;
        match op {
            EditOp::GapInQuery => s_start += n,
            _ => q_start += n,
        }
        edit_script.ops.remove(0);
    }
    while let Some(&(op, n)) = edit_script.ops.last() {
        if op == EditOp::Sub { break; }
        total_score += gap_open + n as i32 * gap_extend;
        match op {
            EditOp::GapInQuery => s_end -= n,
            _ => q_end -= n,
        }
        edit_script.ops.pop();
    }

    if total_score < 0 {
        return None;
    }

    Some((total_score, q_start, q_end, s_start, s_end, edit_script))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::AlignStats;
    use crate::search::gapped::AlignWorkspace;

    // ── helpers shared by ncbi-ported tests ──────────────────────────────────

    fn parse_fasta_seq(bytes: &[u8]) -> Vec<u8> {
        bytes.split(|&b| b == b'\n')
             .filter(|l| !l.is_empty() && !l.starts_with(b">"))
             .flat_map(|l| l.iter().filter(|&&b| b != b'\r').copied())
             .collect()
    }

    // Convert raw IUPAC DNA to BLASTNA with leading and trailing sentinel (15).
    fn dna_to_blastna_sentinels(dna: &[u8]) -> Vec<u8> {
        let mut v = Vec::with_capacity(dna.len() + 2);
        v.push(15u8);
        for &b in dna {
            v.push(match b {
                b'A' | b'a' => 0, b'C' | b'c' => 1,
                b'G' | b'g' => 2, b'T' | b't' => 3,
                _ => 14,
            });
        }
        v.push(15u8);
        v
    }

    // +1/-3 scoring matrix (BLASTN defaults).
    fn blastn_1m3_matrix() -> ScoreMatrix {
        let mut scores = [[0i32; 16]; 16];
        for i in 0..4usize {
            for j in 0..4usize {
                scores[i][j] = if i == j { 1 } else { -3 };
            }
        }
        ScoreMatrix { scores, freqs: [0.0f64; 16], lambda: 0.0, name: "blastn_1m3".to_string(), karlin: None }
    }

    fn make_prelim(score: i32, q_start: u32, q_end: u32, s_start: u32, s_end: u32) -> PrelimHsp {
        PrelimHsp { q_seed: q_start, s_seed: s_start, q_start, q_end, s_start, s_end, score, strand: Strand::Plus }
    }

    fn make_result(score: i32, q_start: u32, q_end: u32, s_start: u32, s_end: u32) -> AlignResult {
        AlignResult {
            hsp: Hsp {
                score, q_start, q_end, q_len: q_end,
                s_start, s_end, s_len: s_end,
                strand: Strand::Plus,
                edit_script: EditScript::new(),
                q_seq: Vec::new(),
                s_seq: Vec::new(),
            },
            query_id: String::new(),
            subject_id: String::new(),
            stats: AlignStats::default(),
        }
    }

    // Port of NCBI's testCheckHSPCommonEndpoints (blasthits_unit_test.cpp).
    // Verifies that purge_prelim_common_endpoints removes HSPs sharing a
    // query-start or query-end point with a higher-scoring HSP.
    // ── bug #47 regression tests ─────────────────────────────────────────────
    // NCBI's BlastGetStartForGappedAlignmentNucl can snap the traceback seed
    // into an N-x-N identity run (blastna N==N counts as a byte-equality
    // match).  Two behaviors around the forced seed must match NCBI:
    // an empty left half keeps the HSP, and leading/trailing gap ops are
    // pruned with their cost refunded (BLAST_GappedAlignmentWithTraceback).

    fn n_penalty_matrix(match_score: i32, mismatch: i32) -> ScoreMatrix {
        let mut scores = [[0i32; 16]; 16];
        for i in 0..16usize {
            for j in 0..16usize {
                scores[i][j] = if i < 4 && j < 4 {
                    if i == j { match_score } else { mismatch }
                } else {
                    -1 // N rows/cols: -1, like the RM nt matrix
                };
            }
        }
        ScoreMatrix { scores, freqs: [0.0f64; 16], lambda: 0.0, name: "n_pen".to_string(), karlin: None }
    }

    #[test]
    fn test_bug47_empty_left_half_keeps_hsp() {
        // Seed placed mid N-run: every column left of the seed scores -1, so
        // the left extension is empty.  NCBI keeps the HSP (alignment starts
        // at seed+1); the old code returned None and lost the hit.
        let seq = dna_to_blastna_sentinels(b"NNNNNNNNNNNNAAAAAAAAAAAAAAAAAAAA");
        let matrix = n_penalty_matrix(1, -3);
        let mut ws = AlignWorkspace::new();
        let r = gapped_extend_bidirectional(&seq, &seq, 6, 6, 20, 5, 250, &matrix, &mut ws);
        let (score, q_start, q_end, s_start, s_end, es) =
            r.expect("HSP with empty left half must be kept (NCBI sl=0 case)");
        // right half: 5 N,N columns (-1 each) then 20 A,A (+1 each) = 15
        assert_eq!(score, 15);
        assert_eq!((q_start, q_end, s_start, s_end), (7, 32, 7, 32));
        assert_eq!(es.ops, vec![(EditOp::Sub, 25)]);
    }

    #[test]
    fn test_bug47_trailing_gap_pruned() {
        // Left half's optimal path leaves the forced seed corner through a
        // 1-base subject gap (mismatch at the corner costs more than a gap).
        // NCBI prunes the trailing gap op, refunds gap_open + gap_extend, and
        // pulls in the subject boundary.
        let qa = dna_to_blastna_sentinels(b"AAAAAAAAAAAAAAAAAAAAA"); // 21 A
        let sa = dna_to_blastna_sentinels(b"AAAAAAAAAAAAAAAAAAAAT"); // 20 A + T
        let mut matrix = n_penalty_matrix(9, -3);
        matrix.scores[0][3] = -40; // A vs T worse than a gap (open 20 + ext 5)
        matrix.scores[3][0] = -40;
        let mut ws = AlignWorkspace::new();
        let r = gapped_extend_bidirectional(&qa, &sa, 20, 20, 20, 5, 250, &matrix, &mut ws);
        let (score, q_start, q_end, s_start, s_end, es) =
            r.expect("alignment expected");
        // Unpruned: 20 A,A matches (+180) plus trailing gap (-25) ending at
        // s_end 21.  After NCBI-style pruning: score 180, s_end 20.
        assert_eq!(score, 180);
        assert_eq!((q_start, q_end, s_start, s_end), (1, 21, 0, 20));
        assert_eq!(es.ops, vec![(EditOp::Sub, 20)]);
    }

    #[test]
    fn test_purge_prelim_common_endpoints_ncbi() {
        let scores: [i32; 9] =    [1044, 995, 965, 219, 160, 125, 110, 107, 103];
        let q_starts: [u32; 9] =  [2, 2, 2, 236, 88, 259, 278, 259, 278];
        let q_ends: [u32; 9] =    [322, 336, 300, 322, 182, 322, 341, 341, 341];
        let s_starts: [u32; 9] =  [7, 7, 7, 194, 2, 194, 197, 194, 197];
        let s_ends: [u32; 9] =    [292, 293, 301, 292, 96, 292, 260, 260, 266];

        let mut hsps: Vec<PrelimHsp> = (0..9)
            .map(|i| make_prelim(scores[i], q_starts[i], q_ends[i], s_starts[i], s_ends[i]))
            .collect();

        purge_prelim_common_endpoints(&mut hsps, 400);

        // NCBI kSurvivingIndices = {4, 0, 6}: scores {1044, 160, 110}.
        assert_eq!(hsps.len(), 3);
        let mut got: Vec<(i32, u32, u32, u32, u32)> = hsps.iter()
            .map(|h| (h.score, h.q_start, h.q_end, h.s_start, h.s_end))
            .collect();
        got.sort_unstable_by(|a, b| b.0.cmp(&a.0));
        assert_eq!(got[0], (1044, 2,   322, 7,   292));
        assert_eq!(got[1], (160,  88,  182, 2,   96));
        assert_eq!(got[2], (110,  278, 341, 197, 260));
    }

    // Port of NCBI's testCheckHSPCommonEndpoints applied to AlignResult objects.
    // Verifies purge_common_endpoints with the same test data as above.
    #[test]
    fn test_purge_common_endpoints_ncbi() {
        let scores: [i32; 9] =    [1044, 995, 965, 219, 160, 125, 110, 107, 103];
        let q_starts: [u32; 9] =  [2, 2, 2, 236, 88, 259, 278, 259, 278];
        let q_ends: [u32; 9] =    [322, 336, 300, 322, 182, 322, 341, 341, 341];
        let s_starts: [u32; 9] =  [7, 7, 7, 194, 2, 194, 197, 194, 197];
        let s_ends: [u32; 9] =    [292, 293, 301, 292, 96, 292, 260, 260, 266];

        let mut results: Vec<AlignResult> = (0..9)
            .map(|i| make_result(scores[i], q_starts[i], q_ends[i], s_starts[i], s_ends[i]))
            .collect();

        purge_common_endpoints(&mut results, 400);

        // NCBI kSurvivingIndices = {4, 0, 6}: scores {1044, 160, 110}.
        assert_eq!(results.len(), 3);
        let mut got: Vec<(i32, u32, u32, u32, u32)> = results.iter()
            .map(|r| (r.hsp.score, r.hsp.q_start, r.hsp.q_end, r.hsp.s_start, r.hsp.s_end))
            .collect();
        got.sort_unstable_by(|a, b| b.0.cmp(&a.0));
        assert_eq!(got[0], (1044, 2,   322, 7,   292));
        assert_eq!(got[1], (160,  88,  182, 2,   96));
        assert_eq!(got[2], (110,  278, 341, 197, 260));
    }

    // Port of NCBI's testHSPResultsApplyMasklevel (rmblast_blasthits_unit_test.cpp).
    // Verifies apply_mask_level removes HSPs with >mask_level% query overlap
    // with a higher-scoring HSP.
    // Input: 26 HSPs (9×score=93 interleaved with 17×score=85).
    // After mask_level=80: 17 remain; after mask_level=25 on those: 9 remain.
    #[test]
    fn test_apply_mask_level_ncbi() {
        let s_offsets: [u32; 26] =
            [10,382,14,83,4,382,1000,203,54,32,64,382,89,183,813,132,14,1344,321,224,34,8341,1344,254,861,834];
        let q_offsets: [u32; 26] =
            [1,4,14,21,24,34,41,44,54,61,64,74,81,84,94,101,104,114,121,124,134,141,144,154,161,164];
        let lengths: [u32; 26] =
            [17,9,9,17,9,9,17,9,9,17,9,9,17,9,9,17,9,9,17,9,9,17,9,9,17,9];
        let scores: [i32; 26] =
            [93,85,85,93,85,85,93,85,85,93,85,85,93,85,85,93,85,85,93,85,85,93,85,85,93,85];

        let mut results: Vec<AlignResult> = (0..26)
            .map(|i| make_result(
                scores[i],
                q_offsets[i],
                q_offsets[i] + lengths[i],
                s_offsets[i],
                s_offsets[i] + lengths[i],
            ))
            .collect();

        // masklevel 80: 17 remain (per NCBI comment: q_starts 1,14,21,34,41,54,61,74,81,94,101,114,121,134,141,154,161)
        apply_mask_level(&mut results, 80, &[]);
        assert_eq!(results.len(), 17);

        // masklevel 25 applied to those 17: 9 remain (per NCBI comment: q_starts 1,21,41,61,81,101,121,141,161)
        apply_mask_level(&mut results, 25, &[]);
        assert_eq!(results.len(), 9);

        // All 9 survivors are the score=93 HSPs.
        let q_start_set: std::collections::HashSet<u32> =
            results.iter().map(|r| r.hsp.q_start).collect();
        let expected: std::collections::HashSet<u32> =
            [1u32, 21, 41, 61, 81, 101, 121, 141, 161].iter().copied().collect();
        assert_eq!(q_start_set, expected);
    }

    // Port of NCBI's nuclwordfinder_unit_test.cpp (testExtend).
    //
    // Runs Phase 1 (scan + exact extension + ungapped DP) against real sequences:
    //   query  = gi|3090   X65526.1  N.patriciarum xylanase A  (2338 bp)
    //   subject = gi|33383640 AY131336.1 Neocallimastix xylanase xyn11B (1011 bp)
    //
    // Parameters from the NCBI test: word_size=19, lut_word_length=8 (scan_step=12),
    // xdrop=11, ungapped cutoff=14 (→ min_raw_gapped_score=28), +1/-3 scoring.
    //
    // Expected 5 ungapped hits (from checkResults in the NCBI test):
    //   q_start  s_start  length  score
    //      233       0      31     31
    //     1037      66     945    853
    //      263      27     685    541
    //     1911     811     101     73
    //     1782     940      63     51
    #[test]
    fn test_nuclwordfinder_phase1_ncbi() {
        let query_fa = include_bytes!("../../tests/data/gi3090_xylanase_query.fa");
        let subj_fa  = include_bytes!("../../tests/data/gi33383640_xylanase_subject.fa");

        let query_dna = parse_fasta_seq(query_fa);
        let subj_dna  = parse_fasta_seq(subj_fa);

        let query   = dna_to_blastna_sentinels(&query_dna);
        let subject = dna_to_blastna_sentinels(&subj_dna);
        let q_len   = query_dna.len() as u32;
        let s_len   = subj_dna.len() as u32;

        let packed = {
            let n = subject.len().saturating_sub(2);
            let mut blk = seqblk_from_blastna(&subject[1..1 + n]);
            blast_compress_blastna_sequence(&mut blk);
            blk.packed
        };

        let query_bases = &query[1..1 + q_len as usize];
        let lut = NaLookup::Small(blast_na_lookup_table_new(
            query_bases,
            &[(0, q_len as i32 - 1)],
            19, // word_length
            8,  // lut_word_length
        ));

        let params = SearchParams {
            word_size: 19,
            xdrop_ungap: 11,
            min_raw_gapped_score: 28, // ungapped_cutoff = 28/2 = 14
            dust: false,
            ..Default::default()
        };
        let matrix = blastn_1m3_matrix();

        let mut ungapped = Vec::new();
        collect_ungapped(
            &query, &query, q_len, &subject, &packed[3..], s_len,
            Strand::Plus, &lut, &params, &matrix, &mut DiscardSeeds, &mut ungapped,
        );

        let exp_q: [u32; 5] = [233, 1037, 263, 1911, 1782];
        let exp_s: [u32; 5] = [0, 66, 27, 811, 940];
        let exp_len: [u32; 5] = [31, 945, 685, 101, 63];
        let exp_score: [i32; 5] = [31, 853, 541, 73, 51];

        assert_eq!(ungapped.len(), 5,
            "expected 5 ungapped hits, got {}; hits: {:?}",
            ungapped.len(),
            ungapped.iter().map(|u| (u.q_start, u.s_start, u.q_end - u.q_start, u.score)).collect::<Vec<_>>());

        let mut got: Vec<(u32, u32, u32, i32)> = ungapped.iter()
            .map(|u| (u.q_start, u.s_start, u.q_end - u.q_start, u.score))
            .collect();
        got.sort_unstable_by_key(|t| t.0);

        let mut exp: Vec<(u32, u32, u32, i32)> = (0..5)
            .map(|i| (exp_q[i], exp_s[i], exp_len[i], exp_score[i]))
            .collect();
        exp.sort_unstable_by_key(|t| t.0);

        for (i, (g, e)) in got.iter().zip(exp.iter()).enumerate() {
            assert_eq!(*g, *e,
                "hit {i}: got (q_start={} s_start={} len={} score={}) \
                 expected (q_start={} s_start={} len={} score={})",
                g.0, g.1, g.2, g.3, e.0, e.1, e.2, e.3);
        }
    }

    // Port of NCBI's blastextend_unit_test.cpp (testGapAlignment), ALIGN_EX path.
    //
    // Sequences:
    //   query   = gi|2655203 AF026469.1 range 20001-35000 (15000 bp, mouse Unp region)
    //   subject = gi|2516238 AB004664.1 (3873 bp, mouse Rab33B)
    //
    // Parameters (BLASTN defaults): reward=1, penalty=-3, gap_open=5, gap_extend=2,
    // xdrop_gap=30.
    //
    // Note: NCBI's BLAST_GetGappedScore uses s_BlastAlignPackedNucl for the
    // score-only pass (4-base boundary alignment) before calling ALIGN_EX for
    // traceback.  Our gapped_extend_bidirectional uses ALIGN_EX throughout
    // (matching BLAST_GappedAlignmentWithTraceback).  For long, high-identity
    // alignments the packed-DP and ALIGN_EX boundaries agree; for seeds at the
    // edge of a homologous region they may differ.  We test seed 7 (20-base
    // alignment) where the algorithms agree.
    //
    //   Seed 7: q_off=3212 s_off=3640, ungapped [3201..3221) [3629..3649) score=16
    //           → expected q_start=3201 q_end=3221 s_start=3629 s_end=3649
    #[test]
    fn test_blastextend_gapped_ncbi() {
        let query_fa = include_bytes!("../../tests/data/gi2655203_unp_20001_35000.fa");
        let subj_fa  = include_bytes!("../../tests/data/gi2516238_rab33b.fa");

        let query_dna = parse_fasta_seq(query_fa);
        let subj_dna  = parse_fasta_seq(subj_fa);

        let query   = dna_to_blastna_sentinels(&query_dna);
        let subject = dna_to_blastna_sentinels(&subj_dna);

        // Slices without sentinels — for improve_seed.
        let qa = &query[1..query.len() - 1];
        let sa = &subject[1..subject.len() - 1];

        let matrix = blastn_1m3_matrix();
        let gap_open   = 5i32;
        let gap_extend = 2i32;
        let xdrop      = 30i32;
        let mut ws     = AlignWorkspace::new();

        // Seed 7: q_off=3212, s_off=3640; ungapped q_start=3201 length=20 s_start=3629
        {
            let (new_q, new_s) = improve_seed(
                qa, sa, 3212, 3640,
                3201, 3201 + 20, 3629, 3629 + 20, false,
            );
            let result = gapped_extend_bidirectional(
                &query, &subject, new_q, new_s,
                gap_open, gap_extend, xdrop, &matrix, &mut ws,
            );
            let (_, q_start, q_end, s_start, s_end, _) =
                result.expect("seed (3212,3640) should produce a gapped alignment");
            assert_eq!(q_start, 3201, "seed7 q_start");
            assert_eq!(q_end,   3221, "seed7 q_end");
            assert_eq!(s_start, 3629, "seed7 s_start");
            assert_eq!(s_end,   3649, "seed7 s_end");
        }
    }
}

#[cfg(test)]
mod hit_list_order_tests {
    use super::sort_hit_list_order;
    use crate::hits::{EditScript, Hsp, Strand};
    use crate::output::AlignResult;
    use crate::stats::AlignStats;

    fn hit(subject: &str, score: i32, q_start: u32, q_end: u32, s_start: u32, s_end: u32, strand: Strand) -> AlignResult {
        AlignResult {
            hsp: Hsp {
                score, q_start, q_end, q_len: 1000, s_start, s_end, s_len: 500, strand,
                edit_script: EditScript::new(), q_seq: Vec::new(), s_seq: Vec::new(),
            },
            query_id: "q".to_string(),
            subject_id: subject.to_string(),
            stats: AlignStats::default(),
        }
    }

    fn key(r: &AlignResult) -> (String, i32, u32) {
        (r.subject_id.clone(), r.hsp.score, r.hsp.q_start)
    }

    /// Subjects come out grouped and ranked by their top score, tied top
    /// scores by higher oid first, and each subject's HSPs in
    /// ScoreCompareHSPs order (NCBI measures the minus-strand query offset in
    /// reverse-complement coordinates).
    #[test]
    fn groups_by_subject_and_ranks_like_ncbi() {
        let names: Vec<String> = ["s0", "s1", "s2"].iter().map(|s| s.to_string()).collect();
        let mut results = vec![
            hit("s1", 300, 10, 50, 0, 40, Strand::Plus),
            hit("s0", 900, 100, 200, 0, 100, Strand::Plus),
            hit("s2", 900, 300, 400, 0, 100, Strand::Plus),
            hit("s1", 700, 500, 600, 0, 100, Strand::Plus),
            hit("s0", 100, 700, 720, 0, 20, Strand::Plus),
            // Same score and subject span: minus-strand q_off = 1000-960 = 40,
            // which sorts before the plus-strand q_off of 60.
            hit("s2", 500, 60, 80, 0, 20, Strand::Plus),
            hit("s2", 500, 940, 960, 0, 20, Strand::Minus),
        ];
        sort_hit_list_order(&mut results, &names);
        let got: Vec<_> = results.iter().map(key).collect();
        let want = vec![
            ("s2".to_string(), 900, 300),
            ("s2".to_string(), 500, 940),
            ("s2".to_string(), 500, 60),
            ("s0".to_string(), 900, 100),
            ("s0".to_string(), 100, 700),
            ("s1".to_string(), 700, 500),
            ("s1".to_string(), 300, 10),
        ];
        assert_eq!(got, want);
    }
}
