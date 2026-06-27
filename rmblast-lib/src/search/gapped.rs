//! Gapped alignment with full traceback — faithful Rust port of ALIGN_EX.
//!
//! Reference: blast_gapalign.c line 374.
//!
//! Sequence conventions:
//!   `a` = query BLASTNA bases, 0-indexed, length m.
//!   `b` = subject BLASTNA bases, 0-indexed, length n.
//!
//! NCBI uses 1-indexed A[1..M] and B[1..N] with sentinels at index 0.
//! We use 0-indexed slices and adjust the pointer arithmetic accordingly.
//!
//! NCBI edit-script byte layout (EGapAlignOpType + extend bits):
//!   SCRIPT_GAP_IN_A   = eGapAlignDel = 0  (gap in query: subject advances)
//!   SCRIPT_SUB        = eGapAlignSub = 3  (aligned pair)
//!   SCRIPT_GAP_IN_B   = eGapAlignIns = 6  (gap in subject: query advances)
//!   SCRIPT_OP_MASK    = 0x07
//!   SCRIPT_EXTEND_GAP_A = 0x10
//!   SCRIPT_EXTEND_GAP_B = 0x40

use std::sync::atomic::{AtomicU64, Ordering};
use crate::hits::{EditOp, EditScript};
use crate::matrix::ScoreMatrix;

pub static TOTAL_DP_CELLS: AtomicU64 = AtomicU64::new(0);

const SCRIPT_GAP_IN_A: u8 = 0;
const SCRIPT_SUB: u8 = 3;
const SCRIPT_GAP_IN_B: u8 = 6;
const SCRIPT_OP_MASK: u8 = 0x07;
const SCRIPT_EXTEND_GAP_A: u8 = 0x10;
const SCRIPT_EXTEND_GAP_B: u8 = 0x40;

// MININT = i32::MIN/2.  Adding matrix scores to MININT never overflows
// (matrix values are bounded ±100, MININT ≈ -1.07e9).
const MININT: i32 = i32::MIN / 2;

#[derive(Clone, Copy, Default)]
pub struct DpCell {
    best: i32,
    best_gap: i32,
}

/// Result of a gapped alignment.
#[derive(Debug, Clone)]
pub struct GapAlignResult {
    pub score: i32,
    /// Number of query bases consumed (a[0..a_len-1] were aligned).
    pub a_len: usize,
    /// Number of subject bases consumed (b[0..b_len-1] were aligned).
    pub b_len: usize,
    /// Run-length edit script for the alignment.
    pub edit_script: EditScript,
}

/// Reusable workspace for `align_ex` — eliminates per-call allocations.
/// Create once per thread / strand search and pass to every `align_ex` call.
pub struct AlignWorkspace {
    /// DP score array, reused across calls (re-initialized at start of each call).
    dp: Vec<DpCell>,
    /// Flat traceback byte stream, reset (len→0) between calls.
    flat_edit: Vec<u8>,
    /// (byte_offset_in_flat, first_b_of_row) for each DP row.
    row_info: Vec<(usize, usize)>,
    /// Separate dp array for score-only first pass (not shared with traceback dp).
    pub dp_score: Vec<DpCell>,
}

impl AlignWorkspace {
    pub fn new() -> Self {
        Self {
            dp: Vec::new(),
            flat_edit: Vec::new(),
            row_info: Vec::new(),
            dp_score: Vec::new(),
        }
    }

}

/// Score-only gapped alignment — no traceback storage (mirrors score_only=TRUE path).
///
/// Returns (score, a_len, b_len). Much faster than `align_ex` when traceback not needed.
/// Score-only gapped alignment (dispatch wrapper — picks monomorphised variant).
pub fn align_ex_score_only(
    a: &[u8],
    b: &[u8],
    m: usize,
    n: usize,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    matrix: &ScoreMatrix,
    reverse: bool,
    dp: &mut Vec<DpCell>,
    dump: bool,
) -> (i32, usize, usize) {
    if reverse {
        align_ex_score_only_inner::<true>(a, b, m, n, gap_open, gap_extend, x_dropoff, matrix, dp, dump)
    } else {
        align_ex_score_only_inner::<false>(a, b, m, n, gap_open, gap_extend, x_dropoff, matrix, dp, dump)
    }
}

#[inline(never)]
fn align_ex_score_only_inner<const REVERSE: bool>(
    a: &[u8],
    b: &[u8],
    m: usize,
    n: usize,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    matrix: &ScoreMatrix,
    dp: &mut Vec<DpCell>,
    dump: bool,
) -> (i32, usize, usize) {
    let m = m.min(a.len());
    let n = n.min(b.len());
    if m == 0 || n == 0 {
        return (0, 0, 0);
    }
    let gap_open_extend = gap_open + gap_extend;
    let xdrop = x_dropoff.max(gap_open_extend);
    let num_extra = if gap_extend > 0 { (xdrop / gap_extend + 3) as usize } else { n + 3 };
    // Mirror NCBI s_BlastAlignPackedNucl's dynamic DP window (blast_gapalign.c:3109-3257):
    // the band is indexed absolutely by b, and b_size can drift up to n as the alignment
    // accumulates net indels — so the array must be able to reach n+2 cells.  NCBI starts
    // small (num_extra+100) and reallocs by doubling as the band drifts, rather than
    // pre-allocating n+2 (which would waste ~8 MB per workspace for genome-length inner
    // dimensions).  The previous FIXED cap `m + 3*num_extra + 10` under-sized the band for
    // gappy alignments (e.g. the TIFr satellite, whose query band drifts 72 cells past the
    // subject length) → the band was truncated → suboptimal extension on asymmetric
    // matrices (bug #36 under-extension sub-case).
    let max_cap = n + 2;
    let init_cap = (num_extra + 102).min(max_cap);
    if dp.len() < init_cap {
        dp.resize(init_cap, DpCell { best: MININT, best_gap: MININT });
    }

    dp[0].best = 0;
    dp[0].best_gap = -gap_open_extend;
    let mut score = -gap_open_extend;
    let mut b_size = 1usize;
    for j in 1..=n {
        if score < -xdrop || j >= dp.len() { break; }
        dp[j].best = score;
        dp[j].best_gap = score - gap_open_extend;
        score -= gap_extend;
        b_size = j + 1;
    }

    let mut best_score = 0i32;
    let mut best_a = 0usize;
    let mut best_b = 0usize;
    let mut first_b_index = 0usize;
    let mut total_cells_so: u64 = 0;

    for a_idx in 0..m {
        let a_phys = if REVERSE { m - 1 - a_idx } else { a_idx };
        // SAFETY: a_phys < m <= a.len(); BLASTNA values are 0-14, matrix has 16 rows.
        let matrix_row = unsafe {
            let a_base = (*a.get_unchecked(a_phys) & 15) as usize;
            matrix.scores.get_unchecked(a_base)
        };
        let mut score = MININT;
        let mut score_gap_row = MININT;
        let mut last_b_index = first_b_index;
        let prev_best = best_score;

        total_cells_so += (b_size - first_b_index) as u64;

        // Split at the sequence boundary: the sentinel case (b_idx == n, b_base = 0)
        // runs at most once and is handled separately below.  This removes the per-cell
        // `b_idx >= n` branch from the hot loop.
        let inner_end = b_size.min(n);

        // Explicit b pointer increments (FORWARD: +1, REVERSE: -1) so LLVM can keep it in
        // a register without needing to load b.as_ptr() from the stack each cell.
        // black_box on FORWARD prevents LLVM from decomposing b_cur into base+index form
        // (b.as_ptr() + b_idx), which would force it onto the stack.  REVERSE already
        // produces a register pointer naturally due to the decrement direction.
        // SAFETY: first_b_index <= b_idx < inner_end <= n <= b.len(); BLASTNA values 0-14.
        let mut b_cur: *const u8 = if REVERSE {
            unsafe { b.as_ptr().add(n - 1 - first_b_index) }
        } else {
            // black_box makes this pointer opaque to LLVM, preventing the base+index
            // rewrite that would otherwise spill b.as_ptr() to 0x80(%rsp) each cell.
            std::hint::black_box(unsafe { b.as_ptr().add(first_b_index) })
        };

        // Hot inner loop — no sentinel branch.
        let mut b_idx = first_b_index;
        while b_idx < inner_end {
            let b_base = unsafe { *b_cur as usize };
            b_cur = if REVERSE { unsafe { b_cur.sub(1) } } else { unsafe { b_cur.add(1) } };

            let (score_gap_col, next_score) = unsafe {
                let cell = dp.get_unchecked(b_idx);
                (cell.best_gap, cell.best + *matrix_row.get_unchecked(b_base))
            };

            let mut cell_score = score;
            if score_gap_col > cell_score { cell_score = score_gap_col; }
            if score_gap_row > cell_score { cell_score = score_gap_row; }

            if best_score - cell_score > xdrop {
                if b_idx == first_b_index {
                    first_b_index += 1;
                } else {
                    unsafe { dp.get_unchecked_mut(b_idx).best = MININT; }
                }
            } else {
                last_b_index = b_idx;
                if cell_score > best_score { best_b = b_idx; best_score = cell_score; }
                score_gap_row -= gap_extend;
                let open_gap = cell_score - gap_open_extend;
                let score_gap_col_ext = score_gap_col - gap_extend;
                unsafe {
                    let cell = dp.get_unchecked_mut(b_idx);
                    cell.best_gap = if score_gap_col_ext > open_gap { score_gap_col_ext } else { open_gap };
                    cell.best = cell_score;
                }
                score_gap_row = score_gap_row.max(open_gap);
            }
            score = next_score;
            b_idx += 1;
        }

        // Sentinel position (b_idx == n, b_base = 0 = NULLB).
        // NCBI nucleotide BLAST: FENCE_SENTRY=201 but the sequence sentinels are NULLB=0.
        // The fence check (matrix_index==FENCE_SENTRY) NEVER fires.  NCBI processes this
        // cell normally with b_base=0, which also advances first_b_index past n when the
        // cell is x-dropped (enabling outer-loop termination via first_b_index >= b_size).
        if b_size > n {
            let b_idx = n;
            // Sentinel cell (b_base would be NULLB=0): only score_gap_col feeds the
            // recurrence; the per-cell `next_score` is never read after the last
            // column, so neither it nor b_base is needed here.
            let score_gap_col = unsafe { dp.get_unchecked(b_idx).best_gap };
            let mut cell_score = score;
            if score_gap_col > cell_score { cell_score = score_gap_col; }
            if score_gap_row > cell_score { cell_score = score_gap_row; }
            if best_score - cell_score > xdrop {
                if b_idx == first_b_index { first_b_index += 1; }
                else { unsafe { dp.get_unchecked_mut(b_idx).best = MININT; } }
            } else {
                last_b_index = b_idx;
                if cell_score > best_score { best_b = b_idx; best_score = cell_score; }
                score_gap_row -= gap_extend;
                let open_gap = cell_score - gap_open_extend;
                let score_gap_col_ext = score_gap_col - gap_extend;
                unsafe {
                    let cell = dp.get_unchecked_mut(b_idx);
                    cell.best_gap = if score_gap_col_ext > open_gap { score_gap_col_ext } else { open_gap };
                    cell.best = cell_score;
                }
                score_gap_row = score_gap_row.max(open_gap);
            }
        }

        if best_score > prev_best { best_a = a_idx + 1; }
        if dump && !REVERSE {
            eprintln!("ROW_R a={} fbi={} bsz={} prev={} best={} best_b={}",
                a_idx, first_b_index, b_size, prev_best, best_score, best_b);
        }
        if first_b_index >= b_size { break; }

        // Grow the DP window if the band is about to drift past the allocation
        // (NCBI s_BlastAlignPackedNucl realloc, blast_gapalign.c:3249-3257).  The upcoming
        // band extension can add up to ~num_extra cells plus a sentinel.
        if last_b_index + num_extra + 3 >= dp.len() && dp.len() < max_cap {
            let new_cap = (last_b_index + num_extra + 100).max(dp.len() * 2).min(max_cap);
            dp.resize(new_cap, DpCell { best: MININT, best_gap: MININT });
        }

        if last_b_index + 1 < b_size {
            b_size = last_b_index + 1;
        } else {
            while score_gap_row >= best_score - xdrop && b_size <= n && b_size < dp.len() - 1 {
                dp[b_size].best = score_gap_row;
                dp[b_size].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                b_size += 1;
            }
        }
        if b_size <= n && b_size < dp.len() {
            dp[b_size].best = MININT;
            dp[b_size].best_gap = MININT;
            b_size += 1;
        }
    }

    TOTAL_DP_CELLS.fetch_add(total_cells_so, Ordering::Relaxed);
    let reset_end = b_size.min(dp.len());
    for cell in dp[..reset_end].iter_mut() {
        *cell = DpCell { best: MININT, best_gap: MININT };
    }
    (best_score, best_a, best_b)
}

/// Score-only bidirectional gapped alignment (for preliminary extension — no traceback needed).
///
/// Mirrors `gapped_extend_bidirectional` but uses `align_ex_score_only` to skip
/// all traceback bookkeeping.  Returns `(score, q_start, q_end, s_start, s_end)`.
///
/// Loop orientation matches NCBI `s_BlastAlignPackedNucl`: subject is the OUTER loop
/// (a / M dimension) and query is the INNER loop (b / N dimension).  This ensures
/// tie-breaking in endpoint tracking agrees with NCBI when multiple cells share the
/// global-optimum score.
///
/// ORIENTATION (verified faithful 2026-06-25): NCBI's prelim is
/// `s_BlastDynProgNtGappedAlignment` -> `s_BlastAlignPackedNucl(Uint1* B, Uint1* A, ...)`
/// called as `s_BlastAlignPackedNucl(query, subject, q_length, s_length, ...)`, so
/// **B = query, A = subject**; the outer loop runs over `M = s_length` (subject) and
/// `matrix_row = matrix[A[a_index]] = matrix[subject]` — i.e. the prelim scores
/// `matrix[subject][query]`, subject-outer.  This routine matches that exactly
/// (`align_ex_score_only(sa, qa)`, a = subject).  The final/reported score instead
/// comes from the query-outer TRACEBACK (`align_ex(qa, sa)`), which is also faithful.
/// (An earlier "#41" hypothesis that the prelim should be query-outer was a MISREAD of
/// the s_BlastAlignPackedNucl arg order — there is no orientation bug.  The mirs prelim
/// score gaps on asymmetric matrices are the #38 seed-anchor difference, not orientation.)
pub fn gapped_extend_score_only(
    query: &[u8],
    subject: &[u8],
    q_seed: u32,
    s_seed: u32,
    gap_open: i32,
    gap_extend: i32,
    xdrop: i32,
    matrix: &ScoreMatrix,
    dp: &mut Vec<DpCell>,
    dump: bool,
) -> Option<(i32, u32, u32, u32, u32)> {
    let qa = &query[1..query.len() - 1];
    let sa = &subject[1..subject.len() - 1];

    // Split convention matches NCBI `s_BlastDynProgNtGappedAlignment`: the pivot
    // (q_seed, s_seed) = NCBI's (q_length, s_length) is the EXCLUSIVE left bound and
    // the FIRST base of the right extension.  (Rust previously included the pivot in
    // the LEFT extension — `..q_seed+1` — which scored the pivot base in a different
    // running/xdrop context than NCBI, diverging the prelim score on tight-xdrop
    // diverged hits while leaving the endpoints unchanged.  The span is identical
    // either way: moving the pivot from left to right shifts `lb→lb-1`, `rb→rb+1`.)
    //
    // Left extension: outer = subject (a), inner = query (b), pivot EXCLUDED.
    // Returns (score, s_ext_left, q_ext_left).
    let left_m = s_seed as usize;
    let left_n = q_seed as usize;
    let (ls, la, lb) = align_ex_score_only(
        &sa[..left_m], &qa[..left_n],
        left_m, left_n,
        gap_open, gap_extend, xdrop, matrix, true, dp, false,
    );
    let s_start = s_seed - la as u32;
    let q_start = q_seed - lb as u32;

    // Right extension: outer = subject (a), inner = query (b), pivot INCLUDED as first base.
    // Returns (score, s_ext_right, q_ext_right).
    let rq = q_seed as usize;
    let rs = s_seed as usize;
    let right_m = sa.len().saturating_sub(rs);
    let right_n = qa.len().saturating_sub(rq);
    if dump {
        eprintln!("RIGHT_EXT s_seed={} q_seed={} right_m={} right_n={} xdrop={}",
            s_seed, q_seed, right_m, right_n, xdrop);
    }
    let (rs_score, ra, rb) = if right_m > 0 && right_n > 0 {
        align_ex_score_only(
            &sa[rs..], &qa[rq..],
            right_m, right_n,
            gap_open, gap_extend, xdrop, matrix, false, dp, dump,
        )
    } else { (0, 0, 0) };

    let total_score = ls + rs_score;
    if total_score < 0 { return None; }

    let s_end = s_seed + ra as u32;
    let q_end = q_seed + rb as u32;
    Some((total_score, q_start, q_end, s_start, s_end))
}

/// Gapped alignment with full traceback (ALIGN_EX port).
///
/// `a` = query bases (BLASTNA, 0-indexed).
/// `b` = subject bases (BLASTNA, 0-indexed).
/// `gap_open`, `gap_extend` — affine gap penalties (positive values).
/// `x_dropoff` — drop from best score that triggers pruning.
/// `reverse` — if true, scan a and b from the end (left extension).
/// `ws` — caller-supplied workspace; contents are overwritten.
pub fn align_ex(
    a: &[u8],
    b: &[u8],
    m: usize,
    n: usize,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    matrix: &ScoreMatrix,
    reverse: bool,
    ws: &mut AlignWorkspace,
) -> GapAlignResult {
    if reverse {
        align_ex_inner::<true>(a, b, m, n, gap_open, gap_extend, x_dropoff, matrix, ws)
    } else {
        align_ex_inner::<false>(a, b, m, n, gap_open, gap_extend, x_dropoff, matrix, ws)
    }
}

/// Monomorphised inner function — `REVERSE` is a compile-time constant so
/// the branch is eliminated and the compiler generates specialised code.
#[inline(never)]
fn align_ex_inner<const REVERSE: bool>(
    a: &[u8],
    b: &[u8],
    m: usize,
    n: usize,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    matrix: &ScoreMatrix,
    ws: &mut AlignWorkspace,
) -> GapAlignResult {
    let m = m.min(a.len());
    let n = n.min(b.len());

    if m == 0 || n == 0 {
        return GapAlignResult { score: 0, a_len: 0, b_len: 0, edit_script: EditScript::new() };
    }

    let gap_open_extend = gap_open + gap_extend;
    let xdrop = x_dropoff.max(gap_open_extend);
    let num_extra = if gap_extend > 0 { (xdrop / gap_extend + 3) as usize } else { n + 3 };
    let max_cap = n + 2;

    // Faithful port of NCBI ALIGN_EX (blast_gapalign.c:411+).  No score-only pre-pass:
    // the band is discovered in a SINGLE pass, and the two backing buffers grow on
    // demand, so memory is O(total band cells) rather than the old O(dp_cap^2) reserve.
    //
    // Score array (NCBI `dp_mem`): starts at num_extra+100 and is doubled when the
    // band's right edge nears the allocation (line 472/659), capped at n+2 (the
    // largest index ever touched is the sentinel at b_size <= n+1).  `dp` is grown in
    // place; `&mut ws.dp` and the `ws.flat_edit`/`ws.row_info` pushes below are
    // disjoint struct-field borrows, so no raw pointer / refresh dance is needed.
    let mut dp_cap = (num_extra + 101).min(max_cap);
    if ws.dp.len() < dp_cap {
        ws.dp.resize(dp_cap, DpCell { best: MININT, best_gap: MININT });
    }
    // Raw cursor into ws.dp so the score array can be resized (realloc) mid-pass
    // without a long-lived borrow conflicting with the ws.flat_edit/ws.row_info
    // pushes.  Refreshed after every resize (NCBI keeps `score_array = dp_mem`).
    let mut dp_ptr = ws.dp.as_mut_ptr();
    let mut dp_high = dp_cap; // highest index initialised (reset range at the end)

    // Traceback storage (NCBI per-row `state_struct`, allocated by s_GapGetState).
    // flat_edit grows incrementally — reserved per row to the row's worst case, then
    // written through a raw cursor with no per-cell capacity check (this is the hot
    // path RepeatMasker runs constantly).  `flat_len` is the authoritative logical
    // length; `flat_ptr` is refreshed after every (possibly reallocating) reserve.
    // row_info[a] = (flat offset of row a, first_b_index).
    //
    // SAFETY: every `*flat_ptr.add(flat_len)` write is preceded, since the last pointer
    // refresh, by a `reserve` of this row's worst-case cell count, so flat_len is always
    // < capacity at the write.  flat_len advances only by bytes actually written, so the
    // closing `set_len(flat_len)` exposes only initialised bytes.  Resizing ws.dp does
    // not touch ws.flat_edit's allocation, so flat_ptr stays valid across the dp grow.
    ws.flat_edit.clear();
    ws.row_info.clear();
    ws.flat_edit.reserve(num_extra + 104); // row 0: GAP_IN_B + init loop (bounded by dp_cap)
    let mut flat_ptr = ws.flat_edit.as_mut_ptr();
    let mut flat_len = 0usize;

    // ------ Initialise dp[0..b_size) and row 0 of traceback ------
    // SAFETY: index 0 < dp_cap <= ws.dp.len(); init loop indices j < dp_cap.
    unsafe {
        let c = &mut *dp_ptr.add(0);
        c.best = 0;
        c.best_gap = -gap_open_extend;
    }
    ws.row_info.push((0, 0));
    unsafe { *flat_ptr.add(flat_len) = SCRIPT_GAP_IN_B; }
    flat_len += 1;

    let mut score = -gap_open_extend;
    let mut b_size = 1usize;
    for j in 1..=n {
        if score < -xdrop || j >= dp_cap { break; }
        unsafe {
            let c = &mut *dp_ptr.add(j);
            c.best = score;
            c.best_gap = score - gap_open_extend;
        }
        unsafe { *flat_ptr.add(flat_len) = SCRIPT_GAP_IN_A; }
        flat_len += 1;
        score -= gap_extend;
        b_size = j + 1;
    }

    let mut best_score = 0i32;
    let mut best_a = 0usize;
    let mut best_b = 0usize;
    let mut first_b_index = 0usize;

    let mut total_cells: u64 = 0;
    // ------ Main DP loop ------
    for a_idx in 0..m {
        let a_phys = if REVERSE { m - 1 - a_idx } else { a_idx };
        // SAFETY: a_phys < m <= a.len(); a_base masked to 0..=15; matrix has 16 rows.
        let matrix_row = unsafe {
            let a_base = (*a.get_unchecked(a_phys) & 15) as usize;
            matrix.scores.get_unchecked(a_base)
        };
        let orig_first = first_b_index;

        // Record the start of this row, then reserve the row's worst-case cells
        // (main loop writes `b_size - first_b_index`; the right extension adds at most
        // `num_extra`) and refresh the raw cursor (NCBI s_GapGetState allocates per row).
        ws.row_info.push((flat_len, orig_first));
        unsafe { ws.flat_edit.set_len(flat_len); } // commit prior rows so reserve accounts for them
        ws.flat_edit.reserve((b_size - first_b_index) + num_extra + 4);
        flat_ptr = ws.flat_edit.as_mut_ptr();

        let mut score = MININT;
        let mut score_gap_row = MININT;
        let mut last_b_index = first_b_index;

        total_cells += (b_size - first_b_index) as u64;
        // SAFETY: b_idx iterates first_b_index..b_size where b_size <= dp_cap == dp.len().
        //   b_base = (*b.get_unchecked(b_phys) & 15) as usize is always 0..=15;
        //   matrix.scores and each row have exactly BLASTNA_SIZE (16) elements.
        //   b_phys < n <= b.len() whenever b_idx < n.
        for b_idx in first_b_index..b_size {
            // NCBI nucleotide BLAST: FENCE_SENTRY=201 never fires; sentinels are NULLB=0.
            // At b_idx==n (one past sequence end) substitute NULLB rather than OOB-access.
            let b_base = if b_idx >= n {
                0usize
            } else {
                let b_phys = if REVERSE { n - 1 - b_idx } else { b_idx };
                unsafe { (*b.get_unchecked(b_phys) & 15) as usize }
            };

            let (mut score_gap_col, next_score) = unsafe {
                let cell = &*dp_ptr.add(b_idx);
                (cell.best_gap, cell.best + *matrix_row.get_unchecked(b_base))
            };

            let mut cell_op = SCRIPT_SUB;
            let mut cell_score = score;

            if score_gap_col > cell_score { cell_score = score_gap_col; cell_op = SCRIPT_GAP_IN_B; }
            if score_gap_row > cell_score { cell_score = score_gap_row; cell_op = SCRIPT_GAP_IN_A; }

            if best_score - cell_score > xdrop {
                if b_idx == first_b_index {
                    first_b_index += 1;
                } else {
                    unsafe { (*dp_ptr.add(b_idx)).best = MININT; }
                }
            } else {
                last_b_index = b_idx;
                if cell_score > best_score {
                    best_score = cell_score;
                    best_a = a_idx + 1;
                    best_b = b_idx;
                }

                score_gap_row -= gap_extend;
                score_gap_col -= gap_extend;
                let open_gap = cell_score - gap_open_extend;

                unsafe {
                    let cell = &mut *dp_ptr.add(b_idx);
                    if score_gap_col < open_gap {
                        cell.best_gap = open_gap;
                    } else {
                        cell.best_gap = score_gap_col;
                        cell_op |= SCRIPT_EXTEND_GAP_B;
                    }
                    cell.best = cell_score;
                }

                if score_gap_row < open_gap {
                    score_gap_row = open_gap;
                } else {
                    cell_op |= SCRIPT_EXTEND_GAP_A;
                }
            }

            score = next_score;
            unsafe { *flat_ptr.add(flat_len) = cell_op; }
            flat_len += 1;
        }

        if first_b_index >= b_size { break; }

        // Grow the score array by doubling when the band's right edge nears the
        // allocation (NCBI blast_gapalign.c:659), capped at n+2.  After this the
        // right-extension and sentinel writes below can never exceed dp_cap, so the
        // `< dp_cap` guards are safety nets that never bind (no #36 truncation).
        if last_b_index + num_extra + 3 >= dp_cap && dp_cap < max_cap {
            let new_cap = (last_b_index + num_extra + 100).max(dp_cap * 2).min(max_cap);
            if ws.dp.len() < new_cap {
                ws.dp.resize(new_cap, DpCell { best: MININT, best_gap: MININT });
            }
            dp_cap = new_cap;
            if new_cap > dp_high { dp_high = new_cap; }
            dp_ptr = ws.dp.as_mut_ptr(); // realloc may have moved the buffer
        }

        // Right extension.
        if last_b_index + 1 < b_size {
            b_size = last_b_index + 1;
        } else {
            while score_gap_row >= best_score - xdrop && b_size <= n && b_size < dp_cap {
                unsafe {
                    let c = &mut *dp_ptr.add(b_size);
                    c.best = score_gap_row;
                    c.best_gap = score_gap_row - gap_open_extend;
                }
                score_gap_row -= gap_extend;
                unsafe { *flat_ptr.add(flat_len) = SCRIPT_GAP_IN_A; }
                flat_len += 1;
                b_size += 1;
            }
        }
        if b_size <= n && b_size < dp_cap {
            unsafe {
                let c = &mut *dp_ptr.add(b_size);
                c.best = MININT;
                c.best_gap = MININT;
            }
            b_size += 1;
        }

    }

    TOTAL_DP_CELLS.fetch_add(total_cells, Ordering::Relaxed);

    // Commit the bytes written through flat_ptr — all of [0, flat_len) are initialised.
    unsafe { ws.flat_edit.set_len(flat_len); }
    let edit = traceback_flat(&ws.flat_edit, &ws.row_info, best_a, best_b);

    // Reinitialise the dp cells this call touched so the workspace stays clean for
    // the next call.  (NCBI relies on the within-call write-before-read invariant and
    // does not clear dp_mem, but we keep the reset to preserve established behaviour.)
    let reset_to = dp_high.min(ws.dp.len());
    for cell in ws.dp[..reset_to].iter_mut() {
        *cell = DpCell { best: MININT, best_gap: MININT };
    }

    GapAlignResult { score: best_score, a_len: best_a, b_len: best_b, edit_script: edit }
}

/// Walk the flat traceback buffer backward from (a_row, b_col) to (0, 0).
fn traceback_flat(
    flat: &[u8],
    row_info: &[(usize, usize)], // (byte_offset, first_b_of_row)
    end_a: usize,
    end_b: usize,
) -> EditScript {
    let mut result = EditScript::new();
    if end_a == 0 && end_b == 0 { return result; }

    let mut ops_rev: Vec<(EditOp, u32)> = Vec::new();
    let mut a = end_a;
    let mut b = end_b;
    let mut prev_op = SCRIPT_SUB;

    while a > 0 || b > 0 {
        if a >= row_info.len() {
            push_rev(&mut ops_rev, EditOp::GapInQuery);
            if b > 0 { b -= 1; }
            continue;
        }

        let (byte_off, first_b) = row_info[a];
        let row_end = if a + 1 < row_info.len() { row_info[a + 1].0 } else { flat.len() };

        let raw_op = if b >= first_b {
            let byte_idx = byte_off + (b - first_b);
            if byte_idx < row_end { flat[byte_idx] } else { SCRIPT_GAP_IN_A }
        } else {
            SCRIPT_GAP_IN_A
        };

        let op = match prev_op & SCRIPT_OP_MASK {
            SCRIPT_GAP_IN_A => {
                if raw_op & SCRIPT_EXTEND_GAP_A != 0 { SCRIPT_GAP_IN_A }
                else { raw_op & SCRIPT_OP_MASK }
            }
            SCRIPT_GAP_IN_B => {
                if raw_op & SCRIPT_EXTEND_GAP_B != 0 { SCRIPT_GAP_IN_B }
                else { raw_op & SCRIPT_OP_MASK }
            }
            _ => raw_op & SCRIPT_OP_MASK,
        };

        prev_op = op;

        match op {
            SCRIPT_GAP_IN_A => { push_rev(&mut ops_rev, EditOp::GapInQuery);   if b > 0 { b -= 1; } }
            SCRIPT_GAP_IN_B => { push_rev(&mut ops_rev, EditOp::GapInSubject); if a > 0 { a -= 1; } }
            _ =>               { push_rev(&mut ops_rev, EditOp::Sub);
                                  if a > 0 { a -= 1; }
                                  if b > 0 { b -= 1; } }
        }

        if a == 0 && b == 0 { break; }
    }

    ops_rev.reverse();
    for (op, n) in ops_rev { result.push(op, n); }
    result
}

fn push_rev(ops: &mut Vec<(EditOp, u32)>, op: EditOp) {
    if let Some(last) = ops.last_mut() {
        if last.0 == op { last.1 += 1; return; }
    }
    ops.push((op, 1));
}

/// Convenience: extract aligned sequence pair from a and b using the edit script.
pub fn extract_aligned(
    a: &[u8],
    b: &[u8],
    a_start: usize,
    b_start: usize,
    script: &EditScript,
    n_mask: &[u8],
) -> (Vec<u8>, Vec<u8>) {
    // n_mask[si] != 0 → restore the original BLASTNA ambiguity code (14=N,
    // 4-13=R/Y/M/K/W/S/B/D/H/V) overriding the random base used for scanning.
    let subj_base = |si: usize| -> u8 {
        let b = if si < b.len() { b[si] & 15 } else { 14 };
        if si < n_mask.len() && n_mask[si] != 0 { n_mask[si] } else { b }
    };

    let mut qa = Vec::new();
    let mut sa = Vec::new();
    let mut qi = a_start;
    let mut si = b_start;

    for &(op, count) in &script.ops {
        let n = count as usize;
        match op {
            EditOp::Sub => {
                for _ in 0..n {
                    qa.push(if qi < a.len() { a[qi] & 15 } else { 14 });
                    sa.push(subj_base(si));
                    qi += 1; si += 1;
                }
            }
            EditOp::GapInSubject => {
                for _ in 0..n {
                    qa.push(if qi < a.len() { a[qi] & 15 } else { 14 });
                    sa.push(15);
                    qi += 1;
                }
            }
            EditOp::GapInQuery => {
                for _ in 0..n {
                    qa.push(15);
                    sa.push(subj_base(si));
                    si += 1;
                }
            }
        }
    }
    (qa, sa)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::matrix::ScoreMatrix;
    use std::io::Cursor;

    const SIMPLE_MAT: &str = r"
# FREQS A 0.25 C 0.25 G 0.25 T 0.25
   A   C   G   T
A  5  -4  -4  -4
C -4   5  -4  -4
G -4  -4   5  -4
T -4  -4  -4   5
";

    fn mat() -> ScoreMatrix { ScoreMatrix::from_reader("test", Cursor::new(SIMPLE_MAT)).unwrap() }
    fn ws() -> AlignWorkspace { AlignWorkspace::new() }

    #[test]
    fn test_perfect_alignment() {
        let m = mat(); let mut w = ws();
        let a = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let s = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let r = align_ex(&a, &s, 8, 8, 4, 4, 50, &m, false, &mut w);
        assert_eq!(r.score, 5 * 8, "score={}", r.score);
        assert_eq!(r.a_len, 8);
        assert_eq!(r.b_len, 8);
        assert_eq!(r.edit_script.ops.len(), 1);
        assert_eq!(r.edit_script.ops[0], (EditOp::Sub, 8));
    }

    #[test]
    fn test_single_gap_in_subject() {
        let m = mat(); let mut w = ws();
        let a = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let s = vec![0u8, 1, 2, 3, 1, 2, 3];
        let r = align_ex(&a, &s, 8, 7, 4, 4, 50, &m, false, &mut w);
        assert!(r.score > 0, "score={}", r.score);
        assert!(r.a_len > 0);
        assert!(r.b_len > 0);
    }

    #[test]
    fn test_empty_sequence() {
        let m = mat(); let mut w = ws();
        let a: Vec<u8> = vec![];
        let s: Vec<u8> = vec![0, 1, 2, 3];
        let r = align_ex(&a, &s, 0, 4, 4, 4, 50, &m, false, &mut w);
        assert_eq!(r.score, 0);
    }

    #[test]
    fn test_extract_aligned_no_gap() {
        let a = vec![0u8, 1, 2, 3];
        let s = vec![0u8, 1, 2, 3];
        let mut script = EditScript::new();
        script.push(EditOp::Sub, 4);
        let (qa, sa) = extract_aligned(&a, &s, 0, 0, &script, &[]);
        assert_eq!(qa, vec![0, 1, 2, 3]);
        assert_eq!(sa, vec![0, 1, 2, 3]);
    }
}
