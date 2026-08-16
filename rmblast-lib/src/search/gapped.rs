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

/// Diagnostic counter: number of times the REVERSE band pointer had to be clamped
/// because `first_b_index` reached `n` (see `align_ex_score_only_inner`).  The clamp
/// is behaviour-neutral — the inner loop runs zero times in that state — this counter
/// only records how often the degenerate band state occurs.
pub static REVERSE_FBI_CLAMPED: AtomicU64 = AtomicU64::new(0);

/// Diagnostic counter: rows processed by the AVX2 pass kernel.
pub static SCORE_ONLY_SIMD_ROWS: AtomicU64 = AtomicU64::new(0);

/// Diagnostic counter: total score-only DP rows.
pub static SCORE_ONLY_ROWS: AtomicU64 = AtomicU64::new(0);

/// Diagnostic counter: total score-only DP cells in rows of width >= 16.
pub static SCORE_ONLY_WIDE_CELLS: AtomicU64 = AtomicU64::new(0);

/// Diagnostic counter: total score-only DP cells (subset of TOTAL_DP_CELLS).
pub static SCORE_ONLY_CELLS: AtomicU64 = AtomicU64::new(0);

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
#[repr(C)] // guaranteed [best, best_gap] layout — the SIMD kernel de/re-interleaves pairs
pub struct DpCell {
    best: i32,
    best_gap: i32,
}

/// Whether the AVX2 score-only row kernel is available and enabled.
/// `RMBLAST_NO_SIMD=1` forces the scalar reference kernel (for A/B validation).
#[cfg(target_arch = "x86_64")]
fn simd_enabled() -> bool {
    static ENABLED: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ENABLED.get_or_init(|| {
        std::env::var_os("RMBLAST_NO_SIMD").is_none()
            && std::arch::is_x86_feature_detected!("avx2")
    })
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
    // The SIMD row kernel requires AVX2 and positive gap costs (its
    // pruned-cell equivalence argument needs gap_extend > 0 and
    // gap_open_extend > 0); otherwise take the reference loop.
    #[cfg(target_arch = "x86_64")]
    if gap_extend > 0 && gap_open + gap_extend > 0 && simd_enabled() {
        return if reverse {
            align_ex_score_only_inner::<true>(a, b, m, n, gap_open, gap_extend, x_dropoff, matrix, dp, dump)
        } else {
            align_ex_score_only_inner::<false>(a, b, m, n, gap_open, gap_extend, x_dropoff, matrix, dp, dump)
        };
    }
    if reverse {
        align_ex_score_only_inner_scalar::<true>(a, b, m, n, gap_open, gap_extend, x_dropoff, matrix, dp, dump)
    } else {
        align_ex_score_only_inner_scalar::<false>(a, b, m, n, gap_open, gap_extend, x_dropoff, matrix, dp, dump)
    }
}

#[inline(never)]
fn align_ex_score_only_inner_scalar<const REVERSE: bool>(
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
            // `first_b_index` can reach `n` (whole band x-dropped except the b_idx==n
            // sentinel, which keeps b_size at n+1 so the `first_b_index >= b_size`
            // break below does not fire).  The inner loop then runs zero times
            // (b_idx == first_b_index == inner_end), so the pointer is never read —
            // but forming it would underflow `n - 1 - first_b_index` (debug panic,
            // out-of-bounds `add` = UB in release).  Substitute a valid placeholder.
            if first_b_index < n {
                unsafe { b.as_ptr().add(n - 1 - first_b_index) }
            } else {
                REVERSE_FBI_CLAMPED.fetch_add(1, Ordering::Relaxed);
                b.as_ptr()
            }
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
    SCORE_ONLY_CELLS.fetch_add(total_cells_so, Ordering::Relaxed);
    let reset_end = b_size.min(dp.len());
    for cell in dp[..reset_end].iter_mut() {
        *cell = DpCell { best: MININT, best_gap: MININT };
    }
    (best_score, best_a, best_b)
}

/// AVX2 implementation of the fused score-only row kernel.
///
/// One pass per row: de-interleaved dp loads, gathered matrix scores, the
/// always-update row-gap chain as a weighted Kogge-Stone max-scan, the
/// exclusive running-best prefix scan, x-drop classification, and blended
/// dp stores — with the diagonal chain carried in a register between blocks.
/// See `align_ex_score_only_inner` for the equivalence argument; the
/// differential fuzz test checks it against the reference loop.
#[cfg(target_arch = "x86_64")]
mod so_simd {
    use super::{DpCell, MININT};
    use std::arch::x86_64::*;

    /// Fill value for shifted-in scan lanes.  Strictly below every legitimate
    /// score (live cells sit above ~MININT - goe - matrix_max) yet with enough
    /// headroom that subtracting the scan decays (≤ 7·gap_extend) cannot wrap.
    /// The fast-path dispatch bounds gap params, so 2M of slack is ample.
    const NEG_FILL: i32 = MININT - 2_000_000;

    pub struct RowOut {
        pub leading: usize, // row-relative index of first surviving cell (w if none)
        pub seen_np: bool,
        pub last_np: usize, // row-relative index of last surviving cell
        pub p: i32,         // running best after the row
        pub best_b_rel: usize,
        pub improved: bool,
        pub f: i32,         // final always-update row-gap chain value
        pub t_last: i32,    // diagonal chain value after the last cell
    }

    /// lane i ← lane i-1 (i32), lane 0 ← fill (any lane of `fv`).
    #[inline(always)]
    unsafe fn sh1(v: __m256i, fv: __m256i) -> __m256i {
        let t = _mm256_permute2x128_si256(v, fv, 0x03); // [fv.high | v.low]
        _mm256_alignr_epi8(v, t, 12)
    }
    /// lane i ← lane i-2, lanes 0-1 ← fill.
    #[inline(always)]
    unsafe fn sh2(v: __m256i, fv: __m256i) -> __m256i {
        let t = _mm256_permute2x128_si256(v, fv, 0x03);
        _mm256_alignr_epi8(v, t, 8)
    }
    /// lane i ← lane i-4, lanes 0-3 ← fill.
    #[inline(always)]
    unsafe fn sh4(v: __m256i, fv: __m256i) -> __m256i {
        _mm256_permute2x128_si256(v, fv, 0x03)
    }

    /// De-interleave 8 DpCells starting at `p` into (best, best_gap) vectors.
    #[inline(always)]
    unsafe fn load_cells(p: *const DpCell) -> (__m256i, __m256i) {
        let v0 = _mm256_loadu_si256(p as *const __m256i); // b0 g0 b1 g1 | b2 g2 b3 g3
        let v1 = _mm256_loadu_si256(p.add(4) as *const __m256i); // b4 g4 b5 g5 | b6 g6 b7 g7
        let best_ps = _mm256_shuffle_ps(_mm256_castsi256_ps(v0), _mm256_castsi256_ps(v1), 0x88);
        let gap_ps = _mm256_shuffle_ps(_mm256_castsi256_ps(v0), _mm256_castsi256_ps(v1), 0xDD);
        // shuffle_ps yields [b0 b1 b4 b5 | b2 b3 b6 b7]; fix 64-bit block order.
        let best = _mm256_permute4x64_epi64(_mm256_castps_si256(best_ps), 0b11_01_10_00);
        let gap = _mm256_permute4x64_epi64(_mm256_castps_si256(gap_ps), 0b11_01_10_00);
        (best, gap)
    }

    /// Interleave (best, best_gap) vectors back into 8 DpCells at `p`.
    #[inline(always)]
    unsafe fn store_cells(p: *mut DpCell, best: __m256i, gap: __m256i) {
        let lo = _mm256_unpacklo_epi32(best, gap); // b0 g0 b1 g1 | b4 g4 b5 g5
        let hi = _mm256_unpackhi_epi32(best, gap); // b2 g2 b3 g3 | b6 g6 b7 g7
        let out0 = _mm256_permute2x128_si256(lo, hi, 0x20);
        let out1 = _mm256_permute2x128_si256(lo, hi, 0x31);
        _mm256_storeu_si256(p as *mut __m256i, out0);
        _mm256_storeu_si256(p.add(4) as *mut __m256i, out1);
    }

    /// Fused single-pass row kernel: computes the row and commits dp updates
    /// block-by-block.  Requires w >= 8 (callers use a higher threshold).
    ///
    /// # Safety
    /// Caller guarantees: AVX2 available; `fbi + w <= dp.len()`; for FORWARD
    /// `fbi + w <= n <= b.len()`, for REVERSE the row maps to
    /// `b[n-1-(fbi+w-1) ..= n-1-fbi]`, in range for the same reason.
    #[target_feature(enable = "avx2")]
    pub unsafe fn row_fused<const REVERSE: bool>(
        b: &[u8],
        n: usize,
        fbi: usize,
        w: usize,
        dp: &mut [DpCell],
        matrix_row: &[i32; 16],
        e: i32,
        goe: i32,
        xdrop: i32,
        best_in: i32,
    ) -> RowOut {
        let full = w / 8 * 8;
        let low15 = _mm256_set1_epi32(15);
        let rev_idx = _mm256_setr_epi32(7, 6, 5, 4, 3, 2, 1, 0);
        let idxv = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
        let minint_v = _mm256_set1_epi32(MININT);
        let negf_v = _mm256_set1_epi32(NEG_FILL);
        let e_v = _mm256_set1_epi32(e);
        let e2_v = _mm256_set1_epi32(2 * e);
        let e4_v = _mm256_set1_epi32(4 * e);
        let goe_v = _mm256_set1_epi32(goe);
        let xdrop_v = _mm256_set1_epi32(xdrop);
        let ramp = _mm256_setr_epi32(0, e, 2 * e, 3 * e, 4 * e, 5 * e, 6 * e, 7 * e);

        let mut st = RowOut {
            leading: w,
            seen_np: false,
            last_np: 0,
            p: best_in,
            best_b_rel: 0,
            improved: false,
            f: MININT,
            t_last: MININT,
        };
        // prev_t lane 7 supplies the diagonal carry t[i-1] into each block's
        // lane 0 (sh1 shifts in fv's top lane); MININT for the row's first cell.
        let mut prev_t = minint_v;

        let mut blk = 0usize;
        while blk < full {
            let jp = dp.as_ptr().add(fbi + blk);
            let (bestv, gapv) = load_cells(jp);
            let idx = if REVERSE {
                let raw = _mm_loadl_epi64(b.as_ptr().add(n - 1 - (fbi + blk + 7)) as *const __m128i);
                let ext = _mm256_cvtepu8_epi32(raw);
                _mm256_permutevar8x32_epi32(ext, rev_idx)
            } else {
                let raw = _mm_loadl_epi64(b.as_ptr().add(fbi + blk) as *const __m128i);
                _mm256_cvtepu8_epi32(raw)
            };
            let idx = _mm256_and_si256(idx, low15);
            let ms = _mm256_i32gather_epi32(matrix_row.as_ptr(), idx, 4);
            let tv = _mm256_add_epi32(bestv, ms);
            let dv = sh1(tv, prev_t);
            let ng = _mm256_max_epi32(dv, gapv);
            let g = _mm256_sub_epi32(ng, goe_v);

            // H[i] = max over k<=i of (g[k] - (i-k)*e): weighted Kogge-Stone.
            let mut h = g;
            h = _mm256_max_epi32(h, _mm256_sub_epi32(sh1(h, negf_v), e_v));
            h = _mm256_max_epi32(h, _mm256_sub_epi32(sh2(h, negf_v), e2_v));
            h = _mm256_max_epi32(h, _mm256_sub_epi32(sh4(h, negf_v), e4_v));

            // F[i] = max(f_in - i*e, H[i-1]); cell = max(ng, F).
            let fdec = _mm256_sub_epi32(_mm256_set1_epi32(st.f), ramp);
            let fv = _mm256_max_epi32(fdec, sh1(h, negf_v));
            let cellv = _mm256_max_epi32(ng, fv);

            // Exclusive running-best prefix and prune classification.
            let mut pmv = cellv;
            pmv = _mm256_max_epi32(pmv, sh1(pmv, negf_v));
            pmv = _mm256_max_epi32(pmv, sh2(pmv, negf_v));
            pmv = _mm256_max_epi32(pmv, sh4(pmv, negf_v));
            let pexcl = _mm256_max_epi32(_mm256_set1_epi32(st.p), sh1(pmv, negf_v));
            let prunev = _mm256_cmpgt_epi32(_mm256_sub_epi32(pexcl, cellv), xdrop_v);
            let pbits = _mm256_movemask_ps(_mm256_castsi256_ps(prunev)) as u32;
            let sv = !pbits & 0xFF;

            // Row-gap and diagonal carries (needed on every path).
            let h7 = _mm256_extract_epi32(h, 7);
            let f_next = (st.f - 8 * e).max(h7);

            let mut keep_lanes = 0u32;
            if !st.seen_np {
                if sv == 0 {
                    // Whole block still in the leading pruned run: cells stay
                    // untouched, nothing else to track.
                    st.f = f_next;
                    prev_t = tv;
                    blk += 8;
                    continue;
                }
                let lead = sv.trailing_zeros();
                st.leading = blk + lead as usize;
                st.seen_np = true;
                keep_lanes = lead;
            }
            if sv != 0 {
                st.last_np = blk + (31 - sv.leading_zeros()) as usize;
            }

            // Best / argbest (first attainment of the improved value).
            let pm7 = _mm256_extract_epi32(pmv, 7);
            if pm7 > st.p {
                let eq = _mm256_movemask_ps(_mm256_castsi256_ps(
                    _mm256_cmpeq_epi32(cellv, _mm256_set1_epi32(pm7)),
                )) as u32;
                st.best_b_rel = blk + eq.trailing_zeros() as usize;
                st.p = pm7;
                st.improved = true;
            }

            // Commit: pruned cells keep stale best_gap and take MININT best;
            // surviving cells take the standard update; leading lanes (only in
            // the block where the leading run ends) keep their original pair.
            let newgap = _mm256_blendv_epi8(
                _mm256_max_epi32(_mm256_sub_epi32(gapv, e_v), _mm256_sub_epi32(cellv, goe_v)),
                gapv,
                prunev,
            );
            let newbest = _mm256_blendv_epi8(cellv, minint_v, prunev);
            let (newbest, newgap) = if keep_lanes > 0 {
                let keep = _mm256_cmpgt_epi32(_mm256_set1_epi32(keep_lanes as i32), idxv);
                (
                    _mm256_blendv_epi8(newbest, bestv, keep),
                    _mm256_blendv_epi8(newgap, gapv, keep),
                )
            } else {
                (newbest, newgap)
            };
            store_cells(jp as *mut DpCell, newbest, newgap);

            st.f = f_next;
            prev_t = tv;
            blk += 8;
        }

        // Scalar ragged tail (same fused semantics per cell).
        let mut t_prev = _mm256_extract_epi32(prev_t, 7);
        for i in full..w {
            let j = fbi + i;
            let b_phys = if REVERSE { n - 1 - j } else { j };
            let b_base = (*b.get_unchecked(b_phys) & 15) as usize;
            let cellp = dp.get_unchecked(j);
            let told = cellp.best + *matrix_row.get_unchecked(b_base);
            let cval = cellp.best_gap;
            let ng = if cval > t_prev { cval } else { t_prev };
            let cl = if st.f > ng { st.f } else { ng };
            let pruned = st.p - cl > xdrop;
            if pruned {
                if st.seen_np {
                    dp.get_unchecked_mut(j).best = MININT;
                } // else: leading run, cell stays untouched
            } else {
                if !st.seen_np {
                    st.seen_np = true;
                    st.leading = i;
                }
                st.last_np = i;
                if cl > st.p {
                    st.p = cl;
                    st.best_b_rel = i;
                    st.improved = true;
                }
                let open_gap = cl - goe;
                let cext = cval - e;
                let cellp = dp.get_unchecked_mut(j);
                cellp.best_gap = if cext > open_gap { cext } else { open_gap };
                cellp.best = cl;
            }
            let open = ng - goe;
            let fe = st.f - e;
            st.f = if fe > open { fe } else { open };
            t_prev = told;
        }
        st.t_last = t_prev;
        st
    }
}

/// AVX2 variant of [`align_ex_score_only_inner_scalar`] — identical results.
/// Rows with a band of >= 16 cells run the fused SIMD kernel
/// ([`so_simd::row_fused`]); narrower rows run the reference per-cell loop.
///
/// # Equivalence of the SIMD kernel to the reference loop
///
/// The reference loop updates the row-gap chain (`score_gap_row`) only at
/// cells that survive the x-drop test.  The kernel instead computes an
/// **always-update** chain `F'[i+1] = max(F'[i] - e, nogap[i] - goe)` as a
/// weighted prefix-max scan.  With `e > 0`, `goe > 0` (enforced by the
/// dispatch) the two chains agree on every observable outcome:
///
/// * Wherever the chains differ, both values are `< running_best - xdrop`.
///   For the reference chain: a value frozen at a pruned cell k is bounded by
///   that cell's score (`cell[k] = max(nogap[k], F) >= F`), which the prune
///   test placed below `best_at[k] - xdrop`.  For the always-update chain:
///   contributions taken from pruned cells are `nogap[k] - goe`, bounded the
///   same way.  Both bounds persist — the running best only grows, decay only
///   lowers the values, and refresh terms (`cell - goe` at surviving cells)
///   are common to both chains.
/// * A chain value below `running_best - xdrop` can never determine anything:
///   at a surviving cell the winner must be `>= best - xdrop` (else the cell
///   would have been pruned — and if the chain value exceeds `nogap`, the cell
///   IS pruned in both versions); the band-tail extension requires
///   `score_gap_row >= best - xdrop`; and the strict-improvement best/argbest
///   updates only look at surviving cells, whose values agree.
///
/// Col-gap skip-on-prune is replicated directly (pruned cells keep their stale
/// `best_gap`), as are leading-prune `first_b_index` advances (leading cells
/// keep their stale `best`), the strict-improvement `best_b` tiebreak, and the
/// sentinel / band-tail / grow scaffolding (copied verbatim).  The
/// differential fuzz test validates all of this against the reference loop.
#[cfg(target_arch = "x86_64")]
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
    // gap_extend > 0 is guaranteed by the dispatch in align_ex_score_only.
    let num_extra = (xdrop / gap_extend + 3) as usize;
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
        let inner_end = b_size.min(n);
        let w = inner_end.saturating_sub(first_b_index);
        SCORE_ONLY_ROWS.fetch_add(1, Ordering::Relaxed);

        if w >= 16 {
            SCORE_ONLY_SIMD_ROWS.fetch_add(1, Ordering::Relaxed);
            SCORE_ONLY_WIDE_CELLS.fetch_add(w as u64, Ordering::Relaxed);
            // SAFETY: AVX2 checked by the dispatch (simd_enabled());
            // first_b_index + w = inner_end <= min(b_size, n) <= dp.len() and
            // <= b.len() in FORWARD (REVERSE maps to b[n-1-(fbi+w-1)..=n-1-fbi],
            // in range for the same reason).
            let st = unsafe {
                so_simd::row_fused::<REVERSE>(
                    b, n, first_b_index, w, dp, matrix_row,
                    gap_extend, gap_open_extend, xdrop, best_score,
                )
            };
            let fbi0 = first_b_index;
            first_b_index = fbi0 + st.leading;
            if st.seen_np { last_b_index = fbi0 + st.last_np; }
            if st.improved {
                best_score = st.p;
                best_b = fbi0 + st.best_b_rel;
            }
            score = st.t_last;
            score_gap_row = st.f;
        } else if w > 0 {
            // Narrow band: reference per-cell loop (exact skip-on-prune chain).
            let mut b_idx = first_b_index;
            while b_idx < inner_end {
                let b_phys = if REVERSE { n - 1 - b_idx } else { b_idx };
                let b_base = unsafe { (*b.get_unchecked(b_phys) & 15) as usize };
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
        } else if REVERSE && first_b_index >= n {
            // Degenerate band-state diagnostic parity with the reference loop.
            REVERSE_FBI_CLAMPED.fetch_add(1, Ordering::Relaxed);
        }

        // Sentinel position (b_idx == n, b_base = 0 = NULLB) — verbatim from
        // the reference loop.
        if b_size > n {
            let b_idx = n;
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
    SCORE_ONLY_CELLS.fetch_add(total_cells_so, Ordering::Relaxed);
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

    /// Regression: `align_ex_score_only` with `reverse = true` (the LEFT extension of
    /// `gapped_extend_score_only`) used to compute its band pointer as
    /// `b.as_ptr().add(n - 1 - first_b_index)`.  `first_b_index` can reach `n` when the
    /// whole band x-drops except the b_idx == n sentinel cell (which keeps `b_size` at
    /// n+1, so the `first_b_index >= b_size` loop break does not fire).  That underflowed
    /// the usize subtraction: a panic in debug builds, an out-of-bounds `add` in release.
    /// Random low-identity pairs under a tight x-dropoff drive the band into that state.
    #[test]
    fn test_reverse_band_no_underflow() {
        let m = mat();
        let mut dp: Vec<DpCell> = Vec::new();
        // xorshift64* — deterministic, no rand dependency.
        let mut state: u64 = 0x9E3779B97F4A7C15;
        let mut next = move || {
            state ^= state << 13; state ^= state >> 7; state ^= state << 17;
            state
        };
        for trial in 0..4000 {
            let alen = 1 + (next() % 60) as usize;
            let blen = 1 + (next() % 60) as usize;
            let a: Vec<u8> = (0..alen).map(|_| (next() % 4) as u8).collect();
            let b: Vec<u8> = (0..blen).map(|_| (next() % 4) as u8).collect();
            // Tight x-dropoff is what forces aggressive band pruning.
            let xdrop = 1 + (next() % 20) as i32;
            for &reverse in &[true, false] {
                let (score, ..) = align_ex_score_only(
                    &a, &b, alen, blen, 8, 2, xdrop, &m, reverse, &mut dp, false,
                );
                assert!(score >= 0, "trial {} reverse {} score {}", trial, reverse, score);
            }
        }
        // Guard the guard: confirm these inputs really do drive the band into the
        // degenerate state, so the test cannot silently stop covering the bug.
        assert!(REVERSE_FBI_CLAMPED.load(Ordering::Relaxed) > 0,
                "no reverse-band clamp observed — test inputs no longer reach the bug");
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

    // ────────────────────────────────────────────────────────────────────────
    // Differential fuzz: pass-structured score-only kernel vs reference loop
    // ────────────────────────────────────────────────────────────────────────

    struct FuzzRng(u64);
    impl FuzzRng {
        fn next(&mut self) -> u64 {
            // xorshift64* — deterministic, no external deps.
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x.wrapping_mul(0x2545F4914F6CDD1D)
        }
        fn below(&mut self, m: u64) -> u64 { self.next() % m }
        fn range_i32(&mut self, lo: i32, hi: i32) -> i32 {
            lo + (self.below((hi - lo + 1) as u64) as i32)
        }
    }

    fn fuzz_matrix(rng: &mut FuzzRng) -> ScoreMatrix {
        let mut scores = [[0i32; 16]; 16];
        match rng.below(3) {
            0 => {
                // Match/mismatch style: positive diagonal, negative off-diagonal
                // (independent entries — asymmetric).
                for i in 0..16 {
                    for j in 0..16 {
                        scores[i][j] = if i == j { rng.range_i32(1, 10) } else { rng.range_i32(-15, -1) };
                    }
                }
            }
            1 => {
                for i in 0..16 { for j in 0..16 { scores[i][j] = rng.range_i32(-30, 15); } }
            }
            _ => {
                // Extremes: the documented MININT headroom bound (±100).
                for i in 0..16 { for j in 0..16 { scores[i][j] = rng.range_i32(-100, 100); } }
            }
        }
        ScoreMatrix {
            scores,
            freqs: [0.0; 16],
            lambda: 0.0,
            name: "fuzz".into(),
            karlin: None,
        }
    }

    fn fuzz_seq(rng: &mut FuzzRng, len: usize) -> Vec<u8> {
        let mut s = Vec::with_capacity(len);
        match rng.below(4) {
            0 => { for _ in 0..len { s.push(rng.below(15) as u8); } }
            1 => {
                // Homopolymer blocks — low complexity provokes band holes.
                while s.len() < len {
                    let base = rng.below(5) as u8; // ACGT + N
                    let run = 1 + rng.below(30) as usize;
                    for _ in 0..run.min(len - s.len()) { s.push(base); }
                }
            }
            2 => {
                // Dinucleotide repeat with occasional point mutations.
                let (x, y) = (rng.below(4) as u8, rng.below(4) as u8);
                for i in 0..len {
                    let mut b = if i % 2 == 0 { x } else { y };
                    if rng.below(12) == 0 { b = rng.below(15) as u8; }
                    s.push(b);
                }
            }
            _ => {
                // Mostly ACGT with ~10% ambiguity codes.
                for _ in 0..len {
                    if rng.below(10) == 0 { s.push(rng.range_i32(4, 14) as u8); }
                    else { s.push(rng.below(4) as u8); }
                }
            }
        }
        s
    }

    #[test]
    #[cfg(target_arch = "x86_64")]
    fn fuzz_pass_structured_score_only_vs_reference() {
        if !simd_enabled() {
            eprintln!("fuzz: AVX2 unavailable or disabled — nothing to differential-test");
            return;
        }
        let iters: u64 = std::env::var("RMBLAST_FUZZ_ITERS").ok()
            .and_then(|v| v.parse().ok()).unwrap_or(20000);
        let seed: u64 = std::env::var("RMBLAST_FUZZ_SEED").ok()
            .and_then(|v| v.parse().ok()).unwrap_or(0x243F6A8885A308D3);
        let mut rng = FuzzRng(seed);
        for iter in 0u64..iters {
            let matrix = fuzz_matrix(&mut rng);
            let big = iter % 500 == 499;
            let max_len = if big { 400 } else { 120 };
            let la = 1 + rng.below(max_len) as usize;
            let lb = 1 + rng.below(max_len) as usize;
            let a = fuzz_seq(&mut rng, la);
            let b = fuzz_seq(&mut rng, lb);
            let go = rng.range_i32(0, 30);
            let ge = rng.range_i32(1, 8); // fast-path domain: gap_extend > 0
            let xd = rng.range_i32(0, 250);
            let reverse = iter % 2 == 1;

            let mut dp_new: Vec<DpCell> = Vec::new();
            let mut dp_ref: Vec<DpCell> = Vec::new();
            let got = if reverse {
                align_ex_score_only_inner::<true>(&a, &b, la, lb, go, ge, xd, &matrix, &mut dp_new, false)
            } else {
                align_ex_score_only_inner::<false>(&a, &b, la, lb, go, ge, xd, &matrix, &mut dp_new, false)
            };
            let want = if reverse {
                align_ex_score_only_inner_scalar::<true>(&a, &b, la, lb, go, ge, xd, &matrix, &mut dp_ref, false)
            } else {
                align_ex_score_only_inner_scalar::<false>(&a, &b, la, lb, go, ge, xd, &matrix, &mut dp_ref, false)
            };
            assert_eq!(
                got, want,
                "fuzz mismatch at iter={} (la={} lb={} go={} ge={} xd={} rev={})\n a={:?}\n b={:?}",
                iter, la, lb, go, ge, xd, reverse, a, b
            );
            assert_eq!(dp_new.len(), dp_ref.len(), "dp length mismatch at iter={}", iter);
        }
        let simd_rows = SCORE_ONLY_SIMD_ROWS.load(Ordering::Relaxed);
        eprintln!("fuzz: {} cases identical; simd rows exercised: {}", iters, simd_rows);
        assert!(simd_rows > 0, "fuzz never exercised the SIMD kernel — inputs too narrow to trust the run");
    }

    #[test]
    fn score_only_degenerate_gap_params_route_to_reference() {
        // gap_extend == 0 must take the reference loop (num_extra n+3 branch).
        let m = mat();
        let a = vec![0u8, 1, 2, 3, 0, 1, 2, 3];
        let b = vec![0u8, 1, 2, 0, 0, 1, 2, 3];
        let mut dp: Vec<DpCell> = Vec::new();
        let r = align_ex_score_only(&a, &b, 8, 8, 5, 0, 50, &m, false, &mut dp, false);
        let mut dp2: Vec<DpCell> = Vec::new();
        let want = align_ex_score_only_inner_scalar::<false>(&a, &b, 8, 8, 5, 0, 50, &m, &mut dp2, false);
        assert_eq!(r, want);
    }
}
