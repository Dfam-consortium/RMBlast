// blast_gapalign.rs — faithful port of the nucleotide gapped alignment core.
//
// Primary source: blast_gapalign.c (NCBI BLAST 2.17.0+)
//
// Functions ported (in dependency order):
//   s_BlastAlignPackedNucl           (lines 3033–3244)
//   s_BlastDynProgNtGappedAlignment  (lines 2948–3017)
//   ALIGN_EX                         (lines 374–767)
//   BLAST_GappedAlignmentWithTraceback (lines 4595–4700)
//   BlastGetStartForGappedAlignmentNucl (lines 3323–3390)

use super::types::{
    BlastGapAlignStruct, BlastGapDP, BlastHsp, GapAlignOp,
    GapPrelimEditBlock, ScoringParams, MININT,
};
use super::util::ncbi2na_unpack_base;
use crate::matrix::ScoreMatrix;

// ── Script constants (blast_gapalign.c lines 364–370) ────────────────────────

const SCRIPT_SUB:          u8 = 3;    // eGapAlignSub
const SCRIPT_GAP_IN_A:     u8 = 0;    // eGapAlignDel
const SCRIPT_GAP_IN_B:     u8 = 6;    // eGapAlignIns
const SCRIPT_OP_MASK:      u8 = 0x07; // mask for the base operation bits
const SCRIPT_EXTEND_GAP_A: u8 = 0x10; // continue a gap in A (row gap)
const SCRIPT_EXTEND_GAP_B: u8 = 0x40; // continue a gap in B (col gap)

// ── s_BlastAlignPackedNucl (blast_gapalign.c lines 3033–3244) ────────────────

/// No-traceback gapped extension for nucleotides.
///
/// Mirrors `static Int4 s_BlastAlignPackedNucl(...)`.
///
/// Parameters match the C signature exactly:
///   B               – query sequence, BLASTNA 1-byte/base, pointer to the
///                     extension start.  For a LEFT extension this is the
///                     beginning of the query; for a RIGHT extension it is
///                     one position before the pivot (query[pivot-1]).
///   A               – subject sequence, NCBI2NA 4-bases/byte.  For a LEFT
///                     extension this is the beginning of the packed subject
///                     (subject[0]); for a RIGHT extension it is
///                     subject[(s_length+3)/4 - 1].
///   N               – max number of query letters available (length of B).
///   M               – max number of subject letters available (in bases).
///   b_offset        – output: how far the query reached.
///   a_offset        – output: how far the subject reached.
///   gap_align       – mutable workspace; dp_mem may be grown.
///   score_params     – gap penalties.
///   reverse_sequence – TRUE → left extension (reads A/B backward).
///   x_dropoff       – drop-off threshold for this call.
///
/// Returns the best alignment score found.
pub fn s_blast_align_packed_nucl(
    b: &[u8],             // query (BLASTNA)
    a: &[u8],             // subject (NCBI2NA packed)
    n: i32,               // N: max query length
    m: i32,               // M: max subject length (in bases)
    b_offset: &mut i32,
    a_offset: &mut i32,
    gap_align: &mut BlastGapAlignStruct,
    score_params: &ScoringParams,
    matrix: &ScoreMatrix,
    reverse_sequence: bool,
    x_dropoff: i32,
) -> i32 {
    let gap_open         = score_params.gap_open;
    let gap_extend       = score_params.gap_extend;
    let gap_open_extend  = gap_open + gap_extend;
    let mat              = &matrix.scores;

    *a_offset = 0;
    *b_offset = 0;

    let mut x_dropoff = x_dropoff;
    if x_dropoff < gap_open_extend {
        x_dropoff = gap_open_extend;
    }

    if n <= 0 || m <= 0 {
        return 0;
    }

    // Allocate / grow the DP scratch array.
    // C: `num_extra_cells = x_dropoff / gap_extend + 3`  (if gap_extend > 0).
    let num_extra_cells = if gap_extend > 0 {
        x_dropoff / gap_extend + 3
    } else {
        n + 3
    };

    if num_extra_cells > gap_align.dp_mem.len() as i32 {
        let new_alloc = (num_extra_cells + 100).max(2 * gap_align.dp_mem.len() as i32);
        gap_align.dp_mem.resize(new_alloc as usize, BlastGapDP::default());
    }

    // Initialise the first row.
    let mut score = -gap_open_extend;
    gap_align.dp_mem[0].best     = 0;
    gap_align.dp_mem[0].best_gap = -gap_open_extend;

    // C: `for (i = 1; i <= N; i++) { if (score < -x_dropoff) break; ... }`
    let mut b_size: i32 = 1;
    while b_size <= n {
        if score < -x_dropoff {
            break;
        }
        gap_align.dp_mem[b_size as usize].best     = score;
        gap_align.dp_mem[b_size as usize].best_gap = score - gap_open_extend;
        score -= gap_extend;
        b_size += 1;
    }

    let mut best_score = 0i32;
    let mut first_b_index: i32 = 0;
    let b_increment: i32 = if reverse_sequence { -1 } else { 1 };

    // Outer loop: iterate over subject positions a_index = 1 .. M.
    for a_index in 1..=m {
        // Pick the score-matrix row for subject base A[a_index].
        // C: `if (reverse_sequence)  a_base_pair = NCBI2NA_UNPACK_BASE(A[(M-a_index)/4], (a_index-1)%4);`
        //    `else                   a_base_pair = NCBI2NA_UNPACK_BASE(A[1+((a_index-1)/4)], 3-((a_index-1)%4));`
        let a_base_pair: usize = if reverse_sequence {
            let byte_idx = ((m - a_index) / 4) as usize;
            let bit_pos  = ((a_index - 1) % 4) as u32;
            ncbi2na_unpack_base(a[byte_idx], bit_pos) as usize
        } else {
            let byte_idx = (1 + (a_index - 1) / 4) as usize;
            let bit_pos  = (3 - (a_index - 1) % 4) as u32;
            ncbi2na_unpack_base(a[byte_idx], bit_pos) as usize
        };
        let matrix_row = &mat[a_base_pair];

        // b_ptr starting position for this row.
        // C: `if (reverse_sequence) b_ptr = &B[N - first_b_index];`
        //    `else                  b_ptr = B + first_b_index;`
        let mut b_ptr_idx: i32 = if reverse_sequence {
            n - first_b_index
        } else {
            first_b_index
        };

        let mut score: i32      = MININT;
        let mut score_gap_row   = MININT;
        let mut last_b_index    = first_b_index;

        for b_index in first_b_index..b_size {
            b_ptr_idx += b_increment;

            let mut score_gap_col = gap_align.dp_mem[b_index as usize].best_gap;
            // C: `next_score = score_array[b_index].best + matrix_row[*b_ptr];`
            let q_base = b[b_ptr_idx as usize] as usize;
            let next_score = gap_align.dp_mem[b_index as usize].best + matrix_row[q_base];

            if score < score_gap_col {
                score = score_gap_col;
            }
            if score < score_gap_row {
                score = score_gap_row;
            }

            if best_score - score > x_dropoff {
                if b_index == first_b_index {
                    first_b_index += 1;
                } else {
                    gap_align.dp_mem[b_index as usize].best = MININT;
                }
            } else {
                last_b_index = b_index;
                if score > best_score {
                    best_score = score;
                    *a_offset  = a_index;
                    *b_offset  = b_index;
                }
                score_gap_row -= gap_extend;
                score_gap_col -= gap_extend;
                gap_align.dp_mem[b_index as usize].best_gap =
                    (score - gap_open_extend).max(score_gap_col);
                score_gap_row = (score - gap_open_extend).max(score_gap_row);
                gap_align.dp_mem[b_index as usize].best = score;
            }

            score = next_score;
        }

        if first_b_index == b_size {
            break;
        }

        // Grow dp_mem if needed.
        if last_b_index + num_extra_cells + 3 >= gap_align.dp_mem.len() as i32 {
            let new_alloc = (last_b_index + num_extra_cells + 100)
                .max(2 * gap_align.dp_mem.len() as i32);
            gap_align.dp_mem.resize(new_alloc as usize, BlastGapDP::default());
        }

        if last_b_index < b_size - 1 {
            b_size = last_b_index + 1;
        } else {
            // Extend b_size while score_gap_row stays above threshold.
            while score_gap_row >= best_score - x_dropoff && b_size <= n {
                gap_align.dp_mem[b_size as usize].best     = score_gap_row;
                gap_align.dp_mem[b_size as usize].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                b_size += 1;
            }
        }

        if b_size <= n {
            gap_align.dp_mem[b_size as usize].best     = MININT;
            gap_align.dp_mem[b_size as usize].best_gap = MININT;
            b_size += 1;
        }
    }

    best_score
}

// ── s_BlastDynProgNtGappedAlignment (blast_gapalign.c lines 2948–3017) ───────

/// Preliminary gapped alignment (no traceback) for nucleotides.
///
/// Mirrors `static Int2 s_BlastDynProgNtGappedAlignment(...)`.
///
/// Parameters:
///   query_seq    – BLASTNA query, starting at base-0 (C: query_blk->sequence).
///   query_len    – total query length.
///   subject_seq  – NCBI2NA packed subject, starting at byte-0
///                  (C: subject_blk->sequence).
///   subject_len  – total subject length in bases.
///   q_off        – query offset of the seed (from init_hsp->offsets.q_off).
///   s_off        – subject offset of the seed (from init_hsp->offsets.s_off).
///   gap_align    – workspace; receives query_start/stop, subject_start/stop, score.
///   score_params – gap penalties.
///   matrix       – score matrix.
///
/// Returns 0 on success, -1 on error.
pub fn s_blast_dynprog_nt_gapped_alignment(
    query_seq:   &[u8],
    query_len:   i32,
    subject_seq: &[u8],
    subject_len: i32,
    q_off:       i32,
    s_off:       i32,
    gap_align:   &mut BlastGapAlignStruct,
    score_params: &ScoringParams,
    matrix:      &ScoreMatrix,
) -> i32 {
    let x_dropoff = gap_align.gap_x_dropoff;

    // C lines 2972–2984: round s_off up to next 4-base boundary.
    // "If subject offset is not at the start of a full byte, shift the
    //  alignment start to the next multiple of 4 subject letters."
    // C: `offset_adjustment = COMPRESSION_RATIO - (s_off % COMPRESSION_RATIO)`
    // This is always 1..4, never 0 (when s_off%4==0, adjustment=4, not 0).
    let offset_adjustment = 4 - (s_off % 4); // always 1..4

    let mut q_length = q_off + offset_adjustment;
    let mut s_length = s_off + offset_adjustment;

    // C lines 2980–2984: prevent pivot from being past the end.
    if q_length > query_len || s_length > subject_len {
        q_length -= 4;
        s_length -= 4;
    }

    // Left extension (reverse_sequence = TRUE).
    // C: `score_left = s_BlastAlignPackedNucl(query, subject, q_length, s_length, ...)`
    let mut private_q_start: i32 = 0;
    let mut private_s_start: i32 = 0;

    let score_left = s_blast_align_packed_nucl(
        query_seq,   // B = query (BLASTNA)
        subject_seq, // A = subject (NCBI2NA)
        q_length,    // N
        s_length,    // M
        &mut private_q_start,
        &mut private_s_start,
        gap_align,
        score_params,
        matrix,
        true,        // reverse_sequence
        x_dropoff,
    );
    if score_left < 0 {
        return -1;
    }
    gap_align.query_start   = q_length - private_q_start;
    gap_align.subject_start = s_length - private_s_start;

    // Right extension (reverse_sequence = FALSE).
    // C: `s_BlastAlignPackedNucl(query+q_length-1, subject+(s_length+3)/4-1, ...)`
    let score_right;
    if q_length < query_len && s_length < subject_len {
        // Subject pointer for right extension: byte index (s_length+3)/4 - 1.
        // Since s_length is guaranteed to be a multiple of 4 after offset_adjustment,
        // (s_length+3)/4 = s_length/4.
        let s_byte_start = ((s_length + 3) / 4 - 1) as usize;
        let mut qs: i32 = 0;
        let mut ss: i32 = 0;
        score_right = s_blast_align_packed_nucl(
            &query_seq[(q_length - 1) as usize..],
            &subject_seq[s_byte_start..],
            query_len   - q_length,
            subject_len - s_length,
            &mut qs,
            &mut ss,
            gap_align,
            score_params,
            matrix,
            false, // reverse_sequence
            x_dropoff,
        );
        if score_right < 0 {
            return -1;
        }
        gap_align.query_stop   = qs + q_length;
        gap_align.subject_stop = ss + s_length;
    } else {
        score_right = 0;
        gap_align.query_stop   = q_length;
        gap_align.subject_stop = s_length;
    }

    gap_align.score = score_right + score_left;
    0
}

// ── ALIGN_EX (blast_gapalign.c lines 374–767) ─────────────────────────────────

/// Gapped extension WITH traceback for nucleotides.
///
/// Mirrors `static Int4 ALIGN_EX(const Uint1* A, const Uint1* B, ...)`.
///
/// Both A (subject) and B (query) are BLASTNA unpacked (1-byte-per-base).
/// This is the traceback phase function; the matrix lookup is symmetric.
///
/// Parameters (matching C signature):
///   a               – subject slice starting at extension pivot.
///   b               – query slice starting at extension pivot.
///   m               – max extension in subject (from pivot).
///   n               – max extension in query (from pivot).
///   a_offset        – output: bases consumed in subject.
///   b_offset        – output: bases consumed in query.
///   edit_block      – receives the edit operations for this extension.
///   gap_align       – workspace (dp_mem).
///   score_params     – gap penalties.
///   query_offset    – absolute query offset of the pivot (used for PSSM;
///                     for non-PSSM searches pass 0 and ignore).
///   reversed        – is A the reversed sequence? (for matrix row selection)
///   reverse_sequence – TRUE → left extension (read A/B backward).
///   x_dropoff       – x-drop threshold.
///
/// Returns the best alignment score.
/// In NCBI's C code, ALIGN_EX receives the full `gap_align` struct.  In Rust we
/// pass `dp_mem` separately to avoid the borrow conflict that arises when
/// `edit_block` is also a field of `gap_align` (fwd/rev_prelim_tback).
#[allow(clippy::too_many_arguments)]
pub fn align_ex(
    a: &[u8],             // subject (BLASTNA)
    b: &[u8],             // query (BLASTNA)
    m: i32,               // max subject extension
    n: i32,               // max query extension
    a_offset: &mut i32,
    b_offset: &mut i32,
    edit_block: &mut GapPrelimEditBlock,
    dp_mem: &mut Vec<BlastGapDP>,   // gap_align->dp_mem
    score_params: &ScoringParams,
    matrix: &ScoreMatrix,
    _query_offset: i32,
    _reversed: bool,
    reverse_sequence: bool,
    x_dropoff: i32,
) -> i32 {
    let gap_open        = score_params.gap_open;
    let gap_extend      = score_params.gap_extend;
    let gap_open_extend = gap_open + gap_extend;
    let mat             = &matrix.scores;

    *a_offset = 0;
    *b_offset = 0;

    let mut x_dropoff = x_dropoff;
    if x_dropoff < gap_open_extend {
        x_dropoff = gap_open_extend;
    }

    if n <= 0 || m <= 0 {
        return 0;
    }

    // ── Traceback state allocation (mirrors s_GapPurgeState / edit_script) ─

    // In NCBI, the traceback is stored in a 2D array: edit_script[a_index][b_index].
    // Each row is allocated separately from a pool (GapStateArrayStruct).
    // Here we use a Vec of Vecs; state_rows[a] is the script row for subject
    // position a.  edit_start_offset[a] is the b-index corresponding to
    // state_rows[a][0].
    let mut state_rows: Vec<Vec<u8>> = Vec::with_capacity(100);
    let mut edit_start_offset: Vec<i32> = Vec::with_capacity(100);

    // Allocate / grow DP scratch.
    let num_extra_cells = if gap_extend > 0 {
        x_dropoff / gap_extend + 3
    } else {
        n + 3
    };

    if num_extra_cells > dp_mem.len() as i32 {
        let new_alloc = (num_extra_cells + 100).max(2 * dp_mem.len() as i32);
        dp_mem.resize(new_alloc as usize, BlastGapDP::default());
    }

    // Initialise row 0.
    let mut score: i32 = -gap_open_extend;
    dp_mem[0].best     = 0;
    dp_mem[0].best_gap = -gap_open_extend;

    let mut b_size = 1i32;
    while b_size <= n {
        if score < -x_dropoff {
            break;
        }
        dp_mem[b_size as usize].best     = score;
        dp_mem[b_size as usize].best_gap = score - gap_open_extend;
        score -= gap_extend;
        b_size += 1;
    }

    // Allocate row 0: all cells are SCRIPT_GAP_IN_A (gaps in subject along initial row).
    // C: for i=1..N: edit_script_row[i] = SCRIPT_GAP_IN_A; (position 0 is 0-initialised = same)
    {
        // Allocate num_extra_cells extra so the extension loop can write into this row.
        let row_cap = (b_size + num_extra_cells) as usize;
        let row = vec![SCRIPT_GAP_IN_A; row_cap];
        state_rows.push(row);
        edit_start_offset.push(0);
    }

    let mut best_score = 0i32;
    let mut first_b_index = 0i32;
    let b_increment: i32 = if reverse_sequence { -1 } else { 1 };

    let mut best_a_index = 0i32;
    let mut best_b_index = 0i32;

    for a_index in 1i32..=m {
        // Get subject base for this row.
        let a_ptr: i32 = if reverse_sequence { m - a_index } else { a_index };
        let sub_base = a[a_ptr as usize] as usize;
        let matrix_row = &mat[sub_base];

        let mut b_ptr_idx: i32 = if reverse_sequence {
            n - first_b_index
        } else {
            first_b_index
        };

        // Allocate edit-script row for a_index.
        // C allocates `b_size - first_b_index + num_extra_cells` cells so the
        // extension loop (after the inner loop) can write SCRIPT_GAP_IN_A into
        // the new b_size slots without growing the row.
        let row_start = first_b_index;
        let row_cap   = (b_size - first_b_index + num_extra_cells) as usize;
        let mut row   = vec![SCRIPT_GAP_IN_A; row_cap];
        edit_start_offset.push(row_start);

        let mut score         = MININT;
        let mut score_gap_row = MININT;
        let mut last_b_index  = first_b_index;

        for b_index in first_b_index..b_size {
            b_ptr_idx += b_increment;

            // ri is the index into `row` for this b_index.
            let ri = (b_index - row_start) as usize;

            let mut score_gap_col = dp_mem[b_index as usize].best_gap;
            let q_base            = b[b_ptr_idx as usize] as usize;
            let next_score        = dp_mem[b_index as usize].best + matrix_row[q_base];

            // Determine best incoming move (mirrors C lines 588–599).
            let mut script     = SCRIPT_SUB;
            let script_col: u8 = SCRIPT_EXTEND_GAP_B; // 0x40 — gap-in-B continuation flag
            let script_row: u8 = SCRIPT_EXTEND_GAP_A; // 0x10 — gap-in-A continuation flag

            if score < score_gap_col {
                script = SCRIPT_GAP_IN_B;
                score  = score_gap_col;
            }
            if score < score_gap_row {
                script = SCRIPT_GAP_IN_A;
                score  = score_gap_row;
            }

            if best_score - score > x_dropoff {
                // Cell dropped by x-drop.
                if b_index == first_b_index {
                    first_b_index += 1;
                } else {
                    dp_mem[b_index as usize].best = MININT;
                }
                // row[ri] is already SCRIPT_GAP_IN_A from initialisation.
            } else {
                last_b_index = b_index;
                if score > best_score {
                    best_score   = score;
                    best_a_index = a_index;
                    best_b_index = b_index;
                    *a_offset    = a_index;
                    *b_offset    = b_index;
                }

                score_gap_row -= gap_extend;
                score_gap_col -= gap_extend;

                // Update best_gap for column gaps (gap in B / subject).
                // When score_gap_col wins, we are *extending* a gap in B:
                // set the SCRIPT_EXTEND_GAP_B flag in the stored script.
                if score_gap_col < score - gap_open_extend {
                    dp_mem[b_index as usize].best_gap = score - gap_open_extend;
                } else {
                    dp_mem[b_index as usize].best_gap = score_gap_col;
                    script += script_col; // flag: gap-in-B continues
                }

                // Update best_gap for row gaps (gap in A / query).
                // When score_gap_row wins, we are *extending* a gap in A:
                // set the SCRIPT_EXTEND_GAP_A flag in the stored script.
                if score_gap_row < score - gap_open_extend {
                    score_gap_row = score - gap_open_extend;
                } else {
                    script += script_row; // flag: gap-in-A continues
                }

                dp_mem[b_index as usize].best = score;
                row[ri] = script;
            }

            score = next_score;
        }

        state_rows.push(row);

        if first_b_index == b_size {
            break;
        }

        // Grow dp_mem if needed.
        if last_b_index + num_extra_cells + 3 >= dp_mem.len() as i32 {
            let new_alloc = (last_b_index + num_extra_cells + 100)
                .max(2 * dp_mem.len() as i32);
            dp_mem.resize(new_alloc as usize, BlastGapDP::default());
        }

        if last_b_index < b_size - 1 {
            b_size = last_b_index + 1;
        } else {
            // Extend b_size for additional gap-in-A cells.
            // C: edit_script_row[b_size] = SCRIPT_GAP_IN_A (already in row from init).
            while score_gap_row >= best_score - x_dropoff && b_size <= n {
                dp_mem[b_size as usize].best     = score_gap_row;
                dp_mem[b_size as usize].best_gap = score_gap_row - gap_open_extend;
                score_gap_row -= gap_extend;
                // row[(b_size - row_start) as usize] = SCRIPT_GAP_IN_A; // already initialised
                b_size += 1;
            }
        }

        if b_size <= n {
            dp_mem[b_size as usize].best     = MININT;
            dp_mem[b_size as usize].best_gap = MININT;
            b_size += 1;
        }
    }

    // ── Traceback (mirrors C lines 682–727) ───────────────────────────────
    // Walk backward from (best_a_index, best_b_index) to (0, 0).
    // `prev_script` holds the operation taken in the *previous* step; it
    // determines which extension flag to consult in the current cell.
    edit_block.reset();
    let mut ai          = best_a_index;
    let mut bi          = best_b_index;
    let mut prev_script = SCRIPT_SUB; // C initialises `script = SCRIPT_SUB`

    while ai > 0 || bi > 0 {
        // Look up the traceback byte for the current cell.
        let row_idx = ai as usize;
        let row_off = edit_start_offset[row_idx];
        let ri_signed = bi - row_off;

        let next_script =
            if ri_signed >= 0 && (ri_signed as usize) < state_rows[row_idx].len() {
                state_rows[row_idx][ri_signed as usize]
            } else {
                // Shouldn't happen for a valid alignment; use SCRIPT_GAP_IN_A
                // as a safe fallback (row 0 / boundary case).
                SCRIPT_GAP_IN_A
            };

        // Determine the operation for this step.
        // C switch(script) { case SCRIPT_GAP_IN_A: ...; case SCRIPT_GAP_IN_B: ...; default: ... }
        let script = match prev_script {
            SCRIPT_GAP_IN_A => {
                // We were extending a gap in A (deletion).  Check the
                // SCRIPT_EXTEND_GAP_A flag to see if it continues.
                let s = next_script & SCRIPT_OP_MASK;
                if next_script & SCRIPT_EXTEND_GAP_A != 0 { SCRIPT_GAP_IN_A } else { s }
            }
            SCRIPT_GAP_IN_B => {
                // We were extending a gap in B (insertion).  Check the
                // SCRIPT_EXTEND_GAP_B flag to see if it continues.
                let s = next_script & SCRIPT_OP_MASK;
                if next_script & SCRIPT_EXTEND_GAP_B != 0 { SCRIPT_GAP_IN_B } else { s }
            }
            _ => next_script & SCRIPT_OP_MASK,
        };
        prev_script = script;

        // Move backward and emit the edit operation.
        // SCRIPT_GAP_IN_A: gap in subject A → consume one query (B) base → bi--.
        // SCRIPT_GAP_IN_B: gap in query B  → consume one subject (A) base → ai--.
        // SCRIPT_SUB:      substitution    → consume both → ai--, bi--.
        if script == SCRIPT_GAP_IN_A {
            bi -= 1;
            edit_block.add(GapAlignOp::Del, 1);
        } else if script == SCRIPT_GAP_IN_B {
            ai -= 1;
            edit_block.add(GapAlignOp::Ins, 1);
        } else {
            ai -= 1;
            bi -= 1;
            edit_block.add(GapAlignOp::Sub, 1);
        }
    }

    best_score
}

// ── BLAST_GappedAlignmentWithTraceback (blast_gapalign.c lines 4595–4700) ────

/// Full gapped alignment WITH traceback for blastn.
///
/// Mirrors `Int2 BLAST_GappedAlignmentWithTraceback(...)`.
///
/// Both query and subject are BLASTNA unpacked (1-byte-per-base).
/// Uses the larger xdrop_gap_final x-drop.
///
/// Returns 0 on success; stores results in gap_align.
pub fn blast_gapped_alignment_with_traceback(
    query:          &[u8],   // BLASTNA
    subject:        &[u8],   // BLASTNA
    gap_align:      &mut BlastGapAlignStruct,
    score_params:   &ScoringParams,
    matrix:         &ScoreMatrix,
    q_start:        i32,
    s_start:        i32,
    query_length:   i32,
    subject_length: i32,
) -> i32 {
    let x_dropoff = gap_align.gap_x_dropoff;
    gap_align.fwd_prelim_tback.reset();
    gap_align.rev_prelim_tback.reset();

    // Left extension: includes the starting point [q_start, s_start].
    // C: `score_left = ALIGN_EX(query, subject, q_start+1, s_start+1, ...
    //                            rev_prelim_tback, ..., reversed=FALSE,
    //                            reverse_sequence=TRUE)`
    let mut priv_q_len = 0i32;
    let mut priv_s_len = 0i32;
    let score_left = align_ex(
        subject,               // A
        query,                 // B
        q_start + 1,           // M (subject extent; note C passes q/s swapped here)
        s_start + 1,           // N
        &mut priv_q_len,
        &mut priv_s_len,
        &mut gap_align.rev_prelim_tback,
        &mut gap_align.dp_mem,
        score_params,
        matrix,
        q_start,
        false,
        true,  // reverse_sequence
        x_dropoff,
    );
    gap_align.query_start   = q_start - priv_q_len + 1;
    gap_align.subject_start = s_start - priv_s_len + 1;

    // Right extension: does NOT include the starting point.
    // C: `score_right = ALIGN_EX(query+q_start, subject+s_start,
    //                             query_length-q_start, subject_length-s_start, ...
    //                             fwd_prelim_tback, ..., reversed=FALSE,
    //                             reverse_sequence=FALSE)`
    let score_right;
    if q_start < query_length && s_start < subject_length {
        let mut qs = 0i32;
        let mut ss = 0i32;
        score_right = align_ex(
            &subject[s_start as usize..],
            &query[q_start as usize..],
            subject_length - s_start,
            query_length   - q_start,
            &mut ss,
            &mut qs,
            &mut gap_align.fwd_prelim_tback,
            &mut gap_align.dp_mem,
            score_params,
            matrix,
            q_start,
            false,
            false, // forward
            x_dropoff,
        );
        gap_align.query_stop   = qs + q_start;
        gap_align.subject_stop = ss + s_start;
    } else {
        score_right = 0;
        gap_align.query_stop   = q_start;
        gap_align.subject_stop = s_start;
    }

    gap_align.score = score_left + score_right;
    0
}

// ── BlastGetStartForGappedAlignmentNucl (blast_gapalign.c lines 3322–3390) ───

/// Find an improved seed point within the preliminary HSP boundaries.
///
/// Mirrors `void BlastGetStartForGappedAlignmentNucl(const Uint1* query,
///   const Uint1* subject, BlastHSP* hsp)`.
///
/// Both query and subject are BLASTNA unpacked (1-byte-per-base).
///
/// On entry, `hsp.query.gapped_start` and `hsp.subject.gapped_start` hold the
/// preliminary seed (from the prelim phase).  This function may update them to
/// a position with a longer identity run.
pub fn blast_get_start_for_gapped_alignment_nucl(
    query:   &[u8],  // BLASTNA query,   index 0 = first base
    subject: &[u8],  // BLASTNA subject, index 0 = first base
    hsp:     &mut BlastHsp,
) {
    const HSP_MAX_IDENT_RUN: i32 = 10;

    let q_gapped_start = hsp.query.gapped_start;
    let s_gapped_start = hsp.subject.gapped_start;

    // C: `offset = MIN(hsp->subject.gapped_start - hsp->subject.offset,
    //                  hsp->query.gapped_start - hsp->query.offset)`
    let offset = (s_gapped_start - hsp.subject.offset)
        .min(q_gapped_start - hsp.query.offset);

    // ── First check if the old value is already in a long identity run ────
    let q_len = hsp.query.end;
    let mut q = q_gapped_start;
    let mut s = s_gapped_start;
    let mut score: i32 = -1;

    // Forward scan from gapped_start.
    while q < q_len && query[q as usize] == subject[s as usize] {
        score += 1;
        q     += 1;
        s     += 1;
        if score > HSP_MAX_IDENT_RUN {
            return; // already good
        }
    }
    // Backward scan from gapped_start.
    let mut q = q_gapped_start;
    let mut s = s_gapped_start;
    while q >= 0 && s >= 0 && query[q as usize] == subject[s as usize] {
        score += 1;
        if score > HSP_MAX_IDENT_RUN {
            return;
        }
        if q == 0 || s == 0 { break; }
        q -= 1;
        s -= 1;
    }

    // ── Need to find a better starting point ─────────────────────────────
    // C: `hspMaxIdentRun *= 1.5` (integer: 10 → 15)
    let hsp_max_ident_run = (HSP_MAX_IDENT_RUN as f64 * 1.5) as i32;

    let q_start = q_gapped_start - offset;
    let s_start = s_gapped_start - offset;
    let scan_len = (hsp.subject.end - s_start).min(hsp.query.end - q_start);

    let mut max_score  = 0i32;
    let mut max_offset = q_start;
    let mut score      = 0i32;
    let mut match_     = false;
    let mut prev_match = false;

    for index in q_start..q_start + scan_len {
        let qi = index as usize;
        let si = (s_start + (index - q_start)) as usize;
        match_ = query[qi] == subject[si];

        if match_ != prev_match {
            prev_match = match_;
            if match_ {
                score = 1;
            } else if score > max_score {
                max_score  = score;
                max_offset = index - score / 2;
            }
        } else if match_ {
            score += 1;
            if score > hsp_max_ident_run {
                max_offset = index - hsp_max_ident_run / 2;
                hsp.query.gapped_start   = max_offset;
                hsp.subject.gapped_start = max_offset + s_start - q_start;
                return;
            }
        }
    }

    if match_ && score > max_score {
        max_score  = score;
        max_offset = (q_start + scan_len) - score / 2;
    }

    if max_score > 0 {
        hsp.query.gapped_start   = max_offset;
        hsp.subject.gapped_start = max_offset + s_start - q_start;
    }
}
