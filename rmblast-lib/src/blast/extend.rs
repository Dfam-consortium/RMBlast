// blast/extend.rs — faithful port of nucleotide ungapped extension.
//
// Primary source: na_ungapped.c (NCBI BLAST 2.17.0+)
//
// Functions ported:
//   s_NuclUngappedExtendExact   (na_ungapped.c lines 149–243)
//   blast_save_initial_hit      (blast_extend.c lines 325–357)

use super::types::{
    BlastInitHitList, BlastInitHsp, BlastOffsetPair, BlastUngappedData, COMPRESSION_RATIO,
};
use super::util::ncbi2na_unpack_base;
use crate::matrix::ScoreMatrix;

// ── s_NuclUngappedExtendExact (na_ungapped.c lines 149–243) ─────────────────

/// Matrix-based ungapped extension from a word hit.
///
/// This is the RMBlast path (matrix_only_scoring set, or word_length < 11).
/// Mirrors `static void s_NuclUngappedExtendExact(...)`.
///
/// # Parameters
/// * `query_blastna` — full BLASTNA slice **including** the sentinel byte at
///   index 0.  Base i (0-based) lives at index i+1.  This sentinel mimics
///   NCBI's `query_start[0] = 0` which bounds left extension automatically.
/// * `subject_seq`   — NCBI2NA packed bytes.  Byte k covers bases 4k..4k+3
///   with base 4k in bits 7:6 (N=3) and base 4k+3 in bits 1:0 (N=0).
/// * `query_len`     — query length in bases (not counting the sentinel).
/// * `subject_len`   — subject length in bases.
/// * `q_off`         — 0-based query position of the word hit.
/// * `s_off`         — 0-based subject position of the word hit (in bases).
/// * `x_dropoff`     — **negative** x-drop threshold (C passes −cutoff).
/// * `ungapped_data` — filled with extension result.
pub fn s_nucl_ungapped_extend_exact(
    query_blastna: &[u8],
    subject_seq:   &[u8],
    query_len:     i32,
    subject_len:   i32,
    matrix:        &ScoreMatrix,
    q_off:         i32,
    s_off:         i32,
    x_dropoff:     i32,   // negative
    ungapped_data: &mut BlastUngappedData,
) {
    let mat = &matrix.scores;

    // In C: q = query->sequence + q_off  (1 byte/base, sentinel at [-1])
    // In Rust: query_blastna[0] = sentinel; base q_off is at index q_off+1.
    let q_start_ptr: i32 = q_off + 1; // blastna index for base q_off

    // Subject: current byte index and bit position within that byte.
    let mut s_ptr: i32 = s_off / COMPRESSION_RATIO;
    let mut base:  i32 = 3 - (s_off % COMPRESSION_RATIO);

    // Subject left boundary.
    // C: if (q_off < s_off)
    //      start = subject0 + (s_off-q_off)/4;  remainder = 3-((s_off-q_off)%4)
    //    else
    //      start = subject0;                     remainder = 3
    let (start_s_ptr, start_rem): (i32, i32) = if q_off < s_off {
        let back = s_off - q_off;
        (back / COMPRESSION_RATIO, 3 - (back % COMPRESSION_RATIO))
    } else {
        (0, 3)
    };

    let mut score:  i32 = 0;
    let mut sum:    i32 = 0;
    let mut q_ptr:  i32 = q_start_ptr;
    let mut q_beg:  i32 = q_start_ptr; // leftmost good query index
    let mut q_end:  i32 = q_start_ptr; // one-past rightmost good query index

    // ── Left extension ───────────────────────────────────────────────────────
    // C: while ((s > start) || (s == start && base < remainder))
    //      if (base == 3) { s--; base = 0; } else { ++base; }
    //      sum += matrix[*--q][NCBI2NA_UNPACK_BASE(*s, base)];
    loop {
        // Subject boundary check (same logic as C).
        if !(s_ptr > start_s_ptr || (s_ptr == start_s_ptr && base < start_rem)) {
            break;
        }

        // Step subject pointer left to the previous base.
        if base == 3 { s_ptr -= 1; base = 0; } else { base += 1; }

        // Step query left (pre-decrement in C: *--q).
        q_ptr -= 1;
        // Safety guard: should not underflow with correct boundary logic,
        // but prevents Rust panic if they diverge by one step.
        if q_ptr < 0 { break; }

        let ch     = subject_seq[s_ptr as usize];
        let s_base = ncbi2na_unpack_base(ch, base as u32) as usize;
        let q_base = query_blastna[q_ptr as usize] as usize; // [0] = sentinel OK

        sum += mat[q_base][s_base];

        if sum > 0 {
            q_beg   = q_ptr;
            score  += sum;
            sum     = 0;
        } else if sum < x_dropoff {
            break;
        }
    }

    // C: ungapped_data->q_start = q_beg - query->sequence
    // Our q_beg is a blastna index (1-based), so subtract 1 for 0-based.
    ungapped_data.q_start = q_beg - 1;
    ungapped_data.s_start = s_off - (q_off - ungapped_data.q_start);

    // ── Right extension ──────────────────────────────────────────────────────
    // C: q_avail = query->length - q_off;  s_avail = subject->length - s_off
    //    if (q_avail < s_avail)
    //      sf = subject0 + (s_off + q_avail) / 4;
    //      remainder = 3 - ((s_off + q_avail) % 4);
    //    else
    //      sf = subject0 + subject->length / 4;
    //      remainder = 3 - (subject->length % 4);
    let q_avail   = query_len   - q_off;
    let s_avail   = subject_len - s_off;
    let s_end_off = s_off + q_avail.min(s_avail); // exclusive end in bases
    let sf:    i32 = s_end_off / COMPRESSION_RATIO;
    let rem:   i32 = 3 - (s_end_off % COMPRESSION_RATIO);

    // Restart from the hit point.
    q_ptr = q_start_ptr;
    s_ptr = s_off / COMPRESSION_RATIO;
    base  = 3 - (s_off % COMPRESSION_RATIO);
    sum   = 0;

    // Adaptive x-drop: X_current = (-score > X) ? -score : X.
    // Tightens the threshold when the current score is low (early extension).
    let mut x_current: i32 = x_dropoff;

    // C: while (s < sf || (s == sf && base > remainder))
    //      sum += matrix[*q++][NCBI2NA_UNPACK_BASE(*s, base)];
    //      if (sum > 0) { q_end = q; score += sum; X_current = …; sum = 0; }
    //      if (base == 0) { base = 3; s++; } else { base--; }
    loop {
        // Right boundary check.
        if !(s_ptr < sf || (s_ptr == sf && base > rem)) { break; }

        let ch     = subject_seq[s_ptr as usize];
        let s_base = ncbi2na_unpack_base(ch, base as u32) as usize;
        let q_base = query_blastna[q_ptr as usize] as usize;
        q_ptr     += 1; // post-increment (C: *q++)

        sum += mat[q_base][s_base];

        if sum > 0 {
            q_end     = q_ptr;
            score    += sum;
            x_current = (-score).max(x_dropoff); // adaptive threshold
            sum       = 0;
        } else if sum < x_current {
            break;
        }

        // Advance to the next subject base.
        if base == 0 { base = 3; s_ptr += 1; } else { base -= 1; }
    }

    // C: ungapped_data->length = q_end - q_beg;  (both are pointer differences)
    // Our q_end and q_beg are blastna indices; their difference = base count.
    ungapped_data.length = q_end - q_beg;
    ungapped_data.score  = score;
}

// ── BLAST_SaveInitialHit (blast_extend.c lines 325–357) ──────────────────────

/// Save an initial HSP into the hit list.
///
/// Mirrors `Boolean BLAST_SaveInitialHit(BlastInitHitList*, Int4, Int4,
///   BlastUngappedData*)`.
///
/// In C, `ungapped_data` is a heap-allocated pointer (may be NULL).
/// Here we pass `Option<BlastUngappedData>` by value.
///
/// Returns `true` if the hit was saved, `false` if the list is full and
/// reallocation is disabled (not applicable for Vec — always returns true).
pub fn blast_save_initial_hit(
    init_hitlist:  &mut BlastInitHitList,
    q_off:         i32,
    s_off:         i32,
    ungapped_data: Option<BlastUngappedData>,
) -> bool {
    let (has_ungapped, ud) = match ungapped_data {
        Some(d) => (true,  d),
        None    => (false, BlastUngappedData::default()),
    };
    init_hitlist.init_hsp_array.push(BlastInitHsp {
        offsets: BlastOffsetPair { q_off: q_off as u32, s_off: s_off as u32 },
        ungapped_data: ud,
        has_ungapped,
    });
    true
}
