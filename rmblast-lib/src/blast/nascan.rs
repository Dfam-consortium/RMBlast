// blast/nascan.rs — faithful port of nucleotide subject scanning.
//
// Primary source: blast_nascan.c (NCBI BLAST 2.17.0+)
//
// Functions ported:
//   s_BlastLookupGetNumHits   (blast_nascan.c lines 41–49)
//   s_BlastLookupRetrieve     (blast_nascan.c lines 58–82)
//   s_BlastNaScanSubject_8_4  (blast_nascan.c lines 95–134)
//   s_BlastNaScanSubject_Any  (blast_nascan.c lines 147–278, extended to lw 9–12)
//   s_NaChooseScanSubject     (blast_nascan.c lines 284–293)
//
// Extension beyond NCBI source: s_BlastNaScanSubject_Any is extended to handle
// lut_word_length 9–12 using 3-byte (aligned) or 3/4-byte (unaligned) reads.
// NCBI uses eMBLookupTable for those widths; we extend eNaLookupTable instead.

use super::types::{BlastOffsetPair, COMPRESSION_RATIO};
use super::nalookup::{BlastNaLookupTable, NaLookupPayload, NA_HITS_PER_CELL, pv_test};

// ── s_BlastLookupGetNumHits (blast_nascan.c lines 41–49) ─────────────────────

/// Return the number of query hits for the given lookup index.
/// Returns 0 if the PV bit is not set (no hits in that cell).
#[inline]
fn s_blast_lookup_get_num_hits(lookup: &BlastNaLookupTable, index: i32) -> i32 {
    if pv_test(&lookup.pv, index) {
        lookup.thick_backbone[index as usize].num_used
    } else {
        0
    }
}

// ── s_BlastLookupRetrieve (blast_nascan.c lines 58–82) ───────────────────────

/// Copy query offsets from the lookup cell into `offset_pairs`.
/// Fills `offset_pairs[0..num_used]`; caller must ensure capacity.
///
/// Mirrors `static NCBI_INLINE void s_BlastLookupRetrieve(...)`.
#[inline]
fn s_blast_lookup_retrieve(
    lookup:       &BlastNaLookupTable,
    index:        i32,
    offset_pairs: &mut [BlastOffsetPair],
    s_off:        i32,
) {
    let cell      = &lookup.thick_backbone[index as usize];
    let num_hits  = cell.num_used as usize;

    let query_offsets: &[i32] = match &cell.payload {
        NaLookupPayload::Entries(e) => &e[..num_hits.min(NA_HITS_PER_CELL)],
        NaLookupPayload::OverflowCursor(cur) => {
            let start = *cur as usize;
            &lookup.overflow[start..start + num_hits]
        }
    };

    for i in 0..num_hits {
        offset_pairs[i].q_off = query_offsets[i] as u32;
        offset_pairs[i].s_off = s_off as u32;
    }
}

// ── s_BlastNaScanSubject_8_4 (blast_nascan.c lines 95–134) ───────────────────

/// Scan for lut_word_length==8, scan_step==4 (aligned stride: 1 byte/step).
///
/// Subject must be NCBI2NA packed bytes with at least one sentinel byte at the
/// end (same contract as NCBI C code, which relies on the packing allocation).
///
/// `scan_range` — `[start_base, end_base]` inclusive; on return `[0]` is the
///               next unscanned base position.
/// Returns total hits written into `offset_pairs[0..return_value]`.
pub fn s_blast_na_scan_subject_8_4(
    lookup:       &BlastNaLookupTable,
    subject:      &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    max_hits:     i32,
    scan_range:   &mut [i32; 2],
) -> i32 {
    let mut s_byte  = (scan_range[0] / COMPRESSION_RATIO) as usize;
    let s_end       = (scan_range[1] / COMPRESSION_RATIO) as usize;
    let mut total   = 0i32;

    while s_byte <= s_end {
        // Two-byte index: lut_word_length == 8 ↔ 16 bits = 2 subject bytes.
        let index = ((subject[s_byte] as i32) << 8) | (subject[s_byte + 1] as i32);

        let num_hits = s_blast_lookup_get_num_hits(lookup, index);
        if num_hits == 0 {
            s_byte += 1;
            continue;
        }
        if num_hits > max_hits - total {
            break;
        }

        s_blast_lookup_retrieve(
            lookup,
            index,
            &mut offset_pairs[total as usize..],
            (s_byte as i32) * COMPRESSION_RATIO,
        );
        total  += num_hits;
        s_byte += 1;
    }

    scan_range[0] = (s_byte as i32) * COMPRESSION_RATIO;
    total
}

// ── s_BlastNaScanSubject_Any (blast_nascan.c lines 147–278) ──────────────────

/// General scanner for lut_word_length 4–12 at arbitrary stride.
///
/// Handles three sub-cases based on lut_word_length:
///  • lw ≤ 5: two packed bytes per lookup (word always fits in 8 bits of index).
///  • lw 6–8: two bytes (aligned stride) or three bytes (unaligned stride).
///  • lw 9–12: three bytes (aligned stride); three or four bytes (unaligned
///    stride, four when s_off%4 + lut_word_length > 12).
///
/// `scan_range` — `[start_base, end_base]` inclusive; updated to stopping pos.
/// Returns total hits written.
pub fn s_blast_na_scan_subject_any(
    lookup:       &BlastNaLookupTable,
    subject:      &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    max_hits:     i32,
    scan_range:   &mut [i32; 2],
) -> i32 {
    let mask            = lookup.mask;
    let scan_step       = lookup.scan_step;
    let lut_word_length = lookup.lut_word_length;
    let mut total       = 0i32;

    if lut_word_length > 8 {
        // lut_word_length 9–12: a 3-byte window covers 12 bases (sufficient for
        // aligned stride and for unaligned positions where s%4 + lw ≤ 12).
        // When s%4 + lw > 12 the word straddles a 4th byte; use a u32 read.

        if scan_step % COMPRESSION_RATIO == 0 {
            // Aligned: every word starts on a 4-base boundary.
            // 3 bytes = 12 bases; right-justify by shifting 2*(12 - lw).
            let shift     = (2 * (12 - lut_word_length)) as u32;
            let byte_step = (scan_step / COMPRESSION_RATIO) as usize;
            let mut s_byte = (scan_range[0] / COMPRESSION_RATIO) as usize;
            let s_end      = (scan_range[1] / COMPRESSION_RATIO) as usize;

            while s_byte <= s_end {
                let raw   = ((subject[s_byte]     as u32) << 16)
                          | ((subject[s_byte + 1] as u32) << 8)
                          |  (subject[s_byte + 2] as u32);
                let index = (raw >> shift) as i32 & mask;

                let num_hits = s_blast_lookup_get_num_hits(lookup, index);
                if num_hits == 0 {
                    s_byte += byte_step;
                    continue;
                }
                if num_hits > max_hits - total {
                    break;
                }

                s_blast_lookup_retrieve(
                    lookup,
                    index,
                    &mut offset_pairs[total as usize..],
                    (s_byte as i32) * COMPRESSION_RATIO,
                );
                total  += num_hits;
                s_byte += byte_step;
            }
            scan_range[0] = (s_byte as i32) * COMPRESSION_RATIO;

        } else {
            // Unaligned: use 3 bytes when s%4 + lw ≤ 12, else 4 bytes.
            while scan_range[0] <= scan_range[1] {
                let s_off  = scan_range[0];
                let s_byte = (s_off / COMPRESSION_RATIO) as usize;
                let sum    = s_off % COMPRESSION_RATIO + lut_word_length;

                let index = if sum <= 12 {
                    let raw   = ((subject[s_byte]     as i32) << 16)
                              | ((subject[s_byte + 1] as i32) << 8)
                              |  (subject[s_byte + 2] as i32);
                    let shift = 2 * (12 - sum);
                    (raw >> shift) & mask
                } else {
                    // sum 13–15: word straddles 4th byte.
                    let raw   = ((subject[s_byte]     as u32) << 24)
                              | ((subject[s_byte + 1] as u32) << 16)
                              | ((subject[s_byte + 2] as u32) << 8)
                              |  (subject[s_byte + 3] as u32);
                    let shift = (32 - 2 * sum) as u32;
                    (raw >> shift) as i32 & mask
                };

                let num_hits = s_blast_lookup_get_num_hits(lookup, index);
                if num_hits == 0 {
                    scan_range[0] += scan_step;
                    continue;
                }
                if num_hits > max_hits - total {
                    break;
                }

                s_blast_lookup_retrieve(
                    lookup,
                    index,
                    &mut offset_pairs[total as usize..],
                    s_off,
                );
                total         += num_hits;
                scan_range[0] += scan_step;
            }
        }

    } else if lut_word_length > 5 {
        // lut_word_length 6, 7, or 8 — word spans two packed bytes (or three
        // when the word is not 4-base aligned).  Max sum = 3+8 = 11 ≤ 12,
        // so the 3-byte shift is always non-negative.

        if scan_step % COMPRESSION_RATIO == 0 {
            // Aligned case: every word starts on a 4-base (1-byte) boundary.
            // stride in bytes = scan_step / 4.
            let shift     = 2 * (8 - lut_word_length); // right-justify word
            let byte_step = (scan_step / COMPRESSION_RATIO) as usize;
            let mut s_byte = (scan_range[0] / COMPRESSION_RATIO) as usize;
            let s_end      = (scan_range[1] / COMPRESSION_RATIO) as usize;

            while s_byte <= s_end {
                let raw   = ((subject[s_byte] as i32) << 8)
                          | (subject[s_byte + 1] as i32);
                let index = raw >> shift;

                let num_hits = s_blast_lookup_get_num_hits(lookup, index);
                if num_hits == 0 {
                    s_byte += byte_step;
                    continue;
                }
                if num_hits > max_hits - total {
                    break;
                }

                s_blast_lookup_retrieve(
                    lookup,
                    index,
                    &mut offset_pairs[total as usize..],
                    (s_byte as i32) * COMPRESSION_RATIO,
                );
                total  += num_hits;
                s_byte += byte_step;
            }
            scan_range[0] = (s_byte as i32) * COMPRESSION_RATIO;

        } else {
            // Unaligned case: word may straddle three packed bytes.
            // shift = 2*(12 - (s_off%4 + lut_word_length)) extracts the word.
            while scan_range[0] <= scan_range[1] {
                let s_off  = scan_range[0];
                let shift  = 2 * (12 - (s_off % COMPRESSION_RATIO + lut_word_length));
                let s_byte = (s_off / COMPRESSION_RATIO) as usize;

                let raw   = ((subject[s_byte]     as i32) << 16)
                          | ((subject[s_byte + 1] as i32) << 8)
                          |  (subject[s_byte + 2] as i32);
                let index = (raw >> shift) & mask;

                let num_hits = s_blast_lookup_get_num_hits(lookup, index);
                if num_hits == 0 {
                    scan_range[0] += scan_step;
                    continue;
                }
                if num_hits > max_hits - total {
                    break;
                }

                s_blast_lookup_retrieve(
                    lookup,
                    index,
                    &mut offset_pairs[total as usize..],
                    s_off,
                );
                total         += num_hits;
                scan_range[0] += scan_step;
            }
        }

    } else {
        // lut_word_length 4 or 5 — word fits in two packed bytes.
        // Stride is never a multiple of 4 for these widths.
        while scan_range[0] <= scan_range[1] {
            let s_off  = scan_range[0];
            let shift  = 2 * (8 - (s_off % COMPRESSION_RATIO + lut_word_length));
            let s_byte = (s_off / COMPRESSION_RATIO) as usize;

            let raw   = ((subject[s_byte] as i32) << 8)
                      |  (subject[s_byte + 1] as i32);
            let index = (raw >> shift) & mask;

            let num_hits = s_blast_lookup_get_num_hits(lookup, index);
            if num_hits == 0 {
                scan_range[0] += scan_step;
                continue;
            }
            if num_hits > max_hits - total {
                break;
            }

            s_blast_lookup_retrieve(
                lookup,
                index,
                &mut offset_pairs[total as usize..],
                s_off,
            );
            total         += num_hits;
            scan_range[0] += scan_step;
        }
    }

    total
}

// ── s_NaChooseScanSubject (blast_nascan.c lines 284–293) ─────────────────────

/// Dispatch to the appropriate scanner based on lookup parameters.
///
/// Mirrors `s_NaChooseScanSubject`: selects `s_BlastNaScanSubject_8_4` when
/// lut_word_length == 8 && scan_step == 4; otherwise `s_BlastNaScanSubject_Any`.
pub fn blast_na_scan_subject(
    lookup:       &BlastNaLookupTable,
    subject:      &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    max_hits:     i32,
    scan_range:   &mut [i32; 2],
) -> i32 {
    if lookup.lut_word_length == 8 && lookup.scan_step == 4 {
        s_blast_na_scan_subject_8_4(lookup, subject, offset_pairs, max_hits, scan_range)
    } else {
        s_blast_na_scan_subject_any(lookup, subject, offset_pairs, max_hits, scan_range)
    }
}
