// blast/mblookup.rs — faithful port of NCBI MegaBlast lookup table.
//
// Primary sources (NCBI BLAST 2.17.0+):
//   blast_nalookup.h → BlastMBLookupTable struct
//   blast_nalookup.c → s_FillContigMBTable, BlastMBLookupTableNew
//   blast_nascan.c   → s_BlastMBLookupHasHits, s_BlastMBLookupRetrieve,
//                      s_MBScanSubject_Any
//   blast_lookup.h   → PV_SET, PV_TEST (dynamic pv_array_bts variant)
//
// Used when lut_word_length > 8 (eMBLookupTable in NCBI nomenclature).
// Contiguous megablast only (no discontiguous templates).

use super::types::{BlastOffsetPair, COMPRESSION_RATIO};

// ── Constants (mirrors blast_nalookup.c / blast_lookup.h) ────────────────────

const BITS_PER_NUC:      i32 = 2;
const BLAST2NA_MASK:     u8  = 0xfc; // invalid-base detector (all bits except lower 2)
const PV_ARRAY_BYTES:    i32 = 4;    // sizeof(PV_ARRAY_TYPE) == sizeof(u32)
const PV_ARRAY_BTS:      i32 = 5;    // bits-to-shift for fixed-width PV (BlastNaLookupTable)
const K_TARGET_PV_SIZE:  i64 = 131_072; // 128 KiB target for pv_array
const K_SMALL_QUERY_CUTOFF: usize = 15_000;
const K_LARGE_QUERY_CUTOFF: usize = 800_000;

// ── PV helpers with dynamic pv_array_bts ─────────────────────────────────────

/// PV_SET(pv, index, pv_array_bts) from blast_lookup.h line 49.
/// Sets the bit for `index` in the presence vector.
#[inline]
pub fn mb_pv_set(pv: &mut [u32], index: i64, pv_array_bts: u32) {
    pv[(index >> pv_array_bts) as usize] |= 1u32 << (index & 31);
}

/// PV_TEST(pv, index, pv_array_bts) from blast_lookup.h line 55.
#[inline]
pub fn mb_pv_test(pv: &[u32], index: i64, pv_array_bts: u32) -> bool {
    pv[(index >> pv_array_bts) as usize] & (1u32 << (index & 31)) != 0
}

// ── ilog2 (lookup_util.h) ─────────────────────────────────────────────────────

/// Integer floor-log2 for positive i64.  Equivalent to NCBI's `ilog2`.
#[inline]
fn ilog2(x: i64) -> u32 {
    debug_assert!(x > 0);
    63 - (x as u64).leading_zeros()
}

// ── BlastMBLookupTable (blast_nalookup.h lines 189–232) ──────────────────────

/// Contiguous MegaBlast lookup table.
///
/// hashtable and next_pos implement a per-ecode singly-linked list of 1-based
/// query offsets: `hashtable[ecode]` is the most-recently added 1-based offset;
/// `next_pos[q_off]` is the previously added 1-based offset for the same ecode,
/// or 0 (end of chain).  Retrieve converts to 0-based by subtracting 1.
pub struct BlastMBLookupTable {
    pub word_length:     i32,
    pub lut_word_length: i32,
    pub hashsize:        i64,  // = 4^lut_word_length
    pub scan_step:       i32,  // = word_length - lut_word_length + 1
    pub hashtable:       Vec<i32>, // size = hashsize; 0 = empty, else 1-based q_off
    pub next_pos:        Vec<i32>, // size = query_length + 1; next_pos[q1] = prev q1
    pub pv_array:        Vec<u32>,
    pub pv_array_bts:    u32,
    pub longest_chain:   i32,
}

// ── s_BlastMBLookupHasHits / s_BlastMBLookupRetrieve ─────────────────────────

/// True if the presence vector records any hits for `index`.
#[inline]
pub fn mb_has_hits(lookup: &BlastMBLookupTable, index: i64) -> bool {
    mb_pv_test(&lookup.pv_array, index, lookup.pv_array_bts)
}

/// Copy all query offsets for `index` into `offset_pairs`.
/// Returns the number of pairs written.
#[inline]
pub fn mb_retrieve(
    lookup:       &BlastMBLookupTable,
    index:        i64,
    offset_pairs: &mut [BlastOffsetPair],
    s_off:        i32,
) -> i32 {
    let mut i   = 0usize;
    let mut q_off = lookup.hashtable[index as usize];
    while q_off != 0 {
        offset_pairs[i].q_off = (q_off - 1) as u32; // 1-based → 0-based
        offset_pairs[i].s_off = s_off as u32;
        i += 1;
        q_off = lookup.next_pos[q_off as usize];
    }
    i as i32
}

// ── s_FillContigMBTable (blast_nalookup.c lines 949–1111) ────────────────────

/// Index query positions into the MB lookup table.
///
/// `query_sequence` — 0-based BLASTNA (no sentinel).
/// `locations`      — `(left, right)` inclusive 0-based base ranges.
///
/// Mirrors s_FillContigMBTable (contiguous megablast; no discontiguous templates,
/// no db_filter).  All stored indices are 1-based query offsets.
fn s_fill_contig_mb_table(
    query_sequence: &[u8],
    locations:      &[(i32, i32)],
    lookup:         &mut BlastMBLookupTable,
) {
    let lw     = lookup.lut_word_length as usize;
    let wl     = lookup.word_length     as usize;
    let mask   = lookup.hashsize - 1;

    // NCBI uses a compressed helper to estimate longest_chain.
    // helper_array[ecode / K_COMPRESSION] counts the number of "collision"
    // insertions (second-and-beyond hits) in each group.  The max across all
    // groups is a conservative overestimate of the longest single chain.
    // (Matches s_FillContigMBTable in blast_nalookup.c lines 971–973.)
    const K_COMPRESSION: usize = 2048;
    let helper_len = ((lookup.hashsize as usize) / K_COMPRESSION).max(1);
    let mut helper = vec![0u32; helper_len];

    for &(left, right) in locations {
        let left  = left  as usize;
        let right = right as usize;

        // Skip regions too short to contain even one full word.
        if wl > right - left + 1 {
            continue;
        }

        let mut ecode: i64 = 0;
        // Minimum read_pos at which we have accumulated lw valid bases.
        // After ambiguity at position P, valid_from = P + lw (next complete window).
        let mut valid_from: usize = left + lw - 1;

        for read_pos in left..=right {
            let base = query_sequence[read_pos];

            if base & BLAST2NA_MASK != 0 {
                // Ambiguous base: reset rolling hash; next valid window starts lw
                // bases later (same logic as `pos = seq + kLutWordLength` in C).
                ecode      = 0;
                valid_from = read_pos + lw;
                continue;
            }

            ecode = ((ecode << BITS_PER_NUC) & mask) | (base as i64);

            if read_pos < valid_from {
                continue;
            }

            // Complete, unambiguous lw-mer.
            // q_off_1based = read_pos - lw + 2 == 0-based_q_off + 1.
            let q_off_1 = (read_pos - lw + 2) as i32;

            if lookup.hashtable[ecode as usize] == 0 {
                mb_pv_set(&mut lookup.pv_array, ecode, lookup.pv_array_bts);
            } else {
                // Collision: a second (or later) hit for this k-mer.
                // Increment the compressed helper for longest_chain estimation.
                helper[(ecode as usize) / K_COMPRESSION] += 1;
            }
            lookup.next_pos[q_off_1 as usize] = lookup.hashtable[ecode as usize];
            lookup.hashtable[ecode as usize]   = q_off_1;
        }
    }

    let longest = helper.iter().copied().max().unwrap_or(0);
    lookup.longest_chain = (longest as i32).max(2);
}

// ── BlastMBLookupTableNew (blast_nalookup.c lines 1227–1360) ─────────────────

/// Build a BlastMBLookupTable from the query.
///
/// `query_sequence`  — 0-based BLASTNA (no sentinel).
/// `locations`       — `(left, right)` inclusive base ranges to index.
/// `word_length`     — full word length (e.g. 11).
/// `lut_width`       — lookup table word length (9–12); **must be > 8**.
/// `approx_entries`  — approximate number of query positions to be indexed
///                     (used to scale pv_array; typically 2 × query_len for
///                     two-strand search, or 1 × query_len for single strand).
pub fn blast_mb_lookup_table_new(
    query_sequence: &[u8],
    locations:      &[(i32, i32)],
    word_length:    i32,
    lut_width:      i32,
    approx_entries: usize,
) -> BlastMBLookupTable {
    debug_assert!(lut_width >= 9 && lut_width <= 12,
        "BlastMBLookupTable requires lut_width 9–12, got {}", lut_width);

    let hashsize = 1i64 << (BITS_PER_NUC * lut_width);

    // ── PV array sizing (mirrors BlastMBLookupTableNew lines 1286–1305) ──────
    // For lut_word_length ≤ 12:
    //   if hashsize ≤ 8 × kTargetPVSize: pv_size = hashsize >> PV_ARRAY_BTS
    //   else:                             pv_size = kTargetPVSize / PV_ARRAY_BYTES
    // Then halve if query is very small (≤15000) or very large (≥800000).
    let mut pv_size = if hashsize <= 8 * K_TARGET_PV_SIZE {
        (hashsize >> PV_ARRAY_BTS) as i64
    } else {
        K_TARGET_PV_SIZE / PV_ARRAY_BYTES as i64
    };
    if approx_entries <= K_SMALL_QUERY_CUTOFF || approx_entries >= K_LARGE_QUERY_CUTOFF {
        pv_size /= 2;
    }
    let pv_array_bts = ilog2(hashsize / pv_size);
    let pv_len       = pv_size as usize; // each element is one u32

    let query_length = query_sequence.len();

    let mut lookup = BlastMBLookupTable {
        word_length:     word_length,
        lut_word_length: lut_width,
        hashsize,
        scan_step:       word_length - lut_width + 1,
        hashtable:       vec![0i32; hashsize as usize],
        next_pos:        vec![0i32; query_length + 1],
        pv_array:        vec![0u32; pv_len],
        pv_array_bts,
        longest_chain:   0,
    };

    s_fill_contig_mb_table(query_sequence, locations, &mut lookup);

    lookup
}

// ── s_MBScanSubject_Any (blast_nascan.c lines 1482–1626) ─────────────────────

/// Scan the compressed subject sequence for MB lookup hits.
///
/// Handles lut_word_length 9–12 with any scan_step (aligned or unaligned).
///
/// `subject`      — NCBI2NA packed bytes (base 0 in byte 0, MSB-first).
/// `scan_range`   — `[start_base, end_base]` inclusive; on return `[0]` holds
///                  the next unscanned base position.
/// Returns total hits written into `offset_pairs[0..return_value]`.
pub fn blast_mb_scan_subject(
    lookup:       &BlastMBLookupTable,
    subject:      &[u8],
    offset_pairs: &mut [BlastOffsetPair],
    max_hits:     i32,
    scan_range:   &mut [i32; 2],
) -> i32 {
    let lut_word_length = lookup.lut_word_length;
    let scan_step       = lookup.scan_step;
    let mask            = lookup.hashsize - 1;
    let mut total       = 0i32;

    // Pre-subtract longest_chain: NCBI checks total_hits >= max_hits BEFORE
    // calling retrieve, where max_hits has already been reduced by longest_chain.
    let max_hits = max_hits - lookup.longest_chain;

    if scan_step % COMPRESSION_RATIO == 0 {
        // ── Aligned: every word starts on a 4-base boundary ─────────────────
        // 3 packed bytes = 12 bases; right-justify by 2*(12 - lut_word_length).
        let shift     = (2 * (12 - lut_word_length)) as u32;
        let byte_step = (scan_step / COMPRESSION_RATIO) as usize;
        let mut s_byte = (scan_range[0] / COMPRESSION_RATIO) as usize;
        let s_end      = (scan_range[1] / COMPRESSION_RATIO) as usize;

        while s_byte <= s_end {
            let raw   = ((subject[s_byte]     as u32) << 16)
                      | ((subject[s_byte + 1] as u32) << 8)
                      |  (subject[s_byte + 2] as u32);
            let index = (raw >> shift) as i64; // no &mask needed (already fits)

            if mb_has_hits(lookup, index) {
                if total >= max_hits { break; }
                total += mb_retrieve(
                    lookup, index,
                    &mut offset_pairs[total as usize..],
                    (s_byte as i32) * COMPRESSION_RATIO,
                );
            }
            s_byte += byte_step;
        }
        scan_range[0] = (s_byte as i32) * COMPRESSION_RATIO;

    } else if lut_word_length > 9 {
        // ── Unaligned, lut_word_length 10–12: 4-byte read ───────────────────
        // shift = 2*(16 - (s_off%4 + lw)); always non-negative since s_off%4 ≤ 3
        // and lw ≤ 12, so sum ≤ 15 < 16.
        while scan_range[0] <= scan_range[1] {
            let s_off  = scan_range[0];
            let s_byte = (s_off / COMPRESSION_RATIO) as usize;
            let sum    = s_off % COMPRESSION_RATIO + lut_word_length;
            let shift  = (2 * (16 - sum)) as u32;

            let raw   = ((subject[s_byte]     as u32) << 24)
                      | ((subject[s_byte + 1] as u32) << 16)
                      | ((subject[s_byte + 2] as u32) << 8)
                      |  (subject[s_byte + 3] as u32);
            let index = ((raw >> shift) as i64) & mask;

            if mb_has_hits(lookup, index) {
                if total >= max_hits { break; }
                total += mb_retrieve(
                    lookup, index,
                    &mut offset_pairs[total as usize..],
                    s_off,
                );
            }
            scan_range[0] += scan_step;
        }

    } else {
        // ── Unaligned, lut_word_length == 9: 3-byte read ─────────────────────
        // shift = 2*(12 - (s_off%4 + 9)); sum ∈ [9,12], shift ∈ [0,6], always ≥ 0.
        while scan_range[0] <= scan_range[1] {
            let s_off  = scan_range[0];
            let s_byte = (s_off / COMPRESSION_RATIO) as usize;
            let sum    = s_off % COMPRESSION_RATIO + lut_word_length;
            let shift  = (2 * (12 - sum)) as u32;

            let raw   = ((subject[s_byte]     as u32) << 16)
                      | ((subject[s_byte + 1] as u32) << 8)
                      |  (subject[s_byte + 2] as u32);
            let index = ((raw >> shift) as i64) & mask;

            if mb_has_hits(lookup, index) {
                if total >= max_hits { break; }
                total += mb_retrieve(
                    lookup, index,
                    &mut offset_pairs[total as usize..],
                    s_off,
                );
            }
            scan_range[0] += scan_step;
        }
    }

    total
}
