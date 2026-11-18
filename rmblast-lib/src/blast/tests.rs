//! Tests of the lookup tables, the subject scanners and the packing
//! utilities against brute-force oracles.

use super::mblookup::{blast_mb_lookup_table_new, blast_mb_scan_subject};
use super::nalookup::{blast_na_lookup_table_new, choose_lut_width};
use super::nascan::blast_na_scan_subject;
use super::types::BlastOffsetPair;
use super::util::{blast_compress_blastna_sequence, ncbi2na_unpack_base, seqblk_from_blastna};
use std::collections::BTreeSet;

/// xorshift64: deterministic, no dependency.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next() % n as u64) as usize
    }
}

/// A random ACGT query, and a subject built from pieces of it (some
/// mutated) and random filler, so that word matches are common.
fn query_and_subject(rng: &mut Rng, q_len: usize, s_len: usize) -> (Vec<u8>, Vec<u8>) {
    let q: Vec<u8> = (0..q_len).map(|_| rng.below(4) as u8).collect();
    let mut s = Vec::with_capacity(s_len);
    while s.len() < s_len {
        if rng.below(3) == 0 {
            let n = 1 + rng.below(40);
            s.extend((0..n).map(|_| rng.below(4) as u8));
        } else {
            let start = rng.below(q_len);
            let n = (1 + rng.below(60)).min(q_len - start);
            for &b in &q[start..start + n] {
                s.push(if rng.below(25) == 0 { rng.below(4) as u8 } else { b });
            }
        }
    }
    s.truncate(s_len);
    (q, s)
}

/// The packed subject as the scanners read it (packed[3..]).
fn packed(subject: &[u8]) -> Vec<u8> {
    let mut blk = seqblk_from_blastna(subject);
    blast_compress_blastna_sequence(&mut blk);
    blk.packed[3..].to_vec()
}

/// Drive a scanner as the engine does (scan_range [0, len - lw], repeated
/// calls until it passes the end) and collect every (q_off, s_off).
fn scan_all(
    s_len: usize,
    lw: usize,
    longest_chain: i32,
    mut scan: impl FnMut(&mut [BlastOffsetPair], i32, &mut [i32; 2]) -> i32,
) -> Vec<(u32, u32)> {
    let end_range = (s_len - lw) as i32;
    let batch_len = 4096usize.max(longest_chain as usize + 4096);
    let mut batch = vec![BlastOffsetPair::default(); batch_len];
    let mut scan_range = [0i32, end_range];
    let mut out = Vec::new();
    loop {
        let n = scan(&mut batch, batch_len as i32, &mut scan_range);
        out.extend(batch[..n as usize].iter().map(|p| (p.q_off, p.s_off)));
        if scan_range[0] > end_range {
            break;
        }
    }
    out
}

/// Every lut-word match at the subject offsets a scanner with this stride
/// visits (0, step, 2 step, ... up to len - lw).
fn oracle(q: &[u8], s: &[u8], lw: usize, step: usize) -> BTreeSet<(u32, u32)> {
    let mut want = BTreeSet::new();
    for so in (0..=s.len() - lw).step_by(step) {
        for qo in 0..=q.len() - lw {
            if q[qo..qo + lw] == s[so..so + lw] {
                want.insert((qo as u32, so as u32));
            }
        }
    }
    want
}

/// The small (width <= 8) and megablast (9..12) tables with their scanners
/// report exactly the lut-word matches at the visited subject offsets: no
/// false hits, none missed, none twice.  Every width each word size can use,
/// so the aligned 8/4 scanner, the generic scanner's aligned and unaligned
/// strides, and the megablast scanner are all covered.
#[test]
fn scanners_report_exactly_the_word_matches() {
    let mut rng = Rng(0x2545F4914F6CDD1D);
    for word_size in [4usize, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16, 20, 28] {
        for lw in 4..=word_size.min(12) {
            for trial in 0..3 {
                let (q, s) = query_and_subject(&mut rng, 150 + 50 * trial, 1500 + 7 * trial);
                let locs = [(0i32, q.len() as i32 - 1)];
                let sp = packed(&s);
                let (got, step) = if lw <= 8 {
                    let t = blast_na_lookup_table_new(&q, &locs, word_size as i32, lw as i32);
                    let step = t.scan_step as usize;
                    (scan_all(s.len(), lw, t.longest_chain,
                        |b, m, r| blast_na_scan_subject(&t, &sp, b, m, r)), step)
                } else {
                    let t = blast_mb_lookup_table_new(&q, &locs, word_size as i32, lw as i32, q.len());
                    let step = t.scan_step as usize;
                    (scan_all(s.len(), lw, t.longest_chain,
                        |b, m, r| blast_mb_scan_subject(&t, &sp, b, m, r)), step)
                };
                assert_eq!(step, word_size - lw + 1, "ws {word_size} lw {lw}: scan step");
                let got_set: BTreeSet<(u32, u32)> = got.iter().copied().collect();
                assert_eq!(got_set.len(), got.len(), "ws {word_size} lw {lw}: a hit reported twice");
                let want = oracle(&q, &s, lw, step);
                assert!(!want.is_empty(), "ws {word_size} lw {lw}: test data has no matches");
                assert_eq!(got_set, want, "ws {word_size} lw {lw} trial {trial}");
            }
        }
    }
}

/// An ambiguous base (N) in the query keeps every word that contains it out
/// of the table.
#[test]
fn query_words_with_ambiguous_bases_are_not_indexed() {
    let mut rng = Rng(0x9E3779B97F4A7C15);
    let (mut q, s) = query_and_subject(&mut rng, 200, 1500);
    for i in (10..200).step_by(23) {
        q[i] = 14; // N
    }
    let locs = [(0i32, q.len() as i32 - 1)];
    let lw = 8;
    let t = blast_na_lookup_table_new(&q, &locs, 11, lw as i32);
    let sp = packed(&s);
    let got: BTreeSet<(u32, u32)> = scan_all(s.len(), lw, t.longest_chain,
        |b, m, r| blast_na_scan_subject(&t, &sp, b, m, r)).into_iter().collect();
    let want: BTreeSet<(u32, u32)> = oracle(&q, &s, lw, t.scan_step as usize)
        .into_iter()
        .filter(|&(qo, _)| q[qo as usize..qo as usize + lw].iter().all(|&b| b < 4))
        .collect();
    assert!(!want.is_empty());
    assert_eq!(got, want);
}

/// choose_lut_width matches NCBI's BlastChooseNaLookupTable thresholds at
/// their boundaries (blast_nalookup.c), and never exceeds the word size.
#[test]
fn lut_width_thresholds_like_ncbi() {
    let cases: &[(usize, usize, usize)] = &[
        (4, 1, 4), (6, 1_000_000, 6),
        (7, 249, 6), (7, 250, 7),
        (8, 8_499, 7), (8, 8_500, 8),
        (9, 1_249, 7), (9, 1_250, 8), (9, 20_999, 8), (9, 21_000, 9),
        (10, 1_249, 7), (10, 1_250, 8), (10, 8_499, 8), (10, 8_500, 9), (10, 17_999, 9), (10, 18_000, 10),
        (11, 11_999, 8), (11, 12_000, 10), (11, 179_999, 10), (11, 180_000, 11),
        (12, 8_499, 8), (12, 8_500, 9), (12, 17_999, 9), (12, 18_000, 10),
        (12, 59_999, 10), (12, 60_000, 11), (12, 899_999, 11), (12, 900_000, 12),
        (14, 8_499, 8), (14, 8_500, 11), (14, 299_999, 11), (14, 300_000, 12),
    ];
    for &(ws, entries, want) in cases {
        assert_eq!(choose_lut_width(ws, entries), want, "word size {ws}, {entries} entries");
    }
    for ws in 4..=40 {
        for entries in [1, 1_000, 100_000, 10_000_000] {
            assert!(choose_lut_width(ws, entries) <= ws.min(12));
        }
    }
}

/// Packing round-trips every base for every length mod 4, with the layout
/// the scanners assume: 3 zero bytes, then 4 bases per byte MSB first, a
/// left-justified last byte, then zero padding.
#[test]
fn packing_round_trips_for_every_length() {
    let mut rng = Rng(0xD1B54A32D192ED03);
    for len in 1..=64usize {
        let s: Vec<u8> = (0..len).map(|_| rng.below(4) as u8).collect();
        let mut blk = seqblk_from_blastna(&s);
        blast_compress_blastna_sequence(&mut blk);
        assert_eq!(&blk.packed[..3], &[0, 0, 0]);
        let data = &blk.packed[3..];
        for (i, &b) in s.iter().enumerate() {
            assert_eq!(ncbi2na_unpack_base(data[i / 4], 3 - (i % 4) as u32), b, "len {len} base {i}");
        }
        if len % 4 != 0 {
            let used = 2 * (len % 4) as u32;
            assert_eq!(data[len / 4] & (0xFFu8 >> used), 0, "len {len}: last byte not left-justified");
        }
        assert!(data[len.div_ceil(4)..].iter().all(|&b| b == 0), "len {len}: tail pad not zero");
    }
}

/// An ambiguous base packs to one of the bases it stands for, the draw
/// depends only on the sequence (seeded by its length, as makeblastdb does),
/// and unambiguous bases are untouched.
#[test]
fn packing_ambiguous_bases_draws_within_their_sets() {
    // BLASTNA codes 4..=14 and the bases each allows (encoding::BLASTNA_TO_IUPAC order).
    let sets: [&[u8]; 11] = [&[0, 2], &[1, 3], &[0, 1], &[2, 3], &[0, 3], &[1, 2],
                             &[1, 2, 3], &[0, 2, 3], &[0, 1, 3], &[0, 1, 2], &[0, 1, 2, 3]];
    let mut rng = Rng(0x94D049BB133111EB);
    let s: Vec<u8> = (0..400).map(|_| if rng.below(3) == 0 { 4 + rng.below(11) as u8 } else { rng.below(4) as u8 }).collect();
    let pack = |s: &[u8]| { let mut b = seqblk_from_blastna(s); blast_compress_blastna_sequence(&mut b); b.packed };
    let p = pack(&s);
    assert_eq!(p, pack(&s), "same sequence, different draws");
    for (i, &code) in s.iter().enumerate() {
        let got = ncbi2na_unpack_base(p[3 + i / 4], 3 - (i % 4) as u32);
        if code < 4 {
            assert_eq!(got, code, "base {i}");
        } else {
            assert!(sets[code as usize - 4].contains(&got), "base {i}: code {code} packed as {got}");
        }
    }
}
