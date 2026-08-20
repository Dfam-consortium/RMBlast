//! Symmetric DUST masking — faithful port of NCBI CSymDustMasker (symdust.cpp).
//!
//! Default parameters: window=64, level=20, linker=1.
//! Input: BLASTNA-encoded sequence with leading/trailing sentinels (seq[0], seq[n+1]).
//! Only the low 2 bits of each base are used (A=0 C=1 G=2 T=3).

use std::collections::VecDeque;

pub const DUST_WINDOW: usize = 64;
pub const DUST_LEVEL:  u32   = 20;
pub const DUST_LINKER: usize = 1;

// ── NCBI CRandom LFG ─────────────────────────────────────────────────────────

/// Port of NCBI's CRandom Lagged Fibonacci Generator (eGetRand_LFG).
/// Uses the same hard-coded initial state as NCBI so that N-base randomisation
/// produces byte-identical DUST results.
///
/// Public because a single generator must be threaded across several
/// `dust_intervals` calls on one record: NCBI's `CSymDustMasker` owns the
/// converter (and hence the RNG) as a member, so the N-split segments of
/// `GetDustMasks_SkipNs` share one RNG stream.  See [`dustmasker_intervals`].
pub struct NcbiLfg {
    state: [u32; 33],
    rk: isize,  // signed so we can detect < 0 wraparound
    rj: isize,
}

impl Default for NcbiLfg {
    fn default() -> Self { Self::new() }
}

impl NcbiLfg {
    pub fn new() -> Self {
        NcbiLfg {
            state: [
                0xd53f1852, 0xdfc78b83, 0x4f256096, 0x0e643df7,
                0x82c359bf, 0xc7794dfa, 0xd5e9ffaa, 0x2c8cb64a,
                0x2f07b334, 0xad5a7eb5, 0x96dc0cde, 0x6fc24589,
                0xa5853646, 0xe71576e2, 0x0dae30df, 0xb09ce711,
                0x5e56ef87, 0x4b4b0082, 0x6f4f340e, 0xc5bb17e8,
                0xd788d765, 0x67498087, 0x9d7aba26, 0x261351d4,
                0x411ee7ea, 0x0393a263, 0x2c5a5835, 0xc115fcd8,
                0x25e9132c, 0xd0c6e906, 0xc2bc5b2d, 0x6c065c98,
                0x6e37bd55,
            ],
            rk: 32, // kStateSize - 1
            rj: 12, // kStateOffset
        }
    }

    /// Mirror of NCBI CRandom::x_GetRand32Bits() >> 1.
    #[inline]
    fn get_rand(&mut self) -> u32 {
        let r = self.state[self.rk as usize].wrapping_add(self.state[self.rj as usize]);
        self.state[self.rk as usize] = r;
        self.rk -= 1;
        self.rj -= 1;
        if self.rk < 0 { self.rk = 32; }
        else if self.rj < 0 { self.rj = 32; }
        r >> 1
    }
}

// ── triplet encoding ──────────────────────────────────────────────────────────

/// Convert one BLASTNA base to ncbi2na for DUST, mirroring NCBI's
/// `CIupac2Ncbi2na_converter`: C→1, G→2, T→3, N→`CRandom.GetRand()&3`, else →0.
/// Advances `rng` only on N (code 14), matching NCBI's per-N RNG consumption.
#[inline]
fn cvt_base(b: u8, rng: &mut NcbiLfg) -> u8 {
    match b {
        1 => 1, // C
        2 => 2, // G
        3 => 3, // T
        14 => (rng.get_rand() & 3) as u8, // N
        _ => 0, // A and all other ambiguity codes (matches NCBI converter default)
    }
}

// ── perfect interval ─────────────────────────────────────────────────────────

/// Mirror of NCBI's `struct perfect { bounds_, score_, len_ }`.
/// All positions are absolute (0-indexed in the sequence).
#[derive(Clone)]
struct Perfect {
    start: usize,  // bounds_.first  = pos
    stop:  usize,  // bounds_.second = stop_ + 1
    score: u32,    // score_
    len:   u32,    // len_  (= count at insertion time)
}

// ── triplets window state ─────────────────────────────────────────────────────

/// Running window state — mirrors NCBI's `class triplets`.
/// All positions (`start`, `stop`, `l`) are absolute indices into `bases[]`.
struct Triplets {
    list:     VecDeque<u8>,  // front=newest (stop), back=oldest (start)
    c_w:      [u32; 64],
    r_w:      u32,
    c_v:      [u32; 64],
    r_v:      u32,
    num_diff: u32,
    max_size: usize,   // window − 2
    low_k:    u32,     // level / 5
    pub start:    usize,   // absolute position of window left
    pub stop:     usize,   // absolute position of window right (just added)
    pub l:        usize,   // absolute position of suffix left (L in NCBI)
}

impl Triplets {
    fn new(window: usize, low_k: u32, j_start: usize) -> Self {
        Triplets {
            list: VecDeque::new(),
            c_w: [0; 64], r_w: 0,
            c_v: [0; 64], r_v: 0,
            num_diff: 0,
            max_size: window.saturating_sub(2),
            low_k,
            start: j_start, stop: j_start, l: j_start,
        }
    }

    #[inline] fn add_w(&mut self, t: u8) {
        let i = t as usize;
        if self.c_w[i] == 0 { self.num_diff += 1; }
        self.r_w += self.c_w[i];
        self.c_w[i] += 1;
    }
    #[inline] fn rem_w(&mut self, t: u8) {
        let i = t as usize;
        self.c_w[i] -= 1;
        if self.c_w[i] == 0 { self.num_diff -= 1; }
        self.r_w -= self.c_w[i];
    }
    #[inline] fn add_v(&mut self, t: u8) {
        let i = t as usize;
        self.r_v += self.c_v[i];
        self.c_v[i] += 1;
    }
    #[inline] fn rem_v(&mut self, t: u8) {
        let i = t as usize;
        self.c_v[i] -= 1;
        self.r_v -= self.c_v[i];
    }

    /// Mirror of NCBI `shift_high`.
    fn shift_high(&mut self, t: u8, perfect: &mut Vec<Perfect>) -> bool {
        let s = *self.list.back().unwrap();
        self.list.pop_back();
        self.rem_w(s);
        self.start += 1;

        self.list.push_front(t);
        self.add_w(t);
        self.stop += 1;

        if self.num_diff <= 1 {
            perfect.insert(0, Perfect { start: self.start, stop: self.stop + 1, score: 0, len: 0 });
            false
        } else {
            true
        }
    }

    /// Mirror of NCBI `shift_window`. Returns false when window is trivial (all one triplet).
    pub fn shift_window(&mut self, t: u8, perfect: &mut Vec<Perfect>) -> bool {
        if self.list.len() >= self.max_size {
            if self.num_diff <= 1 {
                return self.shift_high(t, perfect);
            }
            let s = *self.list.back().unwrap();
            self.list.pop_back();
            self.rem_w(s);
            if self.l == self.start {
                self.l += 1;
                self.rem_v(s);
            }
            self.start += 1;
        }

        self.list.push_front(t);
        self.add_w(t);
        self.add_v(t);

        // If suffix contains too many copies of t, advance l past the oldest one.
        if self.c_v[t as usize] > self.low_k {
            // off = index of oldest suffix element in deque
            let off = self.list.len() - (self.l - self.start) - 1;
            let mut idx = off;
            loop {
                let s = self.list[idx];
                self.rem_v(s);
                self.l += 1;
                if s == t { break; }
                if idx == 0 { break; }
                idx -= 1;
            }
        }

        self.stop += 1;

        if self.list.len() >= self.max_size && self.num_diff <= 1 {
            perfect.clear();
            perfect.insert(0, Perfect { start: self.start, stop: self.stop + 1, score: 0, len: 0 });
            false
        } else {
            true
        }
    }

    /// Mirror of NCBI `needs_processing` (Proposition 2).
    pub fn needs_processing(&self, thresholds: &[u32]) -> bool {
        let count = self.stop - self.l;          // stop_ − L
        count < self.list.len() && self.r_w * 10 > thresholds[count]
    }

    /// Mirror of NCBI `find_perfect`.
    /// `perfect` is ordered with index-0 = largest `start`, last = smallest `start`.
    pub fn find_perfect(&self, perfect: &mut Vec<Perfect>, thresholds: &[u32]) {
        let mut counts = [0u32; 64];
        counts.copy_from_slice(&self.c_v);
        let mut score: u32 = self.r_v;

        // count = stop_ − L (number of suffix elements to skip, matching NCBI init)
        let mut count = (self.stop - self.l) as u32;

        // perfect_iter starts at the front (index 0, largest start)
        let mut perfect_iter: usize = 0;
        let mut max_perfect_score: u32 = 0;
        let mut max_len: u32 = 0;

        // deque[stop-l .. list.len()-1] are the pre-suffix triplets (oldest = back)
        // deque[k] = absolute position (stop - k)
        let start_idx = (self.stop - self.l) as usize; // = count as usize

        for idx in start_idx..self.list.len() {
            let ti = self.list[idx];
            let cnt = counts[ti as usize];
            score += counts[ti as usize];
            counts[ti as usize] += 1;

            if cnt > 0 && score * 10 > thresholds[count as usize] {
                // Absolute sequence position of the triplet at deque[idx] = stop - idx.
                // Left boundary of new interval = one before that triplet = stop - idx - 1.
                let seq_pos = self.stop - idx;
                let pos = seq_pos.wrapping_sub(1); // matches NCBI size_type underflow when L=0

                // Advance perfect_iter past all existing intervals with start >= pos.
                while perfect_iter < perfect.len() && pos <= perfect[perfect_iter].start {
                    let pi = &perfect[perfect_iter];
                    if max_perfect_score == 0
                        || max_len * pi.score > max_perfect_score * pi.len
                    {
                        max_perfect_score = pi.score;
                        max_len = pi.len;
                    }
                    perfect_iter += 1;
                }

                if max_perfect_score == 0 || score * max_len >= max_perfect_score * count {
                    max_perfect_score = score;
                    max_len = count;
                    perfect.insert(perfect_iter, Perfect {
                        start: pos,
                        stop:  self.stop + 1,
                        score: max_perfect_score,
                        len:   count,
                    });
                    // perfect_iter now points to the newly inserted element; don't advance.
                }
            }
            count += 1;
        }
    }
}

// ── save_masked_regions ───────────────────────────────────────────────────────

/// Mirror of NCBI `save_masked_regions`.
/// Outputs at most ONE interval (P.back() = smallest start) and discards all
/// expired intervals (start < wstart). Merges with linker.
fn save_masked(
    perfect: &mut Vec<Perfect>,
    res: &mut Vec<(usize, usize)>,
    wstart: usize,
    linker: usize,
) {
    if perfect.is_empty() { return; }

    let back_start = perfect.last().unwrap().start;
    let back_stop  = perfect.last().unwrap().stop;

    if back_start < wstart {
        let b0 = back_start;
        let b1 = back_stop;

        let should_merge = res.last().map_or(false, |l| l.1 + linker >= b0);
        if should_merge {
            let last = res.last_mut().unwrap();
            if b1 > last.1 { last.1 = b1; }
        } else {
            res.push((b0, b1));
        }

        while !perfect.is_empty() && perfect.last().unwrap().start < wstart {
            perfect.pop();
        }
    }
}

// ── public API ────────────────────────────────────────────────────────────────

/// Clamp DUST parameters to the ranges enforced by NCBI's `CSymDustMasker`
/// constructor (`symdust.cpp:170`): out-of-range values silently fall back to
/// the defaults rather than being rejected.
pub fn clamp_dust_params(window: usize, level: u32, linker: usize) -> (usize, u32, usize) {
    (
        if (8..=64).contains(&window) { window } else { DUST_WINDOW },
        if (2..=64).contains(&level)  { level }  else { DUST_LEVEL },
        if (1..=32).contains(&linker) { linker } else { DUST_LINKER },
    )
}

/// Compute DUST intervals over a run of BLASTNA bases (NO sentinels).
///
/// Mirrors `CSymDustMasker::operator()(seq, start, stop)` for the whole of
/// `bases`; to dust a sub-range, pass the sub-slice and offset the returned
/// intervals (NCBI's range form is the same code with an added offset — see
/// `save_masked_regions(res, wstart, start)`).
///
/// `rng` is threaded in by the caller because NCBI's masker owns the
/// IUPAC→ncbi2na converter — and therefore the CRandom stream — as a member,
/// so several calls covering one record share one stream.  See
/// [`dustmasker_intervals`].
///
/// Returns merged, inclusive `[start, stop]` intervals relative to `bases[0]`.
pub fn dust_intervals(
    bases: &[u8],
    window: usize,
    level: u32,
    linker: usize,
    rng: &mut NcbiLfg,
) -> Vec<(usize, usize)> {
    let (window, level, linker) = clamp_dust_params(window, level, linker);
    let n = bases.len();
    if n < 3 { return Vec::new(); }

    // Base conversion to ncbi2na for the triplet computation, mirroring NCBI
    // `CSymDustMasker::CIupac2Ncbi2na_converter`: C→1, G→2, T→3, N→`CRandom.GetRand()&3`,
    // everything else (A and all other ambiguity codes) →0.  `NcbiLfg` is a faithful port
    // of NCBI's `CRandom` (same hard-coded state, lags, `>>1`).
    //
    // CRITICAL (bug #31): NCBI calls the converter *inline* while sliding the window, and
    // its outer loop (`while stop>2+start`) re-creates the sequence iterator at `start` on
    // every restart — so bases in the re-visited range are converted AGAIN, advancing the
    // shared CRandom and yielding DIFFERENT random values at N positions on re-visits.
    // Converting the whole sequence once up front (the old approach) desyncs the RNG after
    // the first masked region triggers a restart, producing wrong N values downstream
    // (e.g. HAL1ME poly-T+NN @2422-2428 stayed unmasked → spurious seeds).  We replicate
    // NCBI by keeping `rng` persistent across restarts and re-filling `conv` from `j_start`
    // each outer iteration.
    let mut conv = vec![0u8; n]; // ncbi2na of real bases, (re)filled per outer iteration

    // Threshold table: thresholds[0]=1, thresholds[i]=i*level for i>=1.
    let mut thresholds = vec![1u32; window.saturating_sub(2)];
    for i in 1..thresholds.len() {
        thresholds[i] = (i as u32) * level;
    }

    let low_k = level / 5;
    let mut res: Vec<(usize, usize)> = Vec::new();
    let mut perfect: Vec<Perfect> = Vec::new();

    // NCBI: stop = seq.size() - 1 = n - 1 (last valid base index, 0-indexed).
    // Condition "stop > 2 + start" ↔ n - 1 > 2 + j_start ↔ n >= j_start + 4.
    let mut j_start = 0usize; // "start" in NCBI outer loop

    'outer: loop {
        // Need n - j_start >= 4 (at least one complete triplet with room for rolling window).
        if n < j_start + 4 { break; }

        perfect.clear();
        let mut w = Triplets::new(window, low_k, j_start);

        // Number of triplets in range [j_start, n-3] (inclusive) = n - 2 - j_start.
        let num_local = n - 2 - j_start; // guaranteed >= 1

        // NCBI initialises a partial 2-base triplet for bases [start, start+1],
        // then enters the loop at it.GetPos() = start+2.  In our index scheme,
        // j is the absolute index of the next triplet to process.
        let mut j = j_start;   // absolute triplet index
        let j_end = j_start + num_local; // exclusive end

        // Re-convert bases from j_start this outer iteration (mirrors NCBI re-creating its
        // iterator at `start`).  `conv_from` is the next position to convert; the persistent
        // `rng` advances on each N, so re-visited N's get fresh CRandom values exactly as NCBI.
        let mut conv_from = j_start;

        let mut done = false;

        'inner: loop {
            if j >= j_end { break; }

            // save_masked_regions(*res, w.start(), start)
            save_masked(&mut perfect, &mut res, w.start, linker);

            while conv_from <= j + 2 { conv[conv_from] = cvt_base(bases[conv_from], rng); conv_from += 1; }
            let t = (conv[j] << 4) | (conv[j + 1] << 2) | conv[j + 2];
            j += 1; // mirrors NCBI ++it (advance before calling shift_window)

            if w.shift_window(t, &mut perfect) {
                if w.needs_processing(&thresholds) {
                    w.find_perfect(&mut perfect, &thresholds);
                }
            } else {
                // Secondary loop: window became trivial; run shift_window without find_perfect
                // until the window becomes diverse again (or sequence is exhausted).
                loop {
                    if j >= j_end { break; }
                    save_masked(&mut perfect, &mut res, w.start, linker);
                    while conv_from <= j + 2 { conv[conv_from] = cvt_base(bases[conv_from], rng); conv_from += 1; }
                    let t2 = (conv[j] << 4) | (conv[j + 1] << 2) | conv[j + 2];
                    if w.shift_window(t2, &mut perfect) {
                        done = true;
                        // Do NOT advance j here (mirrors NCBI: no ++it when done=true).
                        break;
                    }
                    j += 1;
                }
                break 'inner;
            }
        }

        // Flush remaining perfect intervals.
        let mut wstart = w.start;
        while !perfect.is_empty() {
            save_masked(&mut perfect, &mut res, wstart, linker);
            wstart += 1;
        }

        // Restart outer loop from the advanced window position (NCBI restart logic).
        let _ = done; // used implicitly via w.start advancing
        if w.start > j_start {
            j_start = w.start;
        } else {
            break 'outer;
        }
    }

    // Sort and merge (save_masked already merges incrementally, but res may be
    // slightly out of order across outer-loop restarts).
    res.sort_unstable();
    merge_intervals(res, linker)
}

/// Mask a BLASTNA-encoded sequence in place.
/// `seq` must have sentinel bytes at positions 0 and n+1 (as produced by the engine).
/// Masked positions are set to `mask_val`.
/// Returns the list of masked intervals (absolute, 0-indexed in the base array).
///
/// This is the BLAST-engine path (`dust_filter.cpp`): one fresh CRandom stream
/// per call, no N-run preprocessing.  The `dustmasker` application path is
/// [`dustmasker_intervals`].
pub fn dust_mask(
    seq: &mut [u8],
    window: usize,
    level: u32,
    linker: usize,
    mask_val: u8,
) -> Vec<(usize, usize)> {
    let n = seq.len().saturating_sub(2); // real base count
    if n < 3 { return Vec::new(); }

    let mut rng = NcbiLfg::new();
    let merged = dust_intervals(&seq[1..=n], window, level, linker, &mut rng);

    // Apply masking: interval (s, e) means mask base positions s..=e (inclusive,
    // matching NCBI's [first, second] convention).
    for &(s, e) in &merged {
        let end = e.min(n - 1);
        for pos in s..=end {
            seq[pos + 1] = mask_val;
        }
    }

    merged
}

// ── dustmasker application path (N-run splitting) ─────────────────────────────

/// Mirror of `s_FindSegmentWithLongNs` (`dust_mask_app.cpp:186`).
///
/// Reports runs of N longer than `max_ns`, plus **any** run that starts at
/// position 0 or that runs to the end of the sequence (both are reported
/// regardless of length — that asymmetry is NCBI's, not ours).
fn find_long_n_runs(bases: &[u8], max_ns: usize) -> Vec<(usize, usize)> {
    let mut out: Vec<(usize, usize)> = Vec::new();
    let mut pos = 0usize;
    let mut ns = 0usize;
    for &b in bases {
        if b == 14 {
            ns += 1;
        } else {
            if ns > 0 {
                if ns > max_ns || pos == 0 {
                    out.push((pos, pos + ns - 1));
                }
                pos += ns;
                ns = 0;
            }
            pos += 1;
        }
    }
    if ns > 0 {
        out.push((pos, pos + ns - 1));
    }
    out
}

/// Mirror of `s_InsertMerge` (`dust_mask_app.cpp:222`).
/// Note the **exact equality** test — this is not the `>=` used inside
/// `save_masked_regions`.
fn insert_merge(list: &mut Vec<(usize, usize)>, new_mask: (usize, usize), linker: usize) {
    if let Some(last) = list.last_mut() {
        if last.1 + linker == new_mask.0 {
            last.1 = new_mask.1;
            return;
        }
    }
    list.push(new_mask);
}

/// Mirror of `CSymDustMasker::operator()(seq, start, stop)`: dust `bases[start..=stop]`
/// and return absolute (offset) intervals.
fn dust_range(
    bases: &[u8],
    start: usize,
    stop: usize,
    window: usize,
    level: u32,
    linker: usize,
    rng: &mut NcbiLfg,
) -> Vec<(usize, usize)> {
    if bases.is_empty() { return Vec::new(); }
    let stop = stop.min(bases.len() - 1);
    let start = start.min(stop);
    dust_intervals(&bases[start..=stop], window, level, linker, rng)
        .into_iter()
        .map(|(s, e)| (s + start, e + start))
        .collect()
}

/// Mirror of `GetDustMasks_SkipNs` (`dust_mask_app.cpp:236`) — the masking the
/// `dustmasker` **application** performs, which differs from the BLAST filter
/// path ([`dust_mask`]): long N-runs are reported as masked and the segments
/// between them are dusted independently, sharing one CRandom stream (NCBI
/// constructs a single `CSymDustMasker` per record).
///
/// `bases` are BLASTNA codes with no sentinels; intervals are inclusive.
pub fn dustmasker_intervals(
    bases: &[u8],
    window: usize,
    level: u32,
    linker: usize,
) -> Vec<(usize, usize)> {
    let mut rng = NcbiLfg::new();
    // NOTE: the N-run threshold and the insert_merge linker use the RAW
    // parameters; only the masker itself clamps them (NCBI passes the raw
    // command-line values to GetDustMasks_SkipNs).
    let ns_ranges = find_long_n_runs(bases, window);
    if ns_ranges.is_empty() {
        return dust_intervals(bases, window, level, linker, &mut rng);
    }

    let mut rv: Vec<(usize, usize)> = Vec::new();
    let mut seq_start = 0usize;

    for &iv in &ns_ranges {
        if iv.0 == 0 {
            seq_start = iv.1 + 1;
            rv.push(iv);
            continue;
        }
        let seg = dust_range(bases, seq_start, iv.0 - 1, window, level, linker, &mut rng);
        if !seg.is_empty() {
            insert_merge(&mut rv, seg[0], linker);
            rv.extend_from_slice(&seg[1..]);
            insert_merge(&mut rv, iv, linker);
        } else {
            rv.push(iv);
        }
        seq_start = iv.1 + 1;
    }

    if seq_start < bases.len() {
        let seg = dust_range(bases, seq_start, bases.len() - 1, window, level, linker, &mut rng);
        if !seg.is_empty() {
            insert_merge(&mut rv, seg[0], linker);
            rv.extend_from_slice(&seg[1..]);
        }
    }

    rv
}

/// Wrapper around `dust_mask` that replicates the off-by-one in NCBI's
/// `TMaskedQueryRegions::RestrictToSeqInt` (`seqlocinfo.cpp`):
/// that function creates `CSeq_interval(id, from, range.GetToOpen())` where
/// `GetToOpen()` is the exclusive end, not the inclusive `to`.  The result is
/// that every restricted interval has its right endpoint incremented by 1,
/// which propagates into the per-chunk lcase_mask.
///
/// This compat applies to the MULTI-CHUNK path only (mask_query_for_alignment).
/// For single-chunk queries, NCBI does not call RestrictToSeqInt, so the
/// off-by-one does not occur.  Use plain `dust_mask` for single-chunk paths.
pub fn dust_mask_ncbi_compat(
    seq: &mut [u8],
    window: usize,
    level: u32,
    linker: usize,
    mask_val: u8,
) -> Vec<(usize, usize)> {
    let intervals = dust_mask(seq, window, level, linker, mask_val);
    let n = seq.len().saturating_sub(2); // number of real bases
    for &(_, e) in &intervals {
        let extra = e + 1;
        if extra < n {
            seq[extra + 1] = mask_val;
        }
    }
    intervals
}

fn merge_intervals(mut ivs: Vec<(usize, usize)>, linker: usize) -> Vec<(usize, usize)> {
    if ivs.is_empty() { return ivs; }
    ivs.sort_unstable();
    let mut merged: Vec<(usize, usize)> = Vec::new();
    let mut cur = ivs[0];
    for (s, e) in ivs.into_iter().skip(1) {
        if s <= cur.1 + linker {
            cur.1 = cur.1.max(e);
        } else {
            merged.push(cur);
            cur = (s, e);
        }
    }
    merged.push(cur);
    merged
}

// ── tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn encode_fasta_region(fa_path: &str, start: usize, end: usize) -> Vec<u8> {
        use std::io::BufRead;
        let f = std::fs::File::open(fa_path).expect("fasta not found");
        let reader = std::io::BufReader::new(f);
        let mut bases: Vec<u8> = Vec::new();
        for line in reader.lines() {
            let l = line.unwrap();
            if l.starts_with('>') { continue; }
            for c in l.bytes() {
                let b: u8 = match c {
                    b'A' | b'a' => 0, b'C' | b'c' => 1, b'G' | b'g' => 2, b'T' | b't' => 3, _ => 14,
                };
                bases.push(b);
                if bases.len() >= end { break; }
            }
            if bases.len() >= end { break; }
        }
        let region = &bases[start..end.min(bases.len())];
        let mut seq = vec![15u8];
        seq.extend_from_slice(region);
        seq.push(15u8);
        seq
    }

    /// Test DUST on chr22 region 10510100-10510600 and verify intervals match NCBI dustmasker.
    /// NCBI dustmasker output (0-based, inclusive): 10510356-10510364, 10510518-10510527
    #[test]
    fn test_chr22_dust_region_10510228() {
        let chr22 = "/home/rhubley/projects/Claude/rmblast-port/attic/rmblast-rs/chr22.fa";
        if !std::path::Path::new(chr22).exists() { return; }
        let region_start = 10510100usize;
        let region_end = 10510600usize;
        let mut seq = encode_fasta_region(chr22, region_start, region_end);
        let masked = dust_mask(&mut seq, DUST_WINDOW, DUST_LEVEL, DUST_LINKER, 14);
        println!("Masked intervals in chr22:{}-{} (window test):", region_start, region_end);
        for &(s, e) in &masked {
            println!("  chr22:{}-{}", region_start + s, region_start + e);
        }
        let expected: Vec<(usize, usize)> = vec![
            (10510356 - region_start, 10510364 - region_start),
            (10510518 - region_start, 10510527 - region_start),
        ];
        for (exp_s, exp_e) in &expected {
            let found = masked.iter().any(|&(s, e)| s == *exp_s && e == *exp_e);
            assert!(found, "Expected masked [{}, {}] not found. Got: {:?}",
                region_start + exp_s, region_start + exp_e,
                masked.iter().map(|&(s,e)| (region_start+s, region_start+e)).collect::<Vec<_>>());
        }
    }

    /// Test DUST on full chr22 prefix (0 to 10511000) to verify that intervals near 10510228
    /// are the same as when computed on a smaller window.
    #[test]
    fn test_chr22_dust_full_prefix_10510228() {
        let chr22 = "/home/rhubley/projects/Claude/rmblast-port/attic/rmblast-rs/chr22.fa";
        if !std::path::Path::new(chr22).exists() { return; }
        let region_end = 10511000usize;
        let mut seq = encode_fasta_region(chr22, 0, region_end);
        let masked = dust_mask(&mut seq, DUST_WINDOW, DUST_LEVEL, DUST_LINKER, 14);
        println!("Masked intervals in chr22:10510300-10510540 (full prefix test):");
        for &(s, e) in &masked {
            if s >= 10510300 && s <= 10510540 {
                println!("  chr22:{}-{}", s, e);
            }
        }
        let expected: &[(usize, usize)] = &[(10510356, 10510364), (10510518, 10510527)];
        for &(exp_s, exp_e) in expected {
            let found = masked.iter().any(|&(s, e)| s == exp_s && e == exp_e);
            assert!(found, "Expected masked chr22:{}-{} not found. Near-region intervals: {:?}",
                exp_s, exp_e,
                masked.iter().filter(|&&(s, _)| s >= 10510300 && s <= 10510540)
                    .collect::<Vec<_>>());
        }
    }

    #[test]
    fn test_low_complexity_aa_repeat() {
        let mut seq = vec![15u8];
        seq.extend(vec![0u8; 100]);
        seq.push(15);
        let masked = dust_mask(&mut seq, 64, 20, 1, 14);
        assert!(!masked.is_empty(), "pure A repeat should be masked");
    }

    /// Encode IUPACNA ASCII to BLASTNA 2-bit (A=0 C=1 G=2 T=3).
    fn encode(s: &[u8]) -> Vec<u8> {
        let mut v = vec![15u8]; // leading sentinel
        for &b in s {
            let enc = match b | 32 { // lowercase
                b'a' => 0, b'c' => 1, b'g' => 2, b't' => 3, _ => 14,
            };
            v.push(enc);
        }
        v.push(15); // trailing sentinel
        v
    }

    /// NCBI dustmasker -in alu.fa -outfmt acclist gives: >aluy  281  310
    /// Verify our implementation matches that interval.
    #[test]
    fn test_aluy_matches_ncbi() {
        let aluy = b"GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGGCGGGCGGATCACGAGGTCAGGAGATCGAGACCATCCTGGCTAACACGGTGAAACCCCGTCTCTACTAAAAATACAAAAAATTAGCCGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCGGGAGGCTGAGGCAGGAGAATGGCGTGAACCCGGGAGGCGGAGCTTGCAGTGAGCCGAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGCGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let mut seq = encode(aluy);
        let masked = dust_mask(&mut seq, 64, 20, 1, 14);
        // NCBI gives interval [281, 310].
        assert!(
            masked.iter().any(|&(s, e)| s <= 281 && e >= 310),
            "expected masked region covering [281,310], got {:?}", masked
        );
    }

    /// Compare Rust DUST intervals against NCBI dustmasker for two chr22 regions.
    /// NCBI dustmasker -outfmt acclist (0-indexed inclusive):
    ///   region1 (chr22:39030802-39031701): [366,387], [612,618], [816,889]
    ///   region2 (chr22:40060802-40061701): [580,592], [800,891]
    #[test]
    fn test_chr22_dust_vs_ncbi() {
        // Read the two region files extracted from chr22.
        let r1_seq = std::fs::read_to_string("/tmp/region1_dust.fa")
            .expect("region1 file missing — run samtools faidx chr22_full.fa chr22:39030802-39031701");
        let r2_seq = std::fs::read_to_string("/tmp/region2_dust.fa")
            .expect("region2 file missing — run samtools faidx chr22_full.fa chr22:40060802-40061701");

        let r1_bases: Vec<u8> = r1_seq.lines()
            .filter(|l| !l.starts_with('>'))
            .flat_map(|l| l.bytes())
            .collect();
        let r2_bases: Vec<u8> = r2_seq.lines()
            .filter(|l| !l.starts_with('>'))
            .flat_map(|l| l.bytes())
            .collect();

        let mut s1 = encode(&r1_bases);
        let mut s2 = encode(&r2_bases);

        let m1 = dust_mask(&mut s1, DUST_WINDOW, DUST_LEVEL, DUST_LINKER, 14);
        let m2 = dust_mask(&mut s2, DUST_WINDOW, DUST_LEVEL, DUST_LINKER, 14);

        println!("region1 Rust masked: {:?}", m1);
        println!("region2 Rust masked: {:?}", m2);

        let ncbi_r1 = vec![(366usize, 387usize), (612, 618), (816, 889)];
        let ncbi_r2 = vec![(580usize, 592usize), (800, 891)];

        assert_eq!(m1, ncbi_r1, "region1 Rust DUST differs from NCBI");
        assert_eq!(m2, ncbi_r2, "region2 Rust DUST differs from NCBI");
    }

    /// The `dustmasker` application path (`GetDustMasks_SkipNs`), which differs
    /// from the BLAST filter path: N-runs longer than `window` — plus any run at
    /// the very start or very end, whatever its length — are reported as masked,
    /// and the segments between them are dusted separately off one shared
    /// CRandom stream.
    ///
    /// Expected values are the output of NCBI dustmasker 2.17.0
    /// (`-in nstest.fa -outfmt acclist`) on exactly this sequence.
    #[test]
    fn test_dustmasker_skip_ns_vs_ncbi() {
        let mut s = String::new();
        s.push_str(&"N".repeat(5));                        // leading short N run
        s.push_str(&"ACGTTGCAAGCTTACGGATCCGTATCAGCTAGCTTAGCATCGGATTACGCA".repeat(3));
        s.push_str("CACACACACACACACACACACACACACACACACACACACA");
        s.push_str(&"GATTACAGGCTTACGGATCCGTAAGCTTACGGATCCGTATCAGCTAGCAT".repeat(2));
        s.push_str(&"N".repeat(30));                       // short internal run (<= window)
        s.push_str("TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");
        s.push_str(&"ACGTTGCAAGCTTACGGATCCGTATCAGCTAGCTTAGCATCGGATTACGCA".repeat(2));
        s.push_str(&"N".repeat(100));                      // long internal run (> window)
        s.push_str("GGCCGGGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGA");
        s.push_str("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        s.push_str("TCAGGATCCGTAAGCTTACGGATCCGTATCAGCTAGCATCGGATTACGCA");
        s.push_str(&"N".repeat(7));                        // trailing short N run

        let enc = encode(s.as_bytes());
        let bases = &enc[1..enc.len() - 1];
        let masks = dustmasker_intervals(bases, DUST_WINDOW, DUST_LEVEL, DUST_LINKER);

        let expected = vec![
            (0usize, 4usize),     // leading N run: reported despite being short
            (156, 197),           // (CA)n
            (328, 367),           // poly-T   (the 30-N run at 298..327 is NOT reported)
            (470, 569),           // long internal N run
            (619, 669),           // poly-A
            (720, 726),           // trailing N run: reported despite being short
        ];
        assert_eq!(masks, expected);
    }

    #[test]
    fn test_high_complexity_not_masked() {
        let diverse: Vec<u8> = vec![
            0,1,2,3, 0,2,1,3, 0,3,1,2, 1,0,3,2, 1,3,0,2,
            2,0,3,1, 2,1,0,3, 3,0,2,1, 3,1,0,2, 3,2,0,1,
            0,0,1,1, 2,2,3,3, 0,1,3,2, 1,2,0,3, 2,3,1,0, 3,3,2,1,
        ];
        let mut seq = vec![15u8];
        seq.extend_from_slice(&diverse);
        seq.push(15);
        let masked = dust_mask(&mut seq, 64, 20, 1, 14);
        assert!(masked.is_empty(), "diverse sequence should not be masked");
    }
}
