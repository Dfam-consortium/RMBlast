//! Fixed-size diagonal deduplication hash table, porting NCBI's BLAST_DiagHash.
//!
//! NCBI uses a 512-bucket hash with a dynamically-grown chain array (starting at
//! 256 entries).  The hash function is the Knuth multiplicative hash
//! `(diag as u32).wrapping_mul(0x9E370001) % 512`.
//!
//! Diagonal convention: `diag = s_off.wrapping_sub(q_off)` (NCBI sign, i32).
//!
//! When all chain slots for a bucket are occupied AND none is stale
//! (stale = current s_off has advanced past the stored s_end), a new slot is
//! appended to the chain Vec (doubling threshold mirrors NCBI's realloc strategy).
//! When a stale slot IS found it is evicted and recycled — this is the intentional
//! "leakiness" that matches NCBI's behavior in highly repetitive regions.

pub const NUM_BUCKETS: usize = 512;
const INIT_CHAIN_CAP: usize = 257; // cell 0 unused (null sentinel) + 256 usable

struct DiagHashCell {
    diag:  i32,
    level: u32, // s_end of last extension (0 = empty slot)
    next:  u32, // next cell index; 0 = end of chain
}

pub struct BlastDiagHash {
    backbone:   Box<[u32; NUM_BUCKETS]>,
    chain:      Vec<DiagHashCell>,
    occupancy:  u32, // next free index (starts at 1; cell 0 is the null sentinel)
}

impl BlastDiagHash {
    pub fn new() -> Self {
        let mut chain: Vec<DiagHashCell> = Vec::with_capacity(INIT_CHAIN_CAP);
        for _ in 0..INIT_CHAIN_CAP {
            chain.push(DiagHashCell { diag: 0, level: 0, next: 0 });
        }
        BlastDiagHash {
            backbone:  Box::new([0u32; NUM_BUCKETS]),
            chain,
            occupancy: 1,
        }
    }

    /// Reset to empty without deallocating.
    pub fn clear(&mut self) {
        self.backbone.fill(0);
        self.occupancy = 1;
        // chain cells are not cleared; they will be overwritten on first use via backbone
    }

    #[inline]
    fn bucket(diag: i32) -> usize {
        (diag as u32).wrapping_mul(0x9E370001) as usize % NUM_BUCKETS
    }

    /// Return the stored s_end for `diag` (0 if no entry).
    #[inline]
    pub fn get(&self, diag: i32) -> u32 {
        let mut idx = self.backbone[Self::bucket(diag)] as usize;
        while idx != 0 {
            let cell = &self.chain[idx];
            if cell.diag == diag {
                return cell.level;
            }
            idx = cell.next as usize;
        }
        0
    }

    /// Record `s_end` as the new high-water mark for `diag`.
    ///
    /// `s_off` is the current seed's subject offset (FWD space).
    /// `window` matches NCBI's effective stale threshold: `window_size + Delta + 1`
    ///   where for rmblastn (window_size=0, scan_range=0, word_length=w):
    ///   `Delta = MIN(0, -w)` → `window = 1 - w` (e.g. -13 for w=14).
    /// Stale condition (NCBI sign convention): `s_off - cell.level > window`
    ///   i.e. (s_off as i32 - cell.level as i32) > window.
    #[inline]
    pub fn insert(&mut self, diag: i32, s_end: u32, s_off: u32, window: i32) {
        let bucket = Self::bucket(diag);
        let mut idx = self.backbone[bucket] as usize;

        while idx != 0 {
            let cell = &mut self.chain[idx];
            if cell.diag == diag {
                cell.level = s_end;
                return;
            }
            // Stale: NCBI condition: (s_off as i32) - (cell.level as i32) > window.
            if (s_off as i32).wrapping_sub(cell.level as i32) > window {
                cell.diag  = diag;
                cell.level = s_end;
                return;
            }
            idx = cell.next as usize;
        }

        // No matching or stale slot; append a new cell.
        let new_idx = self.occupancy as usize;
        if new_idx >= self.chain.len() {
            self.chain.push(DiagHashCell { diag: 0, level: 0, next: 0 });
        }
        {
            let cell = &mut self.chain[new_idx];
            cell.diag  = diag;
            cell.level = s_end;
            cell.next  = self.backbone[bucket];
        }
        self.backbone[bucket] = new_idx as u32;
        self.occupancy += 1;
    }
}

impl Default for BlastDiagHash {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    /// With nothing ever stale (window i32::MAX), the table behaves as a map
    /// from diagonal to the last s_end stored, past the initial chain
    /// capacity; clear() empties it.
    #[test]
    fn behaves_as_a_map_when_nothing_is_stale() {
        let mut h = BlastDiagHash::new();
        let mut want: HashMap<i32, u32> = HashMap::new();
        let mut x = 0x9E3779B97F4A7C15u64;
        for i in 0..20_000u32 {
            x ^= x << 13; x ^= x >> 7; x ^= x << 17;
            let diag = (x % 3000) as i32 - 1500;
            let s_end = 1 + i;
            h.insert(diag, s_end, 0, i32::MAX);
            want.insert(diag, s_end);
        }
        assert!(want.len() > INIT_CHAIN_CAP, "test must grow the chain");
        for d in -1600..1600 {
            assert_eq!(h.get(d), want.get(&d).copied().unwrap_or(0), "diag {d}");
        }
        h.clear();
        for d in -1600..1600 {
            assert_eq!(h.get(d), 0, "diag {d} after clear");
        }
        h.insert(7, 42, 0, i32::MAX);
        assert_eq!(h.get(7), 42);
    }

    /// NCBI's stale rule (s_BlastDiagHashInsert): a cell in the same bucket
    /// whose level the scan has passed by more than `window`
    /// (s_off - level > window) is reused for the new diagonal, which drops
    /// the old one; otherwise a new cell is chained.
    #[test]
    fn stale_cell_in_the_bucket_is_reused() {
        let d1 = 5i32;
        let d2 = (6..).find(|&d| BlastDiagHash::bucket(d) == BlastDiagHash::bucket(d1)).unwrap();
        let window = 1 - 14; // word size 14

        // not stale: s_off - level = 80 - 100 = -20, not > -13
        let mut h = BlastDiagHash::new();
        h.insert(d1, 100, 50, window);
        h.insert(d2, 300, 80, window);
        assert_eq!((h.get(d1), h.get(d2)), (100, 300));

        // stale: 200 - 100 = 100 > -13, so d2 takes d1's cell
        let mut h = BlastDiagHash::new();
        h.insert(d1, 100, 50, window);
        h.insert(d2, 300, 200, window);
        assert_eq!((h.get(d1), h.get(d2)), (0, 300));

        // at the threshold: 87 - 100 = -13, not > -13
        let mut h = BlastDiagHash::new();
        h.insert(d1, 100, 50, window);
        h.insert(d2, 300, 87, window);
        assert_eq!((h.get(d1), h.get(d2)), (100, 300));

        // the same diagonal is updated in place, stale or not
        let mut h = BlastDiagHash::new();
        h.insert(d1, 100, 50, window);
        h.insert(d1, 500, 10_000, window);
        assert_eq!(h.get(d1), 500);
    }
}
