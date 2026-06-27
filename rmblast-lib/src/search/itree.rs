/// Interval tree for O(log N) HSP containment checks.
///
/// Faithful port of NCBI `blast_itree.c` (Jason Papadopoulos) for the rmblastn
/// use case.  Mirrors the **2D `eQueryAndSubject` structure**: a query-axis
/// interval tree where each node's midpoint is *itself* an interval tree indexing
/// subject offsets (the "midpoint tree"); HSPs that straddle the center of a
/// subject node collect into a singly-linked list at that node.  This 2D structure
/// is essential — a flat midpoint list checks *all* straddling HSPs and finds
/// containers that NCBI's subject-axis descent skips, producing over-containment in
/// dense low-complexity regions (see PORTING notes / bug investigation).
///
///   - Single query strand per tree (caller maintains separate Plus/Minus trees).
///   - The containment check includes NCBI's MB_HSP_CLOSE diagonal-proximity test
///     (min_diag_separation = 50 for rmblastn).
///   - `add` additionally performs endpoint-sharing deduplication mirroring NCBI's
///     `eQueryAndSubject` add: before inserting a prelim HSP, the tree checks whether
///     any existing entry already occupies the same left endpoint (q_start, s_start)
///     or right endpoint (q_end, s_end).  Higher/equal existing score → drop the new
///     HSP; lower → logically delete the old entry and insert the new one.
///
/// All coordinates passed to `add`/`add_simple`/`contains` must be in **normalized
/// (plus-strand) subject space**; the caller is responsible for the conversion.
///
/// ## Node layout
///
/// Nodes live in a flat `Vec` (arena).  Node 0 is always the root and is never
/// used as a child, so the value `0` is used as a null pointer (`NIL`).
/// Leaf nodes have `hsp_idx != INTERNAL`; internal nodes have `hsp_idx == INTERNAL`.
/// A node's `leftend`/`rightend` describe a **query** range for query-axis nodes and
/// a **subject** range for subject-axis (midpoint-tree) nodes.  An internal query
/// node's `midptr` points at the root of its subject midpoint tree; a subject node's
/// `midptr` is the head of a LIFO linked list of leaves straddling the subject center.

use std::collections::HashMap;

const INTERNAL: u32 = u32::MAX; // sentinel: this node is an internal (non-leaf) node
const NIL: u32 = 0;             // null pointer; node 0 = root, never a child

/// A single node in the arena.
#[derive(Clone, Default)]
struct ITreeNode {
    leftend:  u32, // left endpoint of the (query or subject) region this node covers
    rightend: u32, // right endpoint of the region this node covers
    leftptr:  u32, // index of left child, or NIL
    midptr:   u32, // subject midpoint-tree root (query node) / list head (subject node), or NIL
    rightptr: u32, // index of right child, or NIL
    hsp_idx:  u32, // INTERNAL for internal nodes; index into `hsps` for leaf nodes
}

/// The HSP payload stored at each leaf.
#[derive(Clone)]
pub struct ITreeHsp {
    pub q_start: u32, // query start (absolute)
    pub q_end:   u32, // query end   (absolute)
    pub s_start: u32, // subject start (normalized, plus-strand)
    pub s_end:   u32, // subject end   (normalized, plus-strand)
    pub score:   i32,
}

/// Centered interval tree indexed by query offset, then subject offset.
pub struct BlastIntervalTree {
    nodes: Vec<ITreeNode>,
    hsps:  Vec<ITreeHsp>,
    /// Parallel deleted-flag vector; `deleted[i]` is true when `hsps[i]` was
    /// superseded by a higher-scoring entry sharing its left or right endpoint.
    deleted: Vec<bool>,
    /// Left-endpoint registry: (q_start, s_start) → (score, hsp_idx).
    left_endpoints:  HashMap<(u32, u32), (i32, u32)>,
    /// Right-endpoint registry: (q_end, s_end) → (score, hsp_idx).
    right_endpoints: HashMap<(u32, u32), (i32, u32)>,
    /// Subject-axis bounds for the nested midpoint trees ([s_min, s_max]).
    s_min: u32,
    s_max: u32,
    /// When true, the query axis uses RC coordinates (q_start = L - FWD_q_end,
    /// q_end = L - FWD_q_start). The diagonal formula is (q - s) instead of (q + s),
    /// which is mathematically equivalent to NCBI's MB_HSP_CLOSE in RC space.
    /// Callers must convert coordinates before calling add_simple/contains.
    pub rc_query: bool,
}

impl BlastIntervalTree {
    /// Create an empty tree whose root covers query range `[q_start, q_end]` and
    /// whose nested subject trees span `[s_min, s_max]`.  Pass `q_start = 0`,
    /// `q_end = query_length + 1`, `s_min = 0`, `s_max = subject_length + 1`
    /// (matching NCBI's `Blast_IntervalTreeInit` convention).
    pub fn new(q_start: u32, q_end: u32, s_min: u32, s_max: u32) -> Self {
        let root = ITreeNode {
            leftend:  q_start,
            rightend: q_end,
            leftptr:  NIL,
            midptr:   NIL,
            rightptr: NIL,
            hsp_idx:  INTERNAL,
        };
        BlastIntervalTree {
            nodes: vec![root],
            hsps:  Vec::new(),
            deleted: Vec::new(),
            left_endpoints:  HashMap::new(),
            right_endpoints: HashMap::new(),
            s_min,
            s_max,
            rc_query: false,
        }
    }

    /// Like `new` but with `rc_query = true` — minus-strand Phase 2b containment tree.
    pub fn new_rc(q_start: u32, q_end: u32, s_min: u32, s_max: u32) -> Self {
        let mut t = Self::new(q_start, q_end, s_min, s_max);
        t.rc_query = true;
        t
    }

    /// Reset to empty without releasing memory (mirrors `Blast_IntervalTreeReset`).
    pub fn reset(&mut self) {
        self.nodes.truncate(1);
        let root = &mut self.nodes[0];
        root.leftptr  = NIL;
        root.midptr   = NIL;
        root.rightptr = NIL;
        root.hsp_idx  = INTERNAL;
        self.hsps.clear();
        self.deleted.clear();
        self.left_endpoints.clear();
        self.right_endpoints.clear();
    }

    /// Allocate an internal node covering `[le, re]`; returns its index.
    #[inline]
    fn alloc_internal(&mut self, le: u32, re: u32) -> u32 {
        let idx = self.nodes.len() as u32;
        self.nodes.push(ITreeNode {
            leftend:  le,
            rightend: re,
            leftptr:  NIL,
            midptr:   NIL,
            rightptr: NIL,
            hsp_idx:  INTERNAL,
        });
        idx
    }

    /// Push the HSP and create its leaf node; returns the leaf node index.
    #[inline]
    fn alloc_leaf(&mut self, hsp: ITreeHsp) -> u32 {
        let hsp_idx = self.hsps.len() as u32;
        self.hsps.push(hsp);
        self.deleted.push(false);
        let new_idx = self.nodes.len() as u32;
        self.nodes.push(ITreeNode {
            leftend:  0,
            rightend: 0,
            leftptr:  NIL,
            midptr:   NIL,
            rightptr: NIL,
            hsp_idx,
        });
        new_idx
    }

    /// Insert `hsp` into the tree (mirrors `BlastIntervalTreeAddHSP`, eQueryAndSubject:
    /// endpoint-sharing dedup + 2D query×subject structural indexing).
    ///
    /// Returns `true` if the HSP was inserted, `false` if it was blocked by a
    /// higher-scoring entry sharing the same left or right endpoint.
    pub fn add(&mut self, hsp: ITreeHsp, plus_strand: bool) -> bool {
        let (left_ep, right_ep) = if plus_strand {
            ((hsp.q_start, hsp.s_start), (hsp.q_end, hsp.s_end))
        } else {
            // Minus strand: NCBI's "left endpoint" is at the alignment's 3'-end of
            // the FWD query (= q_end in FWD coords) paired with FWD subject start.
            ((hsp.q_end, hsp.s_start), (hsp.q_start, hsp.s_end))
        };
        let new_score = hsp.score;
        let new_qlen  = hsp.q_end - hsp.q_start;
        let new_slen  = hsp.s_end - hsp.s_start;

        // Endpoint-sharing dedup, mirroring s_HSPsHaveCommonEndpoint: when the new
        // HSP shares an endpoint with an existing one, keep the higher score; on an
        // equal score keep the SHORTER HSP (shorter query range, then shorter subject
        // range); on a full tie keep the existing (tree) HSP.  NCBI applies the left
        // check first (removing a worse existing entry even if the new HSP is later
        // blocked at the right endpoint), so we delete eagerly before the right check.
        //
        // Left endpoint — s_IntervalTreeHasHSPEndpoint(eIntervalTreeLeft).
        if let Some(&(ex_score, ex_idx)) = self.left_endpoints.get(&left_ep) {
            let (ex_qlen, ex_slen) = {
                let e = &self.hsps[ex_idx as usize];
                (e.q_end - e.q_start, e.s_end - e.s_start)
            };
            if !Self::endpoint_new_wins(new_score, new_qlen, new_slen, ex_score, ex_qlen, ex_slen) {
                return false; // existing entry is better or identical — discard new HSP
            }
            // NCBI removes the losing HSP's whole leaf from the tree; mirror that by
            // dropping BOTH of its endpoint entries, not just the matched one.  Leaving
            // the other endpoint behind lets a deleted HSP spuriously block a later HSP.
            self.remove_endpoints(ex_idx, plus_strand);
        }
        // Right endpoint — s_IntervalTreeHasHSPEndpoint(eIntervalTreeRight).
        if let Some(&(ex_score, ex_idx)) = self.right_endpoints.get(&right_ep) {
            let (ex_qlen, ex_slen) = {
                let e = &self.hsps[ex_idx as usize];
                (e.q_end - e.q_start, e.s_end - e.s_start)
            };
            if !Self::endpoint_new_wins(new_score, new_qlen, new_slen, ex_score, ex_qlen, ex_slen) {
                return false;
            }
            self.remove_endpoints(ex_idx, plus_strand);
        }

        let new_idx = self.alloc_leaf(hsp);
        let hsp_idx = self.nodes[new_idx as usize].hsp_idx;
        self.left_endpoints.insert(left_ep,  (new_score, hsp_idx));
        self.right_endpoints.insert(right_ep, (new_score, hsp_idx));

        self.insert_structural(new_idx);
        true
    }

    /// Mark `idx` deleted and drop both of its endpoint-map entries, mirroring NCBI's
    /// structural removal of a losing leaf in `s_IntervalTreeHasHSPEndpoint`.  Each
    /// entry is removed only if it still points at `idx` (a newer HSP may have taken
    /// over an endpoint key), so we never clobber a live entry.
    fn remove_endpoints(&mut self, idx: u32, plus_strand: bool) {
        self.deleted[idx as usize] = true;
        let e = &self.hsps[idx as usize];
        let (left_ep, right_ep) = if plus_strand {
            ((e.q_start, e.s_start), (e.q_end, e.s_end))
        } else {
            ((e.q_end, e.s_start), (e.q_start, e.s_end))
        };
        if matches!(self.left_endpoints.get(&left_ep), Some(&(_, i)) if i == idx) {
            self.left_endpoints.remove(&left_ep);
        }
        if matches!(self.right_endpoints.get(&right_ep), Some(&(_, i)) if i == idx) {
            self.right_endpoints.remove(&right_ep);
        }
    }

    /// Mirrors NCBI `s_HSPsHaveCommonEndpoint`'s "keep best" rule for two HSPs that
    /// share an endpoint: higher score wins; on an equal score the SHORTER HSP wins
    /// (shorter query range, then shorter subject range); a full tie keeps the
    /// existing (tree) HSP.  Returns `true` iff the new HSP should replace the existing.
    #[inline]
    fn endpoint_new_wins(
        new_score: i32, new_qlen: u32, new_slen: u32,
        ex_score:  i32, ex_qlen:  u32, ex_slen:  u32,
    ) -> bool {
        if new_score != ex_score { return new_score > ex_score; }
        if new_qlen  != ex_qlen  { return new_qlen  < ex_qlen;  }
        if new_slen  != ex_slen  { return new_slen  < ex_slen;  }
        false
    }

    /// Insert `hsp` without endpoint-sharing deduplication (Phase 2b/2c trees).
    pub fn add_simple(&mut self, hsp: ITreeHsp) {
        let new_idx = self.alloc_leaf(hsp);
        self.insert_structural(new_idx);
    }

    /// Structural insertion of an already-allocated leaf node, mirroring the
    /// `BlastIntervalTreeAddHSP` body (eQueryAndSubject): descend the query-axis
    /// tree; when the HSP straddles a query-node center, switch to (and descend) the
    /// node's nested subject-axis tree.  Leaf-leaf collisions split into a new
    /// internal node; a straddling old leaf is re-homed into a fresh subject tree.
    fn insert_structural(&mut self, new_idx: u32) {
        let hsp_idx = self.nodes[new_idx as usize].hsp_idx;
        let q_off = self.hsps[hsp_idx as usize].q_start;
        let q_end = self.hsps[hsp_idx as usize].q_end;
        let s_off = self.hsps[hsp_idx as usize].s_start;
        let s_end = self.hsps[hsp_idx as usize].s_end;

        // region_* is the query range until the HSP first straddles a query center,
        // after which it becomes the subject range (index_subject_range = true).
        let mut region_start = q_off;
        let mut region_end   = q_end;
        let mut index_subject_range = false;
        let mut root: u32 = 0;

        loop {
            let r_le = self.nodes[root as usize].leftend;
            let r_re = self.nodes[root as usize].rightend;
            let middle = (r_le as u64 + r_re as u64) / 2;

            let which_half_left: bool;
            let old_index: u32;

            if (region_end as u64) < middle {
                let left = self.nodes[root as usize].leftptr;
                if left == NIL {
                    self.nodes[root as usize].leftptr = new_idx;
                    return;
                }
                if self.nodes[left as usize].hsp_idx == INTERNAL {
                    root = left;
                    continue;
                }
                old_index = left;
                which_half_left = true;
            } else if (region_start as u64) > middle {
                let right = self.nodes[root as usize].rightptr;
                if right == NIL {
                    self.nodes[root as usize].rightptr = new_idx;
                    return;
                }
                if self.nodes[right as usize].hsp_idx == INTERNAL {
                    root = right;
                    continue;
                }
                old_index = right;
                which_half_left = false;
            } else {
                // The new interval straddles the center of this node.
                if index_subject_range {
                    // Already indexing subject offsets: prepend to the midpoint list.
                    let old_mid = self.nodes[root as usize].midptr;
                    self.nodes[new_idx as usize].midptr = old_mid;
                    self.nodes[root as usize].midptr = new_idx;
                    return;
                } else {
                    // Switch to indexing the subject range: descend into (creating if
                    // needed) the nested subject tree rooted at this node's midptr.
                    index_subject_range = true;
                    if self.nodes[root as usize].midptr == NIL {
                        let mid = self.alloc_internal(self.s_min, self.s_max);
                        self.nodes[root as usize].midptr = mid;
                    }
                    root = self.nodes[root as usize].midptr;
                    region_start = s_off;
                    region_end = s_end;
                    continue;
                }
            }

            // Two leaves want the same subtree: create an internal node covering the
            // chosen half, attach it to the parent, re-home the old leaf inside it,
            // then loop to place the new leaf within the new internal node.
            let (mid_le, mid_re) = if which_half_left {
                (r_le, middle as u32)
            } else {
                (middle as u32 + 1, r_re)
            };
            let mid_index = self.alloc_internal(mid_le, mid_re);
            if which_half_left {
                self.nodes[root as usize].leftptr = mid_index;
            } else {
                self.nodes[root as usize].rightptr = mid_index;
            }

            let old_hsp_idx = self.nodes[old_index as usize].hsp_idx;
            let (old_rs, old_re) = if index_subject_range {
                (self.hsps[old_hsp_idx as usize].s_start, self.hsps[old_hsp_idx as usize].s_end)
            } else {
                (self.hsps[old_hsp_idx as usize].q_start, self.hsps[old_hsp_idx as usize].q_end)
            };

            root = mid_index;
            let mmid = (mid_le as u64 + mid_re as u64) / 2;
            if (old_re as u64) < mmid {
                self.nodes[mid_index as usize].leftptr = old_index;
            } else if (old_rs as u64) > mmid {
                self.nodes[mid_index as usize].rightptr = old_index;
            } else if index_subject_range {
                // Old leaf straddles the subject center: goes in the midpoint list.
                self.nodes[mid_index as usize].midptr = old_index;
            } else {
                // Old leaf straddles the query center while we are still indexing
                // query offsets: allocate a fresh subject tree just for it.
                let os = self.hsps[old_hsp_idx as usize].s_start;
                let oe = self.hsps[old_hsp_idx as usize].s_end;
                let mid2 = self.alloc_internal(self.s_min, self.s_max);
                self.nodes[mid_index as usize].midptr = mid2;
                let m2mid = (self.s_min as u64 + self.s_max as u64) / 2;
                if (oe as u64) < m2mid {
                    self.nodes[mid2 as usize].leftptr = old_index;
                } else if (os as u64) > m2mid {
                    self.nodes[mid2 as usize].rightptr = old_index;
                } else {
                    self.nodes[mid2 as usize].midptr = old_index;
                }
            }
            // Loop again (root = mid_index) to place the new leaf.
        }
    }

    /// Return `true` if any tree HSP contains the candidate (spatial containment +
    /// diagonal proximity), mirroring `BlastIntervalTreeContainsHSP`.
    ///
    /// `plus_strand` selects the + or − diagonal formula for the MB_HSP_CLOSE check.
    /// `min_diag_sep` is `MIN_DIAG_SEP` (50 for rmblastn).
    pub fn contains(
        &self,
        cand_q_start: u32,
        cand_q_end:   u32,
        cand_s_start: u32,
        cand_s_end:   u32,
        cand_score:   i32,
        plus_strand:  bool,
        min_diag_sep: i64,
    ) -> bool {
        let region_start = cand_q_start;
        let region_end   = cand_q_end;
        let mut node: u32 = 0; // start at root

        loop {
            // Internal node: first test its nested subject midpoint tree, then descend.
            let mid = self.nodes[node as usize].midptr;
            if mid != NIL
                && self.midpoint_tree_contains(
                    mid, cand_q_start, cand_q_end, cand_s_start, cand_s_end,
                    cand_score, plus_strand, min_diag_sep,
                )
            {
                return true;
            }

            let leftend  = self.nodes[node as usize].leftend;
            let rightend = self.nodes[node as usize].rightend;
            let middle   = (leftend as u64 + rightend as u64) / 2;

            let next = if (region_end as u64) < middle {
                self.nodes[node as usize].leftptr
            } else if (region_start as u64) > middle {
                self.nodes[node as usize].rightptr
            } else {
                // Candidate straddles this center.  Any container must also straddle
                // it and thus live in the (already-checked) subject midpoint tree.
                return false;
            };

            if next == NIL {
                return false;
            }
            // A leaf reached directly off the query axis (no subject midpoint tree):
            // test it and finish.
            if self.nodes[next as usize].hsp_idx != INTERNAL {
                return self.hsp_is_contained(
                    cand_q_start, cand_q_end, cand_s_start, cand_s_end,
                    cand_score, plus_strand, min_diag_sep, self.nodes[next as usize].hsp_idx,
                );
            }
            node = next;
        }
    }

    /// Descend a nested subject-axis midpoint tree rooted at `root`, mirroring NCBI
    /// `s_MidpointTreeContainsHSP`.  At each node the straddling-leaf linked list is
    /// scanned; descent then follows the candidate's *subject* range.  HSPs in
    /// subject subtrees off the candidate's path are intentionally NOT examined —
    /// this is the behavior that prevents NCBI's over-containment.
    #[allow(clippy::too_many_arguments)]
    fn midpoint_tree_contains(
        &self,
        root:         u32,
        cand_q_start: u32,
        cand_q_end:   u32,
        cand_s_start: u32,
        cand_s_end:   u32,
        cand_score:   i32,
        plus_strand:  bool,
        min_diag_sep: i64,
    ) -> bool {
        let region_start = cand_s_start;
        let region_end   = cand_s_end;
        let mut node = root;

        loop {
            if self.nodes[node as usize].hsp_idx != INTERNAL {
                // Reached a leaf of the subject tree.
                return self.hsp_is_contained(
                    cand_q_start, cand_q_end, cand_s_start, cand_s_end,
                    cand_score, plus_strand, min_diag_sep, self.nodes[node as usize].hsp_idx,
                );
            }

            // Scan the straddling-leaf linked list at this subject node.
            let mut mid = self.nodes[node as usize].midptr;
            while mid != NIL {
                if self.hsp_is_contained(
                    cand_q_start, cand_q_end, cand_s_start, cand_s_end,
                    cand_score, plus_strand, min_diag_sep, self.nodes[mid as usize].hsp_idx,
                ) {
                    return true;
                }
                mid = self.nodes[mid as usize].midptr;
            }

            let leftend  = self.nodes[node as usize].leftend;
            let rightend = self.nodes[node as usize].rightend;
            let middle   = (leftend as u64 + rightend as u64) / 2;

            let next = if (region_end as u64) < middle {
                self.nodes[node as usize].leftptr
            } else if (region_start as u64) > middle {
                self.nodes[node as usize].rightptr
            } else {
                return false;
            };

            if next == NIL {
                return false;
            }
            node = next;
        }
    }

    /// Core containment predicate: does tree HSP `hsp_idx` contain the candidate?
    /// Mirrors `s_HSPIsContained` in blast_itree.c.
    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn hsp_is_contained(
        &self,
        cand_q_start: u32,
        cand_q_end:   u32,
        cand_s_start: u32,
        cand_s_end:   u32,
        cand_score:   i32,
        plus_strand:  bool,
        min_diag_sep: i64,
        hsp_idx:      u32,
    ) -> bool {
        if self.deleted[hsp_idx as usize] {
            return false;
        }
        let h = &self.hsps[hsp_idx as usize];

        if cand_score > h.score {
            return false;
        }
        // Spatial: candidate's [q_start,q_end]×[s_start,s_end] inside tree HSP's box
        // (CONTAINED_IN_HSP for both the start and end endpoints).
        if !(h.q_start <= cand_q_start && h.q_end >= cand_q_end
            && h.s_start <= cand_s_start && h.s_end >= cand_s_end)
        {
            return false;
        }
        // MB_HSP_CLOSE: at least one endpoint pair within min_diag_sep diagonals.
        if plus_strand {
            let d0 = ((cand_q_start as i64 - cand_s_start as i64)
                - (h.q_start as i64 - h.s_start as i64)).abs();
            let d1 = ((cand_q_end as i64 - cand_s_end as i64)
                - (h.q_end as i64 - h.s_end as i64)).abs();
            d0 < min_diag_sep || d1 < min_diag_sep
        } else if self.rc_query {
            // Minus strand with RC query coords (q_start = L - FWD_q_end, etc.).
            let d0 = ((cand_q_start as i64 - cand_s_start as i64)
                - (h.q_start as i64 - h.s_start as i64)).abs();
            let d1 = ((cand_q_end as i64 - cand_s_end as i64)
                - (h.q_end as i64 - h.s_end as i64)).abs();
            d0 < min_diag_sep || d1 < min_diag_sep
        } else {
            // Minus strand with FWD coords: anti-diagonal = q_FWD_end + s_start.
            let d0 = ((cand_q_end as i64 + cand_s_start as i64)
                - (h.q_end as i64 + h.s_start as i64)).abs();
            let d1 = ((cand_q_start as i64 + cand_s_end as i64)
                - (h.q_start as i64 + h.s_end as i64)).abs();
            d0 < min_diag_sep || d1 < min_diag_sep
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Regression: deleting an HSP via one endpoint must purge BOTH of its endpoint
    /// entries, else the stale entry spuriously blocks a later HSP from the tree.
    ///
    /// Reproduces chr20:56,317,895 × 6kbHsap satellite (Phase 2a prelim sequence).
    /// Adding 355b shares the RIGHT endpoint (340,4182) with 355a and wins (shorter),
    /// so 355a is removed.  The 353 prelim shares 355a's now-defunct LEFT endpoint
    /// (289,4132); before the fix the stale entry blocked 353, dropping the container
    /// that suppresses a redundant ungapped hit and yielding a spurious extra HSP.
    #[test]
    fn deleted_hsp_does_not_block_via_stale_endpoint() {
        const MDS: i64 = 50; // MIN_DIAG_SEP
        let mut t = BlastIntervalTree::new(0, 602, 0, 6019);
        let hsp = |qs, qe, ss, se, sc| ITreeHsp { q_start: qs, q_end: qe, s_start: ss, s_end: se, score: sc };

        assert!(t.add(hsp(295, 340, 4125, 4170, 367), true));
        assert!(t.add(hsp(289, 340, 4132, 4182, 355), true)); // 355a
        // 355b shares 355a's right endpoint (340,4182), is shorter, and wins → 355a removed.
        assert!(t.add(hsp(292, 340, 4130, 4182, 355), true)); // 355b
        // 353 shares 355a's stale left endpoint (289,4132). 355a is gone, so 353 must be added.
        assert!(t.add(hsp(289, 342, 4132, 4191, 353), true), "353 prelim must enter the tree");

        // The redundant ungapped hit must now be reported contained within the 353 box.
        assert!(
            t.contains(292, 338, 4138, 4184, 340, true, MDS),
            "ungapped hit q[292,338] s[4138,4184] should be contained in the 353 prelim"
        );
    }
}
