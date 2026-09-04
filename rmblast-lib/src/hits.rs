//! HSP (High-scoring Segment Pair) data structures.

/// Strand of the subject (database) sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Strand {
    Plus,
    Minus,
}

impl Strand {
    pub fn as_str(self) -> &'static str {
        match self {
            Strand::Plus => "+",
            Strand::Minus => "-",
        }
    }
}

/// Edit operation in the traceback — using NCBI's eGapAlign naming conventions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EditOp {
    /// Aligned pair (match or mismatch) — query and subject both advance.
    Sub,
    /// Gap in query (subject advances alone).  NCBI eGapAlignDel.
    /// Alignment column: Query='-', Subject=base.
    GapInQuery,
    /// Gap in subject (query advances alone).  NCBI eGapAlignIns.
    /// Alignment column: Query=base, Subject='-'.
    GapInSubject,
}

/// A run-length encoded edit script (traceback).
#[derive(Debug, Clone, Default)]
pub struct EditScript {
    pub ops: Vec<(EditOp, u32)>,
}

impl EditScript {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn push(&mut self, op: EditOp, count: u32) {
        if count == 0 { return; }
        if let Some(last) = self.ops.last_mut() {
            if last.0 == op {
                last.1 += count;
                return;
            }
        }
        self.ops.push((op, count));
    }

    /// Total aligned columns.
    pub fn align_len(&self) -> u32 {
        self.ops.iter().map(|(_, n)| n).sum()
    }

    /// Number of query bases consumed (Sub + GapInSubject).
    pub fn query_consumed(&self) -> u32 {
        self.ops.iter()
            .filter(|(op, _)| *op != EditOp::GapInQuery)
            .map(|(_, n)| n)
            .sum()
    }

    /// Number of subject bases consumed (Sub + GapInQuery).
    pub fn subject_consumed(&self) -> u32 {
        self.ops.iter()
            .filter(|(op, _)| *op != EditOp::GapInSubject)
            .map(|(_, n)| n)
            .sum()
    }

    /// Reverse the order of ops in place (used to flip a left-extension script).
    pub fn reverse(&mut self) {
        self.ops.reverse();
    }
}

/// A single aligned segment (HSP).
#[derive(Debug, Clone)]
pub struct Hsp {
    pub score: i32,

    /// 0-based start offset in the query (first aligned base).
    pub q_start: u32,
    /// 0-based exclusive end offset in the query.
    pub q_end: u32,
    /// Total query sequence length.
    pub q_len: u32,

    /// 0-based start offset in the subject (on the plus strand).
    pub s_start: u32,
    /// 0-based exclusive end offset in the subject (on the plus strand).
    pub s_end: u32,
    /// Total subject sequence length.
    pub s_len: u32,

    pub strand: Strand,

    /// Full traceback edit script.
    pub edit_script: EditScript,

    /// Query bases in alignment (BLASTNA; gaps encoded as 15).
    pub q_seq: Vec<u8>,
    /// Subject bases in alignment (BLASTNA; gaps encoded as 15).
    pub s_seq: Vec<u8>,
}

impl Hsp {
    /// Query span as NCBI stores it: on the minus strand the offsets are
    /// positions in the reverse-complemented query.
    pub fn ncbi_query_span(&self) -> (u32, u32) {
        match self.strand {
            Strand::Plus => (self.q_start, self.q_end),
            Strand::Minus => (self.q_len - self.q_end, self.q_len - self.q_start),
        }
    }
}

/// NCBI's `ScoreCompareHSPs` (blast_hits.c): score DESC, subject offset ASC,
/// subject end DESC, query offset ASC, query end DESC, with the query offsets
/// taken in NCBI's strand-context coordinates.
pub fn score_compare_hsps(a: &Hsp, b: &Hsp) -> std::cmp::Ordering {
    let (q_off_a, q_end_a) = a.ncbi_query_span();
    let (q_off_b, q_end_b) = b.ncbi_query_span();
    b.score.cmp(&a.score)
        .then_with(|| a.s_start.cmp(&b.s_start))
        .then_with(|| b.s_end.cmp(&a.s_end))
        .then_with(|| q_off_a.cmp(&q_off_b))
        .then_with(|| q_end_b.cmp(&q_end_a))
}
