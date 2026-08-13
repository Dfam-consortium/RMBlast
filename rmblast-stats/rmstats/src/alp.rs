//! Startup-time Gumbel parameter fitting via the ALP library (feature
//! "alp-fit"), the way LAST does it: for a matrix + gap-cost combination not
//! covered by the prefit table, run Sls::AlignmentEvaluer::initGapped at
//! program start and use the fitted parameters.
//!
//! The retry strategy mirrors LAST's LastEvaluer::init: on Sls::error, bump
//! the importance-sampling temperature by 0.01 and retry, up to 21 attempts.

use crate::karlin::{KarlinBlk, BLASTNA_SIZE};

/// Mirrors struct AlpFitRaw in csrc/alp_shim.cpp.
#[repr(C)]
#[derive(Debug, Clone, Copy, Default)]
pub struct AlpFitRaw {
    pub lambda: f64,
    pub lambda_error: f64,
    pub k: f64,
    pub k_error: f64,
    pub c: f64,
    pub c_error: f64,
    pub a_i: f64,
    pub a_i_error: f64,
    pub a_j: f64,
    pub a_j_error: f64,
    pub b_i: f64,
    pub b_i_error: f64,
    pub b_j: f64,
    pub b_j_error: f64,
    pub alpha_i: f64,
    pub alpha_i_error: f64,
    pub alpha_j: f64,
    pub alpha_j_error: f64,
    pub beta_i: f64,
    pub beta_i_error: f64,
    pub beta_j: f64,
    pub beta_j_error: f64,
    pub sigma: f64,
    pub sigma_error: f64,
    pub tau: f64,
    pub tau_error: f64,
    pub gapless_a: f64,
    pub gapless_a_error: f64,
    pub gapless_alpha: f64,
    pub gapless_alpha_error: f64,
    pub calc_time: f64,
}

extern "C" {
    fn rmstats_alp_init_gapped(
        alphabet_size: std::os::raw::c_long,
        scores: *const std::os::raw::c_long,
        freqs1: *const f64,
        freqs2: *const f64,
        gap_open1: std::os::raw::c_long,
        gap_epen1: std::os::raw::c_long,
        gap_open2: std::os::raw::c_long,
        gap_epen2: std::os::raw::c_long,
        insertions_after_deletions: std::os::raw::c_int,
        eps_lambda: f64,
        eps_k: f64,
        max_time: f64,
        max_mem: f64,
        rand_seed: std::os::raw::c_long,
        temperature: f64,
        deterministic: std::os::raw::c_int,
        out: *mut AlpFitRaw,
        err_buf: *mut std::os::raw::c_char,
        err_buf_len: std::os::raw::c_long,
    ) -> std::os::raw::c_int;

    fn rmstats_alp_default_temperature() -> f64;
}

/// Options for an ALP fit. The defaults reproduce the settings used for the
/// baked rmblast_*_values fits (ALP CLI defaults: eps_lambda 0.01, eps_K
/// 0.05, insertions_after_deletions false, seed 1) with a 10 s time budget.
#[derive(Debug, Clone, Copy)]
pub struct AlpFitOptions {
    pub eps_lambda: f64,
    pub eps_k: f64,
    /// Time budget in seconds. In deterministic mode this is handed to
    /// set_gapped_computation_parameters_simplified instead (LAST uses 60).
    pub max_time: f64,
    /// Memory cap in MB.
    pub max_mem: f64,
    pub rand_seed: i64,
    pub insertions_after_deletions: bool,
    /// LAST-style reproducible mode: fixed realization counts rather than a
    /// wall-clock budget.
    pub deterministic: bool,
    /// Maximum importance-sampling temperature retries (LAST uses 21
    /// attempts total).
    pub max_attempts: u32,
}

impl Default for AlpFitOptions {
    fn default() -> Self {
        AlpFitOptions {
            eps_lambda: 0.01,
            eps_k: 0.05,
            max_time: 10.0,
            max_mem: 2000.0,
            rand_seed: 1,
            insertions_after_deletions: false,
            deterministic: false,
            max_attempts: 21,
        }
    }
}

/// A completed ALP fit.
#[derive(Debug, Clone, Copy)]
pub struct AlpFit {
    pub raw: AlpFitRaw,
    /// Number of the temperature attempt that succeeded (0-based).
    pub attempt: u32,
}

impl AlpFit {
    /// The Karlin block in the NCBI mapping used by the RMBlast tables:
    /// alpha = a_J + a_I, H = lambda/alpha (so that downstream
    /// alpha = Lambda/H recovers ALP's a).
    pub fn to_karlin_blk(&self) -> KarlinBlk {
        KarlinBlk {
            lambda: self.raw.lambda,
            k: self.raw.k,
            log_k: self.raw.k.ln(),
            h: self.raw.lambda / (self.raw.a_j + self.raw.a_i),
        }
    }

    /// beta in the NCBI mapping: b_J + b_I. (The prefit table carries this
    /// per matrix; BLAST's Mode-1 length adjustment approximates beta = 0,
    /// so use this only if you deliberately want the fitted intercept.)
    pub fn ncbi_beta(&self) -> f64 {
        self.raw.b_j + self.raw.b_i
    }
}

/// Errors from a fit.
#[derive(Debug)]
pub enum AlpFitError {
    /// All temperature attempts raised Sls::error; message of the last one.
    SlsError(String),
    /// A non-Sls exception escaped ALP.
    UnknownException,
    /// initGapped returned but the evaluer was not in a good state.
    NotGood,
}

impl std::fmt::Display for AlpFitError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AlpFitError::SlsError(m) => write!(f, "ALP error: {m}"),
            AlpFitError::UnknownException => write!(f, "ALP: unknown exception"),
            AlpFitError::NotGood => write!(f, "ALP: evaluer not in a good state"),
        }
    }
}

impl std::error::Error for AlpFitError {}

/// Fit Gumbel/FSC parameters for an `n x n` scoring system with per-sequence
/// backgrounds `freqs1`/`freqs2` and NCBI-convention affine gap costs
/// (a length-k gap costs open + k*extend; 1 = insertions, 2 = deletions).
pub fn fit_gumbel_gapped(
    scores: &[Vec<i64>],
    freqs1: &[f64],
    freqs2: &[f64],
    gap_open1: i32,
    gap_extend1: i32,
    gap_open2: i32,
    gap_extend2: i32,
    opts: &AlpFitOptions,
) -> Result<AlpFit, AlpFitError> {
    let n = scores.len();
    assert!(n > 0 && freqs1.len() == n && freqs2.len() == n);
    let mut flat: Vec<std::os::raw::c_long> = Vec::with_capacity(n * n);
    for row in scores {
        assert_eq!(row.len(), n);
        flat.extend(row.iter().map(|&v| v as std::os::raw::c_long));
    }

    // ALP uses process-global state (RNG, importance-sampling scratch);
    // concurrent fits corrupt it (observed: njn_localmaxstatutil assertion
    // abort).  Serialize every fit process-wide.
    static ALP_FFI_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());
    let _guard = ALP_FFI_LOCK.lock().unwrap_or_else(|p| p.into_inner());

    let base_temp = unsafe { rmstats_alp_default_temperature() };
    let mut last_err = AlpFitError::UnknownException;
    for attempt in 0..opts.max_attempts {
        let temperature = base_temp + 0.01 * attempt as f64;
        let mut out = AlpFitRaw::default();
        let mut err_buf = [0i8; 512];
        let status = unsafe {
            rmstats_alp_init_gapped(
                n as std::os::raw::c_long,
                flat.as_ptr(),
                freqs1.as_ptr(),
                freqs2.as_ptr(),
                gap_open1 as std::os::raw::c_long,
                gap_extend1 as std::os::raw::c_long,
                gap_open2 as std::os::raw::c_long,
                gap_extend2 as std::os::raw::c_long,
                opts.insertions_after_deletions as std::os::raw::c_int,
                opts.eps_lambda,
                opts.eps_k,
                opts.max_time,
                opts.max_mem,
                opts.rand_seed as std::os::raw::c_long,
                temperature,
                opts.deterministic as std::os::raw::c_int,
                &mut out,
                err_buf.as_mut_ptr() as *mut std::os::raw::c_char,
                err_buf.len() as std::os::raw::c_long,
            )
        };
        match status {
            0 => return Ok(AlpFit { raw: out, attempt }),
            1 => {
                let msg = unsafe {
                    std::ffi::CStr::from_ptr(err_buf.as_ptr() as *const std::os::raw::c_char)
                }
                .to_string_lossy()
                .into_owned();
                last_err = AlpFitError::SlsError(msg);
                // Sls::error: retry with a higher temperature (LAST-style).
            }
            3 => return Err(AlpFitError::NotGood),
            _ => return Err(AlpFitError::UnknownException),
        }
    }
    Err(last_err)
}

/// Convenience for rmblastn: fit from a 16x16 blastna-encoded matrix and its
/// `# FREQS` background (blastna order: A=0, C=1, G=2, T=3), extracting the
/// 4x4 A/C/G/T core and using the same frequencies for both sequences —
/// exactly how the prefit rmblast_*_values were produced.
pub fn fit_gumbel_blastna(
    matrix: &[[i32; BLASTNA_SIZE]; BLASTNA_SIZE],
    freqs: &[f64; BLASTNA_SIZE],
    gap_open: i32,
    gap_extend: i32,
    opts: &AlpFitOptions,
) -> Result<AlpFit, AlpFitError> {
    let scores: Vec<Vec<i64>> = (0..4)
        .map(|i| (0..4).map(|j| matrix[i][j] as i64).collect())
        .collect();
    let freqs4: Vec<f64> = freqs[..4].to_vec();
    fit_gumbel_gapped(
        &scores,
        &freqs4,
        &freqs4,
        gap_open,
        gap_extend,
        gap_open,
        gap_extend,
        opts,
    )
}
