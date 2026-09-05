//! Build-time diagnostics: stage counters and dump hooks.
//!
//! Both macros compile to nothing unless the crate is built with the
//! `diagnostics` feature (`cargo build --release --features diagnostics`).
//! The counters are global atomics touched once per seed, and the dump hooks
//! used to read the environment once per seed. At 32 threads that cost more
//! than the alignment itself, so the default build leaves both out.

/// Add 1 (or `n`) to a diagnostic counter declared as a `static AtomicU64`.
#[macro_export]
macro_rules! diag_count {
    ($counter:path) => { $crate::diag_count!($counter, 1u64) };
    ($counter:path, $n:expr) => {{
        #[cfg(feature = "diagnostics")]
        { $counter.fetch_add(($n) as u64, ::std::sync::atomic::Ordering::Relaxed); }
        #[cfg(not(feature = "diagnostics"))]
        { let _ = &$n; }
    }};
}

/// True when the named environment variable is set. The variable is read once
/// per process, not once per call.
#[macro_export]
macro_rules! diag_enabled {
    ($name:literal) => {{
        #[cfg(feature = "diagnostics")]
        {
            static FLAG: ::std::sync::OnceLock<bool> = ::std::sync::OnceLock::new();
            *FLAG.get_or_init(|| ::std::env::var_os($name).is_some())
        }
        #[cfg(not(feature = "diagnostics"))]
        { false }
    }};
}
