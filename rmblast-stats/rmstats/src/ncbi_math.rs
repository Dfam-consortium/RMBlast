//! Ports of helpers from ncbi_math.c used by the statistics code.

/// ln(2) as defined in ncbi_math.h (NCBIMATH_LN2).
pub const NCBIMATH_LN2: f64 = 0.69314718055994530941723212145818;

/// BLAST_Gcd (ncbi_math.c): greatest common divisor.
pub fn blast_gcd(a: i32, b: i32) -> i32 {
    let mut a = a;
    let mut b = b.abs();
    if b > a {
        std::mem::swap(&mut a, &mut b);
    }
    while b != 0 {
        let c = a % b;
        a = b;
        b = c;
    }
    a
}

/// BLAST_Nint (ncbi_math.c): round to nearest integer, half away from zero.
pub fn blast_nint(x: f64) -> i64 {
    let x = if x >= 0.0 { x + 0.5 } else { x - 0.5 };
    x as i64
}

/// BLAST_Powi (ncbi_math.c): integer power by repeated squaring.
pub fn blast_powi(x: f64, n: i32) -> f64 {
    if n == 0 {
        return 1.0;
    }
    let mut x = x;
    let mut n = n;
    if x == 0.0 {
        if n < 0 {
            return f64::INFINITY; // HUGE_VAL
        }
        return 0.0;
    }
    if n < 0 {
        x = 1.0 / x;
        n = -n;
    }
    let mut y = 1.0;
    while n > 0 {
        if n & 1 != 0 {
            y *= x;
        }
        n /= 2;
        x *= x;
    }
    y
}

/// BLAST_Expm1 (ncbi_math.c): exp(x) - 1, accurate for small |x|.
pub fn blast_expm1(x: f64) -> f64 {
    let absx = x.abs();
    if absx > 0.33 {
        return x.exp() - 1.0;
    }
    if absx < 1.0e-16 {
        return x;
    }
    x * (1.0
        + x * (1.0 / 2.0
            + x * (1.0 / 6.0
                + x * (1.0 / 24.0
                    + x * (1.0 / 120.0
                        + x * (1.0 / 720.0
                            + x * (1.0 / 5040.0
                                + x * (1.0 / 40320.0
                                    + x * (1.0 / 362880.0
                                        + x * (1.0 / 3628800.0
                                            + x * (1.0 / 39916800.0
                                                + x * (1.0 / 479001600.0
                                                    + x / 6227020800.0))))))))))))
}
