//! Special functions for the t-distribution used by eBayes moderation.
//!
//! Student-t tail uses the regularized incomplete beta via Lentz continued
//! fraction. digamma/trigamma/trigamma_inverse live in rsomics-ebayes-core.

pub fn ln_gamma(x: f64) -> f64 {
    // Lanczos g=607/128, n=15
    const G: f64 = 607.0 / 128.0;
    const C: [f64; 15] = [
        0.999_999_999_999_997_1,
        57.156_235_665_862_92,
        -59.597_960_355_475_49,
        14.136_097_974_741_746,
        -0.491_913_816_097_620_2,
        3.399_464_998_481_189e-5,
        4.652_362_892_704_858e-5,
        -9.837_447_530_487_956e-5,
        1.580_887_032_249_125e-4,
        -2.102_644_417_241_049e-4,
        2.174_396_181_152_126e-4,
        -1.643_181_065_367_639e-4,
        8.441_822_398_385_275e-5,
        -2.610_724_774_107_823e-5,
        3.689_918_265_953_162e-6,
    ];
    if x < 0.5 {
        // reflection
        std::f64::consts::PI.ln() - (std::f64::consts::PI * x).sin().ln() - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut a = C[0];
        let t = x + G + 0.5;
        for (i, &c) in C.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

/// Two-sided Student-t p-value: P(|T_df| >= |t|).
pub fn t_pvalue_two_sided(t: f64, df: f64) -> f64 {
    if df.is_infinite() {
        return 2.0 * normal_sf(t.abs());
    }
    let x = df / (df + t * t);
    betai(df / 2.0, 0.5, x)
}

/// Upper-tail Student-t quantile: the q such that P(T_df > q) = p, p in (0,1).
pub fn t_quantile_upper(p: f64, df: f64) -> f64 {
    if p <= 0.0 {
        return f64::INFINITY;
    }
    if p >= 1.0 {
        return f64::NEG_INFINITY;
    }
    // upper tail prob of t>0 is 0.5*two_sided. Bisection on |t|.
    let target = p;
    let mut lo = 0.0f64;
    let mut hi = 1.0f64;
    while 0.5 * t_pvalue_two_sided(hi, df) > target {
        hi *= 2.0;
        if hi > 1e12 {
            break;
        }
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        let pm = 0.5 * t_pvalue_two_sided(mid, df);
        if pm > target {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-13 * hi.max(1.0) {
            break;
        }
    }
    0.5 * (lo + hi)
}

fn normal_sf(x: f64) -> f64 {
    0.5 * erfc(x / std::f64::consts::SQRT_2)
}

fn erfc(x: f64) -> f64 {
    // Numerical Recipes erfcc, fractional error < 1.2e-7; refined by one
    // Newton-free rational pass is unnecessary at the precision we diff.
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let ans = t
        * (-z * z - 1.265_512_23
            + t * (1.000_023_68
                + t * (0.374_091_96
                    + t * (0.096_784_18
                        + t * (-0.186_288_06
                            + t * (0.278_868_07
                                + t * (-1.135_203_98
                                    + t * (1.488_515_87
                                        + t * (-0.822_152_23 + t * 0.170_872_77)))))))))
            .exp();
    if x >= 0.0 { ans } else { 2.0 - ans }
}

/// Regularized incomplete beta I_x(a,b) via Lentz continued fraction.
fn betai(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let bt = (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    if x < (a + 1.0) / (a + b + 2.0) {
        bt * betacf(a, b, x) / a
    } else {
        1.0 - bt * betacf(b, a, 1.0 - x) / b
    }
}

fn betacf(a: f64, b: f64, x: f64) -> f64 {
    const FPMIN: f64 = 1e-300;
    let qab = a + b;
    let qap = a + 1.0;
    let qam = a - 1.0;
    let mut c = 1.0;
    let mut d = 1.0 - qab * x / qap;
    if d.abs() < FPMIN {
        d = FPMIN;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..=300 {
        let m = m as f64;
        let m2 = 2.0 * m;
        let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        h *= d * c;
        let aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = 1.0 + aa / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-15 {
            break;
        }
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;
    use rsomics_ebayes_core::{digamma, trigamma, trigamma_inverse};

    #[test]
    fn digamma_known() {
        assert!((digamma(1.0) + 0.577_215_664_901_532_9).abs() < 1e-10);
        assert!((digamma(0.5) + 1.963_510_026_021_423).abs() < 1e-9);
    }

    #[test]
    fn trigamma_known() {
        assert!((trigamma(1.0) - std::f64::consts::PI.powi(2) / 6.0).abs() < 1e-10);
    }

    #[test]
    fn trigamma_inverse_roundtrip() {
        for &y in &[0.7, 1.3, 4.0, 20.0, 100.0] {
            let x = trigamma(y);
            let yi = trigamma_inverse(x);
            assert!((yi - y).abs() / y < 1e-6, "y={y} got {yi}");
        }
    }

    #[test]
    fn t_quantile_roundtrip() {
        for &(p, df) in &[(0.025, 10.0), (0.125, 17.89), (0.001, 30.0)] {
            let q = t_quantile_upper(p, df);
            let back = 0.5 * t_pvalue_two_sided(q, df);
            assert!(
                (back - p).abs() / p < 1e-6,
                "p={p} df={df} q={q} back={back}"
            );
        }
    }

    #[test]
    fn t_pvalue_matches_known() {
        // P(|T_10| >= 2.228) ~ 0.05
        let p = t_pvalue_two_sided(2.228, 10.0);
        assert!((p - 0.05).abs() < 1e-3, "p={p}");
    }
}
