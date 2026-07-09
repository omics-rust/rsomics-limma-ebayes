//! Empirical-Bayes moderation (Smyth 2004, Stat Appl Genet Mol Biol 3:1,
//! doi:10.2202/1544-6115.1027).
//!
//! squeezeVar: moment-estimate the scaled-inverse-chisquare prior (d0, s0^2)
//! from the per-gene residual variances, then shrink each toward s0^2.
//! Moderated t = beta / (unscaled_sd * sqrt(s2.post)); p two-sided on df.total.
//! var.prior for the B-statistic comes from a t-mixture estimator over the
//! largest moderated statistics, bounded by stdev.coef.lim^2 / median(s2.prior).

use crate::fit::Fit;
use crate::special::{t_pvalue_two_sided, t_quantile_upper};
use rsomics_ebayes_core::fit_f_dist;

pub struct Moderated {
    /// [gene][coef]
    pub t: Vec<Vec<f64>>,
    pub p: Vec<Vec<f64>>,
    pub lods: Vec<Vec<f64>>,
    pub df_total: f64,
    pub df_prior: f64,
    pub s2_prior: f64,
}

fn squeeze_var(sigma2: &[f64], df: f64) -> (Vec<f64>, f64, f64) {
    let (s20, df0) = fit_f_dist(sigma2, df);
    let s2_post: Vec<f64> = if df0.is_infinite() {
        sigma2.iter().map(|_| s20).collect()
    } else {
        sigma2
            .iter()
            .map(|&s2| (df0 * s20 + df * s2) / (df0 + df))
            .collect()
    };
    (s2_post, df0, s20)
}

/// var.prior for one coefficient column (limma tmixture.vector). df is the
/// shared df.total; equal across genes here, so the MaxDF rescale is a no-op.
/// Each top-probe v0 is bounded by v0_lim before averaging.
fn tmixture_column(
    tstat: &[f64],
    stdev_unscaled: f64,
    df: f64,
    proportion: f64,
    v0_lim: (f64, f64),
) -> f64 {
    let ngenes = tstat.len();
    let ntarget = (proportion / 2.0 * ngenes as f64).ceil() as usize;
    if ntarget < 1 {
        return f64::NAN;
    }
    let p = (ntarget as f64 / ngenes as f64).max(proportion);

    let mut at: Vec<f64> = tstat.iter().map(|t| t.abs()).collect();
    at.sort_by(|a, b| b.partial_cmp(a).unwrap());

    let mut v0_sum = 0.0;
    let v1 = stdev_unscaled * stdev_unscaled;
    for (i, &t) in at.iter().take(ntarget).enumerate() {
        let r = (i + 1) as f64;
        let p0 = t_pvalue_two_sided(t, df); // 2*pt(t, df, upper) = two-sided p
        let ptarget = ((r - 0.5) / ngenes as f64 - (1.0 - p) * p0) / p;
        let mut v0 = 0.0;
        if ptarget > p0 {
            let qtarget = t_quantile_upper(ptarget / 2.0, df);
            v0 = (v1 * ((t / qtarget).powi(2) - 1.0)).max(0.0);
        }
        v0_sum += v0.clamp(v0_lim.0, v0_lim.1);
    }
    v0_sum / ntarget as f64
}

pub fn ebayes(fit: &Fit, xtx_diag_unscaled: &[f64], proportion: f64) -> Moderated {
    let sigma2: Vec<f64> = fit.sigma.iter().map(|s| s * s).collect();
    let (s2_post, df_prior, s2_prior) = squeeze_var(&sigma2, fit.df_residual);

    let ng = fit.coefficients.len();
    let q = fit.coef_names.len();

    let df_pooled = ng as f64 * fit.df_residual;
    let df_total = (fit.df_residual + df_prior).min(df_pooled);

    let mut t = vec![vec![0.0; q]; ng];
    let mut p = vec![vec![0.0; q]; ng];
    for gi in 0..ng {
        let post_sd = s2_post[gi].sqrt();
        for cj in 0..q {
            let tv = fit.coefficients[gi][cj] / (xtx_diag_unscaled[cj] * post_sd);
            t[gi][cj] = tv;
            p[gi][cj] = t_pvalue_two_sided(tv, df_total);
        }
    }

    // var.prior.lim = stdev.coef.lim^2 / median(s2.prior); s2.prior is scalar
    // here (no trend), so its median is itself.
    let v0_lim = (0.1f64 * 0.1 / s2_prior, 4.0f64 * 4.0 / s2_prior);

    let mut lods = vec![vec![0.0; q]; ng];
    let const_term = (proportion / (1.0 - proportion)).ln();
    for cj in 0..q {
        let col_t: Vec<f64> = (0..ng).map(|gi| t[gi][cj]).collect();
        let mut v0 = tmixture_column(&col_t, xtx_diag_unscaled[cj], df_total, proportion, v0_lim);
        if v0.is_nan() {
            v0 = 1.0 / s2_prior;
        }

        let u2 = xtx_diag_unscaled[cj].powi(2);
        let r = (u2 + v0) / u2;
        for gi in 0..ng {
            let t2 = t[gi][cj].powi(2);
            let kernel = if df_prior > 1e6 {
                t2 * (1.0 - 1.0 / r) / 2.0
            } else {
                (1.0 + df_total) / 2.0 * ((t2 + df_total) / (t2 / r + df_total)).ln()
            };
            lods[gi][cj] = const_term - r.ln() / 2.0 + kernel;
        }
    }

    Moderated {
        t,
        p,
        lods,
        df_total,
        df_prior,
        s2_prior,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rsomics_ebayes_core::{digamma, trigamma, trigamma_inverse};

    #[test]
    fn fit_f_dist_infinite_prior_uses_arithmetic_mean() {
        // Equal gene variances -> evar<=0 -> df.prior=Inf. limma sets
        // s2.prior = mean(sigma^2). The old exp(emean) link gave ~12.04 here.
        let s2 = vec![3.38_f64; 20];
        let (s20, df2) = fit_f_dist(&s2, 1.0);
        assert!(df2.is_infinite());
        assert!((s20 - 3.38).abs() < 1e-12, "s20={s20}");
    }

    #[test]
    fn fit_f_dist_finite_prior_link_unchanged() {
        // Spread of variances -> evar>0 -> finite df.prior, digamma-link scale.
        let s2 = vec![1.0, 2.0, 4.0, 0.5, 8.0, 3.0, 6.0, 1.5, 0.7, 5.0];
        let (s20, df2) = fit_f_dist(&s2, 4.0);
        assert!(df2.is_finite() && df2 > 0.0, "df2={df2}");
        assert!(s20 > 0.0, "s20={s20}");
    }

    #[test]
    fn fit_f_dist_clamps_zero_variance_gene() {
        // A single exact-zero variance among a well-conditioned set must be
        // floored (1e-5*median), not dropped, so the moment fit still has spread
        // and df.prior stays finite. Dropping it collapsed df.prior to Inf.
        let mut s2 = vec![
            0.2000665, 0.2039878, 0.2678043, 0.2933900, 0.3335739, 0.55, 0.6, 0.7, 0.8, 0.9, 1.0,
            1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
        ];
        s2[0] = 0.0;
        let (s20, df2) = fit_f_dist(&s2, 6.0);
        assert!(
            df2.is_finite() && df2 > 0.0,
            "df2 must be finite, got {df2}"
        );
        assert!(s20 > 0.0 && s20 < 1.0, "s20={s20}");
    }

    #[test]
    fn fit_f_dist_single_usable_gene_gives_ordinary_t() {
        // nok==1 -> df2=0 (scale=that variance): moderation vanishes, moderated
        // t reduces to the ordinary t. Matches limma on a one-gene fit.
        let (s20, df2) = fit_f_dist(&[0.42], 6.0);
        assert_eq!(df2, 0.0);
        assert!((s20 - 0.42).abs() < 1e-15, "s20={s20}");
    }

    #[test]
    fn fit_f_dist_clamp_is_noop_on_well_conditioned() {
        // Every variance already exceeds 1e-5*median, so the floor changes
        // nothing and the finite-prior link is bit-identical to the unclamped fit.
        let s2 = vec![1.0, 2.0, 4.0, 0.5, 8.0, 3.0, 6.0, 1.5, 0.7, 5.0];
        let (s20, df2) = fit_f_dist(&s2, 4.0);
        let e: Vec<f64> = s2
            .iter()
            .map(|v| v.ln() - digamma(2.0) + 2.0f64.ln())
            .collect();
        let n = e.len() as f64;
        let emean = e.iter().sum::<f64>() / n;
        let evar = e.iter().map(|v| (v - emean).powi(2)).sum::<f64>() / (n - 1.0) - trigamma(2.0);
        let df2_ref = 2.0 * trigamma_inverse(evar);
        let s20_ref = (emean + digamma(df2_ref / 2.0) - (df2_ref / 2.0).ln()).exp();
        assert!((df2 - df2_ref).abs() < 1e-12, "df2={df2} ref={df2_ref}");
        assert!((s20 - s20_ref).abs() < 1e-12, "s20={s20} ref={s20_ref}");
    }
}
