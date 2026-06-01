//! Per-gene least-squares fit (limma lmFit, method="ls") and contrasts.fit.
//!
//! Householder QR of the design X once, then solve for every gene. The QR also
//! yields R, from which (X'X)^-1 = R^-1 R^-T gives the unscaled coefficient
//! covariance shared by all genes (unweighted case).

use rsomics_common::{Result, RsomicsError};

pub struct Qr {
    n: usize,
    p: usize,
    /// packed Householder vectors + R, column-major-ish row-major [n][p]
    qr: Vec<Vec<f64>>,
    rdiag: Vec<f64>,
}

impl Qr {
    pub fn new(x: &[Vec<f64>]) -> Result<Qr> {
        let n = x.len();
        let p = x[0].len();
        if n < p {
            return Err(RsomicsError::InvalidInput(format!(
                "design has {n} samples < {p} coefficients (rank-deficient)"
            )));
        }
        let mut qr: Vec<Vec<f64>> = x.to_vec();
        let mut rdiag = vec![0.0; p];
        for k in 0..p {
            let mut nrm = 0.0f64;
            for row in qr.iter().take(n).skip(k) {
                nrm = nrm.hypot(row[k]);
            }
            if nrm == 0.0 {
                return Err(RsomicsError::InvalidInput(
                    "design matrix is rank-deficient".into(),
                ));
            }
            if qr[k][k] < 0.0 {
                nrm = -nrm;
            }
            for row in qr.iter_mut().take(n).skip(k) {
                row[k] /= nrm;
            }
            qr[k][k] += 1.0;
            for j in (k + 1)..p {
                let mut s = 0.0;
                for row in qr.iter().take(n).skip(k) {
                    s += row[k] * row[j];
                }
                s = -s / qr[k][k];
                for row in qr.iter_mut().take(n).skip(k) {
                    let add = s * row[k];
                    row[j] += add;
                }
            }
            rdiag[k] = -nrm;
        }
        Ok(Qr { n, p, qr, rdiag })
    }

    /// Apply Q' to y in place (length n), returns nothing — y mutated.
    #[allow(clippy::needless_range_loop)]
    fn qty(&self, y: &mut [f64]) {
        for k in 0..self.p {
            let mut s = 0.0;
            for i in k..self.n {
                s += self.qr[i][k] * y[i];
            }
            s = -s / self.qr[k][k];
            for i in k..self.n {
                y[i] += s * self.qr[i][k];
            }
        }
    }

    /// Solve for coefficients given y; returns (beta[p], rss).
    pub fn solve(&self, y: &[f64]) -> (Vec<f64>, f64) {
        let mut qty = y.to_vec();
        self.qty(&mut qty);
        let rss: f64 = qty[self.p..].iter().map(|&e| e * e).sum();
        let mut beta = vec![0.0; self.p];
        for j in (0..self.p).rev() {
            beta[j] = qty[j];
            for k in (j + 1)..self.p {
                beta[j] -= self.qr[j][k] * beta[k];
            }
            beta[j] /= self.rdiag[j];
        }
        (beta, rss)
    }

    /// (X'X)^-1 = R^-1 R^-T, the p×p unscaled coefficient covariance.
    #[allow(clippy::needless_range_loop)]
    pub fn xtx_inv(&self) -> Vec<Vec<f64>> {
        let p = self.p;
        // R is upper triangular: diag = rdiag, off-diag R[i][j] = qr[i][j] (i<j)
        let r_at =
            |i: usize, j: usize| -> f64 { if i == j { self.rdiag[i] } else { self.qr[i][j] } };
        // invert R (upper triangular) -> rinv
        let mut rinv = vec![vec![0.0; p]; p];
        for i in 0..p {
            rinv[i][i] = 1.0 / r_at(i, i);
        }
        for j in 0..p {
            for i in (0..j).rev() {
                let mut s = 0.0;
                for k in (i + 1)..=j {
                    s += r_at(i, k) * rinv[k][j];
                }
                rinv[i][j] = -s / r_at(i, i);
            }
        }
        // (X'X)^-1 = rinv * rinv'
        let mut cov = vec![vec![0.0; p]; p];
        for i in 0..p {
            for j in 0..p {
                let mut s = 0.0;
                for (ri, rj) in rinv[i].iter().zip(&rinv[j]) {
                    s += ri * rj;
                }
                cov[i][j] = s;
            }
        }
        cov
    }
}

pub struct Fit {
    pub coef_names: Vec<String>,
    /// [gene][coef]
    pub coefficients: Vec<Vec<f64>>,
    /// per-coef unscaled sd = sqrt(diag (X'X)^-1) — shared across genes
    pub stdev_unscaled: Vec<f64>,
    /// residual sd per gene
    pub sigma: Vec<f64>,
    /// residual df = n - p (shared)
    pub df_residual: f64,
    pub amean: Vec<f64>,
    pub genes: Vec<String>,
    pub samples_n: usize,
}

pub fn lm_fit(
    y: &[Vec<f64>],
    genes: &[String],
    x: &[Vec<f64>],
    coef_names: &[String],
) -> Result<Fit> {
    let n = x.len();
    let p = x[0].len();
    if y.iter().any(|row| row.len() != n) {
        return Err(RsomicsError::InvalidInput(
            "expression samples do not match design rows".into(),
        ));
    }
    let df_residual = (n - p) as f64;
    if df_residual < 1.0 {
        return Err(RsomicsError::InvalidInput(
            "residual degrees of freedom < 1 (need more samples than coefficients)".into(),
        ));
    }
    let qr = Qr::new(x)?;
    let cov = qr.xtx_inv();
    let stdev_unscaled: Vec<f64> = (0..p).map(|j| cov[j][j].sqrt()).collect();

    let mut coefficients = Vec::with_capacity(y.len());
    let mut sigma = Vec::with_capacity(y.len());
    let mut amean = Vec::with_capacity(y.len());
    for row in y {
        let (beta, rss) = qr.solve(row);
        coefficients.push(beta);
        sigma.push((rss / df_residual).sqrt());
        amean.push(row.iter().sum::<f64>() / n as f64);
    }

    Ok(Fit {
        coef_names: coef_names.to_vec(),
        coefficients,
        stdev_unscaled,
        sigma,
        df_residual,
        amean,
        genes: genes.to_vec(),
        samples_n: n,
    })
}

/// contrasts.fit: transform a coefficient-space fit into contrast space.
/// new coef = C' beta; new unscaled var diag = diag(C' (X'X)^-1 C).
#[allow(clippy::needless_range_loop)]
pub fn contrasts_fit(fit: &Fit, contrast: &crate::matrix::Contrast, xtx_inv: &[Vec<f64>]) -> Fit {
    let p = fit.coef_names.len();
    let q = contrast.names.len();
    let cmat = &contrast.c; // [p][q]

    let mut new_coef = Vec::with_capacity(fit.coefficients.len());
    for beta in &fit.coefficients {
        let mut nb = vec![0.0; q];
        for (col, item) in nb.iter_mut().enumerate() {
            let mut s = 0.0;
            for (i, &b) in beta.iter().enumerate() {
                s += cmat[i][col] * b;
            }
            *item = s;
        }
        new_coef.push(nb);
    }

    // C' (X'X)^-1 C, take diagonal -> sqrt = new stdev.unscaled per contrast
    let mut stdev = vec![0.0; q];
    for (col, sd) in stdev.iter_mut().enumerate() {
        let mut acc = 0.0;
        for i in 0..p {
            for j in 0..p {
                acc += cmat[i][col] * xtx_inv[i][j] * cmat[j][col];
            }
        }
        *sd = acc.sqrt();
    }

    Fit {
        coef_names: contrast.names.clone(),
        coefficients: new_coef,
        stdev_unscaled: stdev,
        sigma: fit.sigma.clone(),
        df_residual: fit.df_residual,
        amean: fit.amean.clone(),
        genes: fit.genes.clone(),
        samples_n: fit.samples_n,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_group_means() {
        // design: intercept + group; samples A A B B
        let x = vec![
            vec![1.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![1.0, 1.0],
        ];
        let names = vec!["Int".to_string(), "Grp".to_string()];
        let y = vec![vec![1.0, 3.0, 5.0, 7.0]];
        let genes = vec!["g".to_string()];
        let f = lm_fit(&y, &genes, &x, &names).unwrap();
        // group A mean = 2 (intercept), group B mean = 6, diff = 4
        assert!((f.coefficients[0][0] - 2.0).abs() < 1e-9);
        assert!((f.coefficients[0][1] - 4.0).abs() < 1e-9);
    }
}
