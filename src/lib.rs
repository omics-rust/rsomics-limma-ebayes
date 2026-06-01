//! lmFit + eBayes + topTable(sort.by="none") for a log-expression matrix.
//!
//! Clean-room reimplementation of limma's moderated t-statistic pipeline:
//! Smyth (2004), Stat Appl Genet Mol Biol 3(1):3, doi:10.2202/1544-6115.1027.
//! No limma (GPL) source was consulted; the method follows the published
//! paper and is validated black-box against the binary.

mod ebayes;
mod fit;
mod matrix;
mod special;

use std::io::{BufWriter, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

pub use matrix::{read_contrast, read_design, read_expr};

pub struct Options<'a> {
    pub expr: &'a Path,
    pub design: &'a Path,
    pub contrast: Option<&'a Path>,
    /// 1-based coefficient (or contrast) to tabulate; defaults to the last.
    pub coef: Option<usize>,
    pub proportion: f64,
}

/// topTable(sort.by="none") row, in input gene order.
pub struct Row {
    pub gene: String,
    pub logfc: f64,
    pub ave_expr: f64,
    pub t: f64,
    pub p_value: f64,
    pub adj_p_val: f64,
    pub b: f64,
}

pub struct Results {
    pub coef_name: String,
    pub rows: Vec<Row>,
    pub df_total: f64,
    pub df_prior: f64,
    pub s2_prior: f64,
}

pub fn run(opts: &Options) -> Result<Results> {
    let expr = read_expr(opts.expr)?;
    let design = read_design(opts.design)?;
    if design.x.len() != expr.samples.len() {
        return Err(RsomicsError::InvalidInput(format!(
            "design has {} rows, expression has {} samples",
            design.x.len(),
            expr.samples.len()
        )));
    }

    let base = fit::lm_fit(&expr.y, &expr.genes, &design.x, &design.coef_names)?;
    let xtx_inv = fit::Qr::new(&design.x)?.xtx_inv();

    let fit = if let Some(cpath) = opts.contrast {
        let contrast = read_contrast(cpath, &design.coef_names)?;
        fit::contrasts_fit(&base, &contrast, &xtx_inv)
    } else {
        base
    };

    let m = ebayes::ebayes(&fit, &fit.stdev_unscaled, opts.proportion);

    let nc = fit.coef_names.len();
    let coef = opts.coef.unwrap_or(nc);
    if coef < 1 || coef > nc {
        return Err(RsomicsError::InvalidInput(format!(
            "--coef {coef} out of range 1..={nc}"
        )));
    }
    let k = coef - 1;

    let pvals: Vec<f64> = (0..fit.genes.len()).map(|gi| m.p[gi][k]).collect();
    let adj = bh_adjust(&pvals);

    let mut rows = Vec::with_capacity(fit.genes.len());
    for (gi, gene) in fit.genes.iter().enumerate() {
        rows.push(Row {
            gene: gene.clone(),
            logfc: fit.coefficients[gi][k],
            ave_expr: fit.amean[gi],
            t: m.t[gi][k],
            p_value: m.p[gi][k],
            adj_p_val: adj[gi],
            b: m.lods[gi][k],
        });
    }

    Ok(Results {
        coef_name: fit.coef_names[k].clone(),
        rows,
        df_total: m.df_total,
        df_prior: m.df_prior,
        s2_prior: m.s2_prior,
    })
}

/// Benjamini-Hochberg adjusted p-values, returned in input order.
fn bh_adjust(p: &[f64]) -> Vec<f64> {
    let n = p.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| p[b].partial_cmp(&p[a]).unwrap());
    let mut adj = vec![0.0; n];
    let mut cummin = f64::INFINITY;
    for (rank, &i) in idx.iter().enumerate() {
        let m = (n - rank) as f64;
        let v = (n as f64 / m * p[i]).min(1.0);
        cummin = cummin.min(v);
        adj[i] = cummin;
    }
    adj
}

pub fn write_results(res: &Results, out: &mut dyn Write) -> Result<()> {
    let mut w = BufWriter::with_capacity(1 << 20, out);
    writeln!(w, "gene\tlogFC\tAveExpr\tt\tP.Value\tadj.P.Val\tB").map_err(RsomicsError::Io)?;
    let mut fmt = ryu::Buffer::new();
    let mut line = String::with_capacity(128);
    for row in &res.rows {
        line.clear();
        line.push_str(&row.gene);
        for v in [
            row.logfc,
            row.ave_expr,
            row.t,
            row.p_value,
            row.adj_p_val,
            row.b,
        ] {
            line.push('\t');
            line.push_str(fmt.format(v));
        }
        line.push('\n');
        w.write_all(line.as_bytes()).map_err(RsomicsError::Io)?;
    }
    w.flush().map_err(RsomicsError::Io)?;
    Ok(())
}
