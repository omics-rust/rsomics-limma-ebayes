//! Differential compat against limma lmFit + eBayes + topTable.
//!
//! - `golden_*` always runs: ours vs a committed R-captured topTable.
//! - `live_r_*` runs only when an Rscript with limma is found; it regenerates
//!   the oracle and diffs against ours (loud-skip otherwise).

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::process::Command;

const EPS: f64 = 1e-6; // relative

fn ours() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-limma-ebayes"))
}

fn golden(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn manifest(rel: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel)
}

type Table = (Vec<String>, BTreeMap<String, Vec<f64>>);

fn parse(text: &str) -> Table {
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(str::to_string)
        .collect();
    let mut rows = BTreeMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let gene = f.next().unwrap().to_string();
        let vals: Vec<f64> = f.map(|s| s.trim().parse().unwrap()).collect();
        rows.insert(gene, vals);
    }
    (header, rows)
}

fn assert_close(a: &Table, b: &Table, label: &str) {
    assert_eq!(a.0, b.0, "{label}: header mismatch");
    assert_eq!(a.1.len(), b.1.len(), "{label}: row count mismatch");
    let mut max_rel = 0.0f64;
    for (gene, x) in &a.1 {
        let y =
            b.1.get(gene)
                .unwrap_or_else(|| panic!("{label}: missing gene {gene}"));
        assert_eq!(x.len(), y.len(), "{label}: {gene} width mismatch");
        for (vx, vy) in x.iter().zip(y) {
            let rel = (vx - vy).abs() / vy.abs().max(1e-9);
            max_rel = max_rel.max(rel);
            assert!(rel < EPS, "{label}: {gene} ours={vx} ref={vy} rel={rel:e}");
        }
    }
    eprintln!("{label}: max relative deviation = {max_rel:e}");
}

fn run_ours(coef: usize) -> String {
    let out = Command::new(ours())
        .arg(golden("expr.tsv"))
        .args(["--design", golden("design.tsv").to_str().unwrap()])
        .args(["--coef", &coef.to_string()])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "ours failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

#[test]
fn golden_coef2() {
    let ours_out = run_ours(2);
    let expected = std::fs::read_to_string(golden("top.coef2.expected.tsv")).unwrap();
    assert_close(
        &parse(&ours_out),
        &parse(&expected),
        "topTable coef2 (golden)",
    );
}

/// Locate an Rscript that has limma installed. Prefers the project's r-bioc
/// conda env (direct binary, no `conda run`), then falls back to PATH.
fn rscript() -> Option<String> {
    let conda = format!(
        "{}/miniconda3/envs/r-bioc/bin/Rscript",
        std::env::var("HOME").unwrap_or_default()
    );
    for cand in [conda.as_str(), "Rscript"] {
        let ok = Command::new(cand)
            .args([
                "-e",
                "if(!requireNamespace('limma',quietly=TRUE)) quit(status=1)",
            ])
            .output()
            .map(|o| o.status.success())
            .unwrap_or(false);
        if ok {
            return Some(cand.to_string());
        }
    }
    None
}

#[test]
fn live_r_coef2() {
    let Some(rs) = rscript() else {
        eprintln!("SKIP live_r_coef2: no Rscript with limma found");
        return;
    };
    let scratch = std::env::temp_dir();
    let r_out = scratch.join(format!("ebayes_r_{}.tsv", std::process::id()));
    let oracle = Command::new(&rs)
        .arg(manifest("tests/ebayes_oracle.R"))
        .arg(golden("expr.tsv"))
        .arg(golden("design.tsv"))
        .arg("2")
        .arg(&r_out)
        .output()
        .unwrap();
    assert!(
        oracle.status.success(),
        "oracle failed: {}",
        String::from_utf8_lossy(&oracle.stderr)
    );
    let ours_out = run_ours(2);
    let r_ref = std::fs::read_to_string(&r_out).unwrap();
    let _ = std::fs::remove_file(&r_out);
    assert_close(&parse(&ours_out), &parse(&r_ref), "topTable coef2 (live R)");
}
