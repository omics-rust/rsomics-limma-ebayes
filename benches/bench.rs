use std::hint::black_box;
use std::path::PathBuf;

use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_limma_ebayes::{Options, run};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn bench_ebayes(c: &mut Criterion) {
    let expr = fixture("counts.tsv");
    let design = fixture("design.tsv");
    if !expr.exists() {
        return;
    }
    c.bench_function("lmfit_ebayes", |b| {
        b.iter(|| {
            let opts = Options {
                expr: &expr,
                design: &design,
                contrast: None,
                coef: Some(2),
                proportion: 0.01,
            };
            black_box(run(&opts).unwrap());
        })
    });
}

criterion_group!(benches, bench_ebayes);
criterion_main!(benches);
