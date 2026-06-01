# rsomics-limma-ebayes

Fit a per-gene linear model to a log-expression matrix and compute
empirical-Bayes moderated t-statistics — a single-binary Rust reimplementation
of limma's `lmFit` + `eBayes` + `topTable`.

Given a log-expression matrix (e.g. microarray intensities or voom `E`) and a
design matrix, it fits least-squares coefficients per gene, shrinks the residual
variances toward a fitted prior, and reports moderated t, moderated p,
BH-adjusted p, and the B-statistic (log-odds of differential expression).

## Usage

```
rsomics-limma-ebayes expr.tsv --design design.tsv [--contrast c.tsv] [--coef N] [-o top.tsv]
```

- `expr.tsv` — header row of sample ids, first column gene ids, log-expression values.
- `--design` — header row of coefficient names, first column sample ids (the model matrix).
- `--contrast` — optional contrast matrix (first column = coefficient names);
  applies `contrasts.fit` before moderation.
- `--coef N` — 1-based coefficient (or contrast) to tabulate; defaults to the
  last column.
- `--proportion` — expected proportion of DE genes for the B-statistic prior (default 0.01).

Output columns match `topTable(sort.by="none")`: `gene logFC AveExpr t P.Value
adj.P.Val B`, in input gene order.

```
rsomics-limma-ebayes E.tsv --design design.tsv --coef 2 -o top.tsv
rsomics-limma-ebayes E.tsv --design design.tsv --contrast c.tsv > top.tsv
```

## Origin

This crate is an independent Rust reimplementation of limma's moderated
t-statistic pipeline (`lmFit` + `eBayes` + `topTable`) based on:

- The published method: Smyth, G.K. (2004), "Linear models and empirical Bayes
  methods for assessing differential expression in microarray experiments",
  Statistical Applications in Genetics and Molecular Biology 3(1):3,
  doi:10.2202/1544-6115.1027 — the moment-based prior estimator (`fitFDist`),
  posterior variance shrinkage, moderated t/p, and the B-statistic.
- Black-box behaviour testing against the limma binary via an R oracle
  (`Rscript` + limma), diffed field-by-field in `tests/compat.rs`.

No source code from limma (GPL) was used as reference during implementation.
Output is value-exact against limma `topTable` (relative deviation < 1e-6 for
coefficients, moderated t/p, adjusted p, and B) across intercept, treatment-
effect, and contrast coefficients.

License: MIT OR Apache-2.0.
Upstream credit: limma (https://bioconductor.org/packages/limma/), GPL (>=2).
