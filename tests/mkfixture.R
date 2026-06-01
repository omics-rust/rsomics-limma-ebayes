#!/usr/bin/env Rscript
# Generate a log-expression matrix + design for compat/perf fixtures.
# Usage: mkfixture.R <ngenes> <nsamples_per_group> <expr_out> <design_out> [seed]
args <- commandArgs(trailingOnly = TRUE)
ng <- as.integer(args[1])
nper <- as.integer(args[2])
expr_out <- args[3]
design_out <- args[4]
seed <- if (length(args) >= 5) as.integer(args[5]) else 1L
set.seed(seed)

n <- 2L * nper
# baseline log-expression per gene, with a real signal in some genes
mu <- rnorm(ng, mean = 6, sd = 2)
sd_gene <- sqrt(0.2 + rgamma(ng, shape = 2, rate = 4))
effect <- rep(0, ng)
de <- sample(ng, size = max(1L, round(0.1 * ng)))
effect[de] <- rnorm(length(de), 0, 1.5)

group <- c(rep(0L, nper), rep(1L, nper))
E <- matrix(0, nrow = ng, ncol = n)
for (i in seq_len(ng)) {
  E[i, ] <- mu[i] + effect[i] * group + rnorm(n, 0, sd_gene[i])
}
rownames(E) <- sprintf("g%06d", seq_len(ng))
colnames(E) <- sprintf("s%03d", seq_len(n))

design <- cbind(Intercept = 1, group = group)
rownames(design) <- colnames(E)

write_tsv <- function(mat, path, corner) {
  con <- file(path, "w")
  writeLines(paste(c(corner, colnames(mat)), collapse = "\t"), con)
  g <- function(x) formatC(x, digits = 10, format = "g", flag = "")
  for (i in seq_len(nrow(mat))) {
    writeLines(paste(c(rownames(mat)[i], g(mat[i, ])), collapse = "\t"), con)
  }
  close(con)
}
write_tsv(E, expr_out, "gene")
write_tsv(design, design_out, "sample")
