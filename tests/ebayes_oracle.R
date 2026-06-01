#!/usr/bin/env Rscript
# eBayes oracle: read a log-expression matrix TSV (header = sample ids, col 1 =
# gene ids) and a design matrix TSV (header = coefficient names, col 1 = sample
# ids), run limma lmFit + eBayes, and write topTable(coef, sort.by="none") with
# columns gene,logFC,AveExpr,t,P.Value,adj.P.Val,B.
#
# Usage: ebayes_oracle.R <expr.tsv> <design.tsv> <coef> <out.tsv> [contrast.tsv]
suppressMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
expr_path <- args[1]
design_path <- args[2]
coef <- as.integer(args[3])
out_path <- args[4]
contrast_path <- if (length(args) >= 5) args[5] else NA

E <- as.matrix(read.delim(expr_path, row.names = 1, check.names = FALSE))
design <- as.matrix(read.delim(design_path, row.names = 1, check.names = FALSE))

fit <- lmFit(E, design)
if (!is.na(contrast_path)) {
  cmat <- as.matrix(read.delim(contrast_path, row.names = 1, check.names = FALSE))
  fit <- contrasts.fit(fit, cmat)
}
fit <- eBayes(fit)

tt <- topTable(fit, coef = coef, sort.by = "none", number = Inf)

con <- file(out_path, "w")
writeLines("gene\tlogFC\tAveExpr\tt\tP.Value\tadj.P.Val\tB", con)
g <- function(x) formatC(x, digits = 10, format = "g", flag = "")
for (i in seq_len(nrow(tt))) {
  writeLines(paste(c(
    rownames(tt)[i],
    g(tt$logFC[i]), g(tt$AveExpr[i]), g(tt$t[i]),
    g(tt$P.Value[i]), g(tt$adj.P.Val[i]), g(tt$B[i])
  ), collapse = "\t"), con)
}
close(con)
