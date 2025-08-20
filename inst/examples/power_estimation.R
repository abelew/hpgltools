if (is.null(pombe_se)) {
  source(system.file("examples", "create_se.R", package = "hpgltools"))
}

if (is.null(pombe_deseq)) {
  source(system.file("examples", "pombe_deseq.R", package = "hpgltools"))
}

pombe_mtrx <- assay(pombe_se)
keepers <- list(
  "wt_120_30" = c("wt.120", "wt.0"),
  "wt_60_30" = c("wt.60", "wt.0"))
pombe_all <- all_pairwise(pombe_se, keepers = keepers)

pombe_proper <- simple_proper(pombe_all, mtrx = pombe_mtrx, de_method = "deseq2",
                              mean_gene_length = 1500)
