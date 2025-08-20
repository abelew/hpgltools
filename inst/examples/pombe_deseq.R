if (is.null(pombe_se)) {
  source(system.file("examples", "create_se.R", package = "hpgltools"))
}

## Let us just compare two timepoints and ignore the rest of the data for time
pombe_simple <- subset_se(pombe_se,
                          subset = "condition=='wt.120'|condition=='wt.60'")

pombe_deseq <- deseq_pairwise(pombe_simple, filter = TRUE)
