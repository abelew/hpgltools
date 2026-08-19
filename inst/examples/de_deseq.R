library(hpgltools)
pombe_se <- get0("pombe_se")
if (is.null(pombe_se)) {
  source(system.file("examples", "create_pombe.R", package = "hpgltools"))
}

## Let us use DESeq2 to compare the pombe wt/mutant samples
## deseq_lrt performs a likelihood ratio test among a set of interesting factors
## taking into account factors which may interact.  DESeq traditionally puts the
## factor of interest last in the set.
pombe_lrt <- deseq_lrt(pombe_se, interactor_column = "replicate", interest_column = "strain")
head(as.data.frame(pombe_lrt[["deseq_table"]]))
summary(as.data.frame(pombe_lrt[["deseq_table"]]))

## In contrast, one may perform a pairwise differential expression analysis with DESeq2.
## The following invocation compares the strains (wt/mutant)
## with time and replicate in the model.
pombe_deseq_noint <- deseq_pairwise(pombe_se, model_fstring = "~ strain + minute + replicate")
summary(pombe_deseq_noint[["all_tables"]])
head(pombe_deseq_noint[["all_tables"]][[1]])

## We should get identical results via the following model:
pombe_deseq_zero_int <- deseq_pairwise(pombe_se, model_fstring = "~ 0 + strain + minute + replicate")
summary(pombe_deseq_zero_int[["all_tables"]])
head(pombe_deseq_zero_int[["all_tables"]][[1]])
