library(hpgltools)
pombe_se <- get0("pombe_se")
if (is.null(pombe_se)) {
  source(system.file("examples", "create_pombe.R", package = "hpgltools"))
}

## Let us use limma to compare the pombe wt/mutant samples
limma_de_noint <- limma_pairwise(pombe_se, model_fstring = "~ strain + minute + replicate")
summary(limma_de_noint[["all_tables"]])
head(limma_de_noint[["all_tables"]][[1]])

## We should get identical results via the following model:
limma_zero_int <- limma_pairwise(pombe_se, model_fstring = "~ 0 + strain + minute + replicate")
summary(limma_zero_int[["all_tables"]])
head(limma_zero_int[["all_tables"]][[1]])

cor.test(limma_de_noint[["all_tables"]][[1]][["logFC"]], limma_zero_int[["all_tables"]][[1]][["logFC"]])
