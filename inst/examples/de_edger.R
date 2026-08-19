library(hpgltools)
pombe_se <- get0("pombe_se")
if (is.null(pombe_se)) {
  source(system.file("examples", "create_pombe.R", package = "hpgltools"))
}

## Let us use edger to compare the pombe wt/mutant samples
edger_de_noint <- edger_pairwise(pombe_se, model_fstring = "~ strain + minute + replicate")
summary(edger_de_noint[["all_tables"]])
head(edger_de_noint[["all_tables"]][[1]])

## We should get identical results via the following model:
edger_zero_int <- edger_pairwise(pombe_se, model_fstring = "~ 0 + strain + minute + replicate")
summary(edger_zero_int[["all_tables"]])
head(edger_zero_int[["all_tables"]][[1]])

cor.test(edger_de_noint[["all_tables"]][[1]][["logFC"]], edger_zero_int[["all_tables"]][[1]][["logFC"]])

## Now repeat the same analysis using a simpler model but sva to help with unknown batch effects.
edger_sva <- edger_pairwise(pombe_se, model_fstring = "~ 0 + strain",
                            model_svs = "svaseq", filter = "simple")
summary(edger_sva[["all_tables"]])
head(edger_sva[["all_tables"]][[1]])
