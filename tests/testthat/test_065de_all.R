start <- as.POSIXlt(Sys.time())
context("065de_all.R")

## All of these functions depend on a SummarizedExperiment to play with:

## Create the requisite data structure and simplify it to just a few pieces.
pombe_se <- make_pombe_se(annotation = FALSE)
pombe_subset <- subset_se(
  pombe_se,
  subset = "minute == 0 | minute == 15 | minute == 30") |>
  set_batches(fact = "replicate")

## Use basic to compare the three conditions remaining in the data.
testing <- basic_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("basic performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})
test_table <- "wt0_vs_mut0"
test <- testing[["all_tables"]][[test_table]]
basic_table <- test
actual <- sum(test[["logFC"]] > 2)
expected <- 1
test_that("basic got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

## Use deseq to compare the three conditions remaining in the data.
testing <- deseq_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("deseq performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})
test_table <- "wt0_vs_mut0"
test <- testing[["all_tables"]][[test_table]]
deseq_table <- test
actual <- sum(test[["logFC"]] > 2)
expected <- 50
test_that("deseq got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

## Use ebseq to compare the three conditions remaining in the data.
testing <- ebseq_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("ebseq performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})
test_table <- "wt0_vs_mut0"
test <- testing[["all_tables"]][[test_table]]
ebseq_table <- test
actual <- sum(test[["logFC"]] > 2)
expected <- 119
test_that("ebseq got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

## Use edger to compare the three conditions remaining in the data.
testing <- edger_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("edger performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})
test_table <- "wt0_vs_mut0"
test <- testing[["all_tables"]][[test_table]]
edger_table <- test
actual <- sum(test[["logFC"]] > 2)
expected <- 66
test_that("edger got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

## Use dream to compare the three conditions remaining in the data.
testing <- dream_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("dream performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})
test_table <- "wt0_vs_mut0"
test <- testing[["all_tables"]][[test_table]]
dream_table <- test
actual <- sum(test[["logFC"]] > 2)
expected <- 2
test_that("dream got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

## Use limma to compare the three conditions remaining in the data.
testing <- limma_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("limma performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})
## We expect a few genes deemed higher in the wt samples than the mutation at time 0.
test_table <- "wt0_vs_mut0"
test <- testing[["all_tables"]][[test_table]]
limma_table <- test
actual <- sum(test[["logFC"]] > 2)
expected <- 8
test_that("limma got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

## Repeat that process using noiseq, keeping in mind that noiseq does not accept models with
## surrogates provided by sva and its ilk.
testing <- noiseq_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("noiseq performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})
test <- testing[["all_tables"]][[test_table]]
noiseq_table <- test
actual <- sum(test[["logFC"]] > 2)
expected <- 1
test_that("noiseq got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

actual <- sum(as.numeric(test[["p"]]) < 0.1)
expected <- 204
test_that("noiseq got some expected results (p)?", {
  expect_equal(expected, actual)
})

## Check that the DESeq2 results and noiseq are reasonably similar.
merged <- merge(limma_table, noiseq_table, by = "row.names")
query <- cor.test(merged[["logFC.x"]], merged[["logFC.y"]])
expected <- 0.9
test_that("Limma and noiseq have similar aresults?", {
  expect_gt(query[["estimate"]], expected)
})

merged <- merge(limma_table, basic_table, by = "row.names")
query <- cor.test(merged[["logFC.x"]], merged[["logFC.y"]])
expected <- 0.85
test_that("Limma and basic have similar aresults?", {
  expect_gt(query[["estimate"]], expected)
})

merged <- merge(limma_table, deseq_table, by = "row.names")
query <- cor.test(merged[["logFC.x"]], merged[["logFC.y"]])
expected <- 0.70
test_that("Limma and deseq have similar aresults?", {
  expect_gt(query[["estimate"]], expected)
})

merged <- merge(limma_table, edger_table, by = "row.names")
query <- cor.test(merged[["logFC.x"]], merged[["logFC.y"]])
expected <- 0.80
test_that("Limma and edger have similar aresults?", {
  expect_gt(query[["estimate"]], expected)
})

merged <- merge(limma_table, dream_table, by = "row.names")
query <- cor.test(merged[["logFC.x"]], merged[["logFC.y"]])
expected <- 0.95
test_that("Limma and dream/varpart have similar aresults?", {
  expect_gt(query[["estimate"]], expected)
})

## Make sure we can write out the noiseq results.
test <- write_noiseq(testing, excel = "test_noiseq_pairwise.xlsx")
test_that("write_noiseq() did something?", {
  expect_true(file.exists("test_noiseq_pairwise.xlsx"))
})

## Now just run a simple model using every method, and perform two
## specific comparisons: wt/mut and wt_30/0
condbatch_keepers <- list("nosig" = c("wt0", "mut0"),
                          "somesig" = c("wt30", "wt0"))

test_cond <- all_pairwise(pombe_subset, model_fstring = "~ 0 + condition",
                          keepers = condbatch_keepers)
## This is a well behaved and simple experiment, we therefore expect the
## various methods to provide extremely similar results.
actual <- min(test_cond[["comparison"]][["comp"]]["deseq_vs_edger", ])
expected <- 0.96
test_that("all_pairwise() provided results reasonably similar (batch in model)?", {
  expect_gt(actual, expected)
})

## In this invocation of all_pairwise(), I did not specify the string which describes
## the model formula (model_fstring), this defaults to "~ 0 + condition + batch".
## The all_pairwise() function assumes the first element (condition) of the formula
## string is the factor of interest and that all following elements (batch in this case)
## are balanced and rank order and therefore quantifiable.
## Note, some tools (dream and maybe noiseq) can support mixed linear models:
## "~ 0 + (1|condition)..." in addition to the various interaction/etc models:
## "~ 0 + condition*batch" or whatever...

## Withe above in mind, this invocation keeps batch in the model.
test_condbatch <- all_pairwise(pombe_subset, keepers = condbatch_keepers)
## We still expect highly similar results between the various methods.
actual <- min(test_condbatch[["comparison"]][["comp"]]["deseq_vs_edger", ])
expected <- 0.82
test_that("all_pairwise() provided results reasonably similar (batch in model)?", {
  expect_gt(actual, expected)
})

## Instead of putting batch in the model, it is possible to use a simplified model
## and then use a tool like sva to append surrogate estimates to it.
## Thus in the following invocation the fstring is just "~ 0 + condition", but
## the model_svs (model surrogate variables) has been set.  Note, that if the data is
## not low-count filtered, this may fail; thus filter was set to TRUE.
simplified <- set_conditions(pombe_subset, fact = "strain")
test_sva <- all_pairwise(simplified, model_svs = "svaseq", filter = TRUE,
                         model_fstring = "~ 0 + condition")
actual <- min(test_sva[["comparison"]][["comp"]])
expected <- 0.81
## When testing in 202503, the minimum was actually 0.87
test_that("all_pairwise() provided results reasonably similar? (svaseq in model)", {
  expect_gt(actual, expected)
})

## Combine the results from the methods and write them to a xlsx file.
## This time, write out just the conditional model.
test_cond_combined <- combine_de_tables(test_cond, excel = "testme.xlsx")
## For the life of me I cannot find where this warning is coming from.
## brought out the source of these warnings when I run 'make test'
test_that("combine_de_tables() gave expected tables?", {
  expect_equal(length(test_cond_combined[["data"]]), 2)
})
removed <- file.remove("testme.xlsx")

## Repeat for the condition+batch model, also write out an addition .rda file which
## may be used by later tests (notably proper)
test_condbatch_combined <- combine_de_tables(test_condbatch, rda = "test_065_combined.rda",
                                             excel = "testme.xlsx")
## For the life of me I cannot find where this warning is coming from.
## brought out the source of these warnings when I run 'make test'
test_that("combine_de_tables() gave expected tables?", {
  expect_equal(length(test_condbatch_combined[["data"]]), 2)
})
removed <- file.remove("testme.xlsx")

## extract_significant_genes() does what it says on the tin, it uses arbitrary definitions
## of 'significant' to pull out genes of interest; the defaults are unsurprisingly a
## |fold change| >= 2.0 and adjusted p-value of <= 0.05.
## Other metrics may include an arbitrary topn by logFC/p, > vs. >=, etc...
cb_de <- test_condbatch
cb_de[["input"]] <- NULL
de_saved <- save(list = "cb_de", file = "test_065_de.rda")
cb_sig <- extract_significant_genes(test_condbatch_combined,
                                    excel = "some_sig.xlsx")
cb_saved <- save(list = "cb_sig", file = "test_065_significant.rda")
test_that("Did we save the result of extract_de_tables?", {
  expect_true(file.exists("test_065_significant.rda"))
})

## Test my plotly writer on the MA plot from deseq.
ggplt_test <- ggplt(test_condbatch_combined[["plots"]][["wt30_vs_wt0"]][["deseq_ma_plots"]])
expected <- "ggplot.html"
actual <- basename(ggplt_test)
test_that("ggplt() returned the filename of a clicky plot?", {
  expect_equal(expected, actual)
})

expected <- 2
actual <- length(test_condbatch_combined[["data"]])
test_that("combine_de_tables() with keepers worked?", {
  expect_equal(expected, actual)
})

## Compare_de_results also does what it says on the tin, it calculates correlation
## coefficients of arbitrary metrics of the results for pairs of tables, defaulting
## to a pearson correlation of the logFC values for two tables.
testing <- compare_de_results(test_condbatch_combined, test_cond_combined)
expected <- 2
actual <- length(testing[["result"]][["limma"]])
test_that("compare_de_results provides some expected output?", {
  expect_equal(expected, actual)
})

## Hey, the limma cor should be > 0.96 for the conditional vs. cond+batch models.
expected <- 0.96
actual <- min(testing[["logfc"]])
test_that("compare_de_results provides some expected logfc comparisons?", {
  expect_gt(actual, expected)
})

## But when adding sva to the mix, that correlation goes down a little.
testing <- correlate_de_tables(test_sva)
actual <- min(testing[["comp"]])
expected <- 0.80
test_that("compare_led_tables provides some expected comparisons?", {
  expect_gt(actual, expected)
})

## get_abundant_genes uses the abundance metric produced by the various methods to
## extract the set of genes which are deemed relatively abundant.  This defaults to using
## the basemean metric from DESeq2 and extracting genes which are > than the IQR;
## but one may choose other metrics of 'abundant' including a static cutoff or
## a multiple of the iqr/etc... (or z-score); it also simultaneously gives the low genes.
testing <- get_abundant_genes(test_sva)
actual <- length(testing[["high"]])
expected <- 2
test_that("Did get_abundant_genes get some stuff?", {
  expect_equal(expected, actual)
})

actual <- names(head(testing[["high"]][["mut"]]))
##expected <- c("SPAC212.09c", "SPAC212.04c", "SPAC977.11",
##              "SPAC977.13c", "SPAC977.15", "SPAC977.16c")
expected <- c("SPRRNA.49", "SPRRNA.01", "SPNCRNA.98",
              "SPRRNA.46", "SPSNRNA.07", "SPBC14F5.04c")
test_that("Did get_abundant_genes get some stuff?", {
  expect_equal(expected, actual)
})

## This takes the result table from the various methods and extracts the expression values
## along with the t statistic, error, and stdev.  This function is crazy old and needs to be
## revisited.
testing <- get_pairwise_gene_abundances(test_sva)
expected <- c(5720, 2)
actual <- dim(testing[["expression_values"]])
test_that("Did get_pairwise_gene_abundances() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## get_sig_genes() does the real work for extract_significant_genes() and works
## at the individual table level.
testing <- get_sig_genes(table = test_sva[["deseq"]][["all_tables"]][[1]])
expected <- c(5, 8)
actual <- dim(testing[["up_genes"]])
test_that("Did get_sig_genes() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

expected <- c(9, 8)
actual <- dim(testing[["down_genes"]])
test_that("Did get_sig_genes() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## The significant barplot is a visualization created (I think) by Laura
## Dillon and re-implemented by me for these results.
testing <- significant_barplots(combined = test_condbatch_combined)
test_that("significant_barplots() gave some plots?", {
  expect_equal(class(testing[["deseq"]])[1], "ggplot2::ggplot")
  expect_equal(class(testing[["limma"]])[1], "ggplot2::ggplot")
  expect_equal(class(testing[["edger"]])[1], "ggplot2::ggplot")
})

## Print out a venn diagram describing the agreement between my three
## favorite methods: deseq/edger/limma; I think this can do arbitrary
## groups of 2 or 3 methods now?
testing <- de_venn(test_condbatch_combined[["data"]][[1]])
test_that("de_venn() gave some plots?", {
  expect_equal(class(testing[["up_noweight"]]), "recordedplot")
  expect_equal(class(testing[["down_noweight"]]), "recordedplot")
})

## These plots are rolling metrics showing how many genes are deemed
## significant as one traverses the possible values for the logFC or
## p-value.  One line for each method employed.
testing <- plot_num_siggenes(test_condbatch_combined[["data"]][[1]])
expected <- "ggplot2::ggplot"
test_that("plot_num_siggenes() gave some plots?", {
  expect_equal(class(testing[["up"]])[1], expected)
  expect_equal(class(testing[["down"]])[1], expected)
  expect_equal(class(testing[["pup"]])[1], expected)
  expect_equal(class(testing[["pdown"]])[1], expected)
})

## Once again, get the abundant genes.
testing <- extract_abundant_genes(test_sva, excel = NULL)
test_that("extract_abundant_genes() gave some stuff?", {
  expect_equal(100, length(testing[["abundances"]][["deseq"]][["high"]][["mut"]]))
})

## Write an individual table instead of a bunch of them.
testing <- write_de_table(data = test_sva, type = "deseq")
test_that("Did write_de_table() write something?", {
  expect_equal(testing, 1)
})

## Test that we can plot the similarities/difference between
## experiments with Steve Christensen's function.
compare <- rank_order_scatter(test_condbatch, test_cond)
test_that("Did we compare the two de results with a rank order plot?", {
  expect_equal(class(compare[["plot"]])[1], expected)
  expect_gt(as.numeric(compare[["correlation"]][["estimate"]]), 0.97)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 065de_all.R in ", elapsed,  " seconds.")
