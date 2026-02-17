start <- as.POSIXlt(Sys.time())
context("065de_all.R")

## All of these functions will depend on an expt to play with:
pombe_se <- make_pombe_se(annotation = FALSE)
pombe_subset <- subset_se(
  pombe_se,
  subset = "minute == 0 | minute == 15 | minute == 30") %>%
  set_batches(fact = "replicate")

testing <- limma_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("limma performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt0_vs_mut0"]]
actual <- sum(test[["logFC"]] > 2)
expected <- 9
test_that("limma got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

testing <- noiseq_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("noiseq performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt0_vs_mut0"]]
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

test <- write_noiseq(testing, excel = "test_noiseq_pairwise.xlsx")
test_that("write_noiseq() did something?", {
  expect_true(file.exists("test_noiseq_pairwise.xlsx"))
})

condbatch_keepers <- list("nosig" = c("wt0", "mut0"),
                          "somesig" = c("wt30", "wt0"))

test_cond <- all_pairwise(pombe_subset, model_fstring = "~ 0 + condition",
                          keepers = condbatch_keepers)
actual <- min(test_cond[["comparison"]][["comp"]]["deseq_vs_edger", ])
expected <- 0.96
test_that("all_pairwise() provided results reasonably similar (batch in model)?", {
  expect_gt(actual, expected)
})

test_condbatch <- all_pairwise(pombe_subset, keepers = condbatch_keepers)
actual <- min(test_condbatch[["comparison"]][["comp"]]["deseq_vs_edger", ])
expected <- 0.82
test_that("all_pairwise() provided results reasonably similar (batch in model)?", {
  expect_gt(actual, expected)
})

simplified <- set_conditions(pombe_subset, fact = "strain")
test_sva <- all_pairwise(simplified, model_svs = "svaseq", filter = TRUE,
                         model_fstring = "~ 0 + condition")
actual <- min(test_sva[["comparison"]][["comp"]])
expected <- 0.81
## When testing in 202503, the minimum was actually 0.87
test_that("all_pairwise() provided results reasonably similar? (svaseq in model)", {
  expect_gt(actual, expected)
})

test_cond_combined <- combine_de_tables(test_cond, excel = "testme.xlsx")
## For the life of me I cannot find where this warning is coming from.
## brought out the source of these warnings when I run 'make test'
test_that("combine_de_tables() gave expected tables?", {
  expect_equal(length(test_cond_combined[["data"]]), 2)
})
removed <- file.remove("testme.xlsx")

test_condbatch_combined <- combine_de_tables(test_condbatch, rda = "test_065_combined.rda",
                                             excel = "testme.xlsx")
## For the life of me I cannot find where this warning is coming from.
## brought out the source of these warnings when I run 'make test'
test_that("combine_de_tables() gave expected tables?", {
  expect_equal(length(test_condbatch_combined[["data"]]), 2)
})
removed <- file.remove("testme.xlsx")

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

testing <- compare_de_results(test_condbatch_combined, test_cond_combined)
expected <- 2
actual <- length(testing[["result"]][["limma"]])
test_that("compare_de_results provides some expected output?", {
  expect_equal(expected, actual)
})

expected <- 0.96
actual <- min(testing[["logfc"]])
test_that("compare_de_results provides some expected logfc comparisons?", {
  expect_gt(actual, expected)
})

testing <- correlate_de_tables(test_sva)
actual <- min(testing[["comp"]])
expected <- 0.80
test_that("compare_led_tables provides some expected comparisons?", {
  expect_gt(actual, expected)
})

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

testing <- get_pairwise_gene_abundances(test_sva)
expected <- c(5720, 2)
actual <- dim(testing[["expression_values"]])
test_that("Did get_pairwise_gene_abundances() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

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

testing <- significant_barplots(combined = test_condbatch_combined)
test_that("significant_barplots() gave some plots?", {
  expect_equal(class(testing[["deseq"]])[1], "ggplot2::ggplot")
  expect_equal(class(testing[["limma"]])[1], "ggplot2::ggplot")
  expect_equal(class(testing[["edger"]])[1], "ggplot2::ggplot")
})

testing <- de_venn(test_condbatch_combined[["data"]][[1]])
test_that("de_venn() gave some plots?", {
  expect_equal(class(testing[["up_noweight"]]), "recordedplot")
  expect_equal(class(testing[["down_noweight"]]), "recordedplot")
})

testing <- plot_num_siggenes(test_condbatch_combined[["data"]][[1]])
expected <- "ggplot2::ggplot"
test_that("plot_num_siggenes() gave some plots?", {
  expect_equal(class(testing[["up"]])[1], expected)
  expect_equal(class(testing[["down"]])[1], expected)
  expect_equal(class(testing[["pup"]])[1], expected)
  expect_equal(class(testing[["pdown"]])[1], expected)
})

testing <- extract_abundant_genes(test_sva, excel = NULL)
test_that("extract_abundant_genes() gave some stuff?", {
  expect_equal(100, length(testing[["abundances"]][["deseq"]][["high"]][["mut"]]))
})

testing <- write_de_table(data = test_sva, type = "deseq")
test_that("Did write_de_table() write something?", {
  expect_equal(testing, 1)
})

## Test that we can plot the similarities/difference between
## experiments with Steve Christensen's function.c
compare <- rank_order_scatter(test_condbatch, test_cond)
test_that("Did we compare the two de results with a rank order plot?", {
  expect_equal(class(compare[["plot"]])[1], expected)
  expect_gt(as.numeric(compare[["correlation"]][["estimate"]]), 0.97)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 065de_all.R in ", elapsed,  " seconds.")
