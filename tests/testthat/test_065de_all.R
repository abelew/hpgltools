start <- as.POSIXlt(Sys.time())
context("065de_all.R")

## All of these functions will depend on an expt to play with:
pombe_expt <- make_pombe_expt(annotation = FALSE)
pombe_subset <- subset_expt(
  pombe_expt,
  subset = "minute == 0 | minute == 15 | minute == 30")

## Well, in the previous test, we created pombe_expt, so let us use it.
testing <- basic_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("Basic performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt0_vs_mut0"]]
actual <- sum(test[["logFC"]] > 2)
expected <- 1
test_that("Basic got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

## Add a little fudge-factor to some of these tests.
actual <- sum(as.numeric(test[["p"]]) < 0.1)
expected <- 358
test_that("Basic got some expected results (p)?", {
  expect_equal(expected, actual, tolerance = 3)
})

test <- write_basic(testing, excel = "test_basic_pairwise.xlsx")
test_that("write_basic() did something?", {
  expect_true(file.exists("test_basic_pairwise.xlsx"))
})

testing <- deseq_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("DESeq performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt0_vs_mut0"]]
actual <- sum(test[["logFC"]] > 2)
expected <- 50
test_that("DESeq got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

actual <- sum(as.numeric(test[["P.Value"]]) < 0.1)
expected <- 319
test_that("DESeq got some expected results (adjp)?", {
  expect_equal(expected, actual)
})

written_test <- write_deseq(testing, excel = "test_deseq_pairwise.xlsx")
test_that("write_deseq() did something?", {
  expect_true(file.exists("test_deseq_pairwise.xlsx"))
})

## edger_pairwise()
testing <- edger_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("edgeR performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt0_vs_mut0"]]
actual <- sum(test[["logFC"]] > 2)
expected <- 61
test_that("edgeR got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

actual <- sum(as.numeric(test[["PValue"]]) < 0.1)
expected <- 328
test_that("edgeR got some expected results (adjp)?", {
  expect_equal(expected, actual)
})

test <- write_edger(testing, excel = "test_edger_pairwise.xlsx")
test_that("write_edger() did something?", {
  expect_true(file.exists("test_edger_pairwise.xlsx"))
})

## hpgl_voomweighted()
## hpgl_voom()
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

actual <- sum(as.numeric(test[["P.Value"]]) < 0.1)
expected <- 421
test_that("limma got some expected results (adjp)?", {
  expect_equal(expected, actual)
})

test <- write_limma(testing, excel = "test_limma_pairwise.xlsx")
test_that("write_limma() did something?", {
  expect_true(file.exists("test_limma_pairwise.xlsx"))
})

test_condbatch <- all_pairwise(pombe_subset)
actual <- min(test_condbatch[["comparison"]][["comp"]]["deseq_vs_edger", ])
expected <- 0.80
test_that("all_pairwise() provided results reasonably similar (batch in model)?", {
  expect_gt(actual, expected)
})

test_cond <- all_pairwise(pombe_subset, model_batch = FALSE)
actual <- min(test_cond[["comparison"]][["comp"]]["deseq_vs_edger", ])
expected <- 0.76
test_that("all_pairwise() provided results reasonably similar (no batch in model)?", {
  expect_gt(actual, expected)
})

test_sva <- all_pairwise(pombe_subset, model_batch = "svaseq", filter = TRUE)
actual <- min(test_sva[["comparison"]][["comp"]])
expected <- 0.63
test_that("all_pairwise() provided results reasonably similar? (svaseq in model)", {
  expect_gt(actual, expected)
})

cond_model <- choose_model(pombe_subset, model_batch = FALSE)
expected <- "~ 0 + condition"
actual <- cond_model[["chosen_string"]]
test_that("choose_model provides expected models?", {
  expect_equal(expected, actual)
})

## I think forcing the expt to keep conditions as specifically ordered factors
## has caused some oddities in how the downstream model is getting made.
## I think I should therefore ensure fully that the conditions of the model
## match the conditions in the original design.
model_df <- as.data.frame(cond_model[["chosen_model"]])
test_df <- data.frame(row.names = names(pombe_subset[["conditions"]]))
for (cond in pombe_subset[["conditions"]]) {
  test_df[[cond]] <- 0
}
for (c in 1:length(pombe_subset[["conditions"]])) {
  name <- names(pombe_subset[["conditions"]])[c]
  value <- as.character(pombe_subset[["conditions"]][c])
  test_df[name, value] <- 1
}
test_that("choose_model provides a model which matches the design?", {
  expect_equal(model_df, test_df)
})

condbatch_model <- choose_model(pombe_subset)
expected <- "~ 0 + condition + batch"
actual <- condbatch_model[["chosen_string"]]
test_that("choose_model provides expected models?", {
  expect_equal(expected, actual)
})

## choose_dataset()
testing <- choose_dataset(pombe_subset, choose_for = "limma")
expected <- c("libsize", "conditions", "batches", "data")
actual <- names(testing)
test_that("choose_dataset provides some expected output?", {
  expect_equal(expected, actual)
})

testing <- choose_dataset(pombe_subset, choose_for = "deseq")
expected <- c("libsize", "conditions", "batches", "data")
actual <- names(testing)
test_that("choose_dataset provides some expected output?", {
  expect_equal(expected, actual)
})

testing <- choose_dataset(pombe_subset, choose_for = "edger")
expected <- c("libsize", "conditions", "batches", "data")
actual <- names(testing)
test_that("choose_dataset provides some expected output?", {
  expect_equal(expected, actual)
})

## we did test_condbatch, test_cond, test_sva
keepers <- list("nosig" = c("wt0", "mut0"),
                "somesig" = c("wt30", "wt0"))
test_condbatch_combined <- combine_de_tables(test_sva,
                                             keepers = keepers,
                                             excel = "testme.xlsx")
## For the life of me I cannot find where this warning is coming from.
## brought out the source of these warnings when I run 'make test'
## 24
test_that("combine_de_tables() gave expected tables?", {
  expect_equal(length(test_condbatch_combined[["data"]]), 2)
})
removed <- file.remove("testme.xlsx")

small_combined <- combine_de_tables(test_condbatch, keepers = keepers)
saved <- save(list = c("small_combined"), file = "065_small_combined.rda")
test_that("Did we save the result of combine_de_tables?", {
  expect_true(file.exists("065_small_combined.rda"))
})

cb_sig <- extract_significant_genes(small_combined,
                                    excel = "some_sig.xlsx")
cb_saved <- save(list = "cb_sig", file = "test_065_significant.rda")
test_that("Did we save the result of extract_de_tables?", {
  expect_true(file.exists("test_065_significant.rda"))
})

## Test my plotly writer on the MA plot from deseq.
ggplt_test <- ggplt(small_combined[["plots"]][["first"]][["deseq_ma_plots"]][["plot"]])
expected <- "ggplot.html"
actual <- basename(ggplt_test)
test_that("ggplt() returned the filename of a clicky plot?", {
  expect_equal(expected, actual)
})

expected <- 2
actual <- length(small_combined[["data"]])
test_that("combine_de_tables() with keepers worked?", {
  expect_equal(expected, actual)
})

## Same query, condition in model
test_cond_combined <- combine_de_tables(test_cond,
                                        keepers = keepers)
test_that("combine_de_tables() gave expected tables?", {
  expect_equal(length(test_cond_combined[["data"]]), 2)
})

testing <- compare_de_results(test_condbatch_combined, test_cond_combined)
expected <- 18
actual <- length(unlist(testing[["result"]]))
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
expected <- 0.70
test_that("compare_led_tables provides some expected comparisons?", {
  expect_gt(actual, expected)
})

## Strange, I got a failure here when running make test
## but running manually everything seems to be working fine...
## 15 compare_significant_contrasts()

## Note that I limited the combine_de_tables() to grabbing just 1 or two tables,
## so a bunch of these tests stopped working.  They should perhaps be redone.

## do_pairwise()
## Saving this so we can use it for ontology searches later.
save(list = "test_condbatch_combined", file = "test_065_combined.rda", compress = TRUE)
save(list = "cb_sig", file = "test_065_significant.rda", compress = TRUE)
## This is done by a bunch of other functions, I am not testing it.

testing <- get_abundant_genes(test_sva)
actual <- length(testing[["high"]])
expected <- 6
test_that("Did get_abundant_genes get some stuff?", {
  expect_equal(expected, actual)
})

actual <- names(head(testing[["high"]][["mut0"]]))
##expected <- c("SPAC212.09c", "SPAC212.04c", "SPAC977.11",
##              "SPAC977.13c", "SPAC977.15", "SPAC977.16c")
expected <- c("SPRRNA.49", "SPRRNA.01", "SPNCRNA.98",
              "SPRRNA.46", "SPSNRNA.07", "SPBC14F5.04c")
test_that("Did get_abundant_genes get some stuff?", {
  expect_equal(expected, actual)
})

testing <- get_pairwise_gene_abundances(test_sva)
expected <- c(5720, 6)
actual <- dim(testing[["expression_values"]])
test_that("Did get_pairwise_gene_abundances() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

testing <- get_sig_genes(table = test_sva[["deseq"]][["all_tables"]][[1]])
expected <- c(199, 6)
actual <- dim(testing[["up_genes"]])
test_that("Did get_sig_genes() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_lt(expected[2], actual[2])
})

expected <- c(183, 6)
actual <- dim(testing[["down_genes"]])
test_that("Did get_sig_genes() get some stuff?", {
  expect_equal(expected[1], actual[1])
  expect_lt(expected[2], actual[2])
})

pombe_model <- choose_model(pombe_subset)
testing <- make_pairwise_contrasts(model = pombe_model[["chosen_model"]],
                                   conditions = pombe_subset[["conditions"]])
actual <- length(names(testing[["all_pairwise"]]))
expected <- 15
test_that("Did make_pairwise_contrasts() get some stuff?", {
  expect_equal(expected, actual)
})

## If we add back some experimental factors, we should get bigger
## models/contrast lists.
pombe_model <- choose_model(pombe_expt)
testing <- make_pairwise_contrasts(model = pombe_model[["chosen_model"]],
                                   conditions = pombe_expt[["conditions"]])
actual <- length(names(testing[["all_pairwise"]]))
expected <- 66
test_that("Did make_pairwise_contrasts() get some stuff?", {
  expect_equal(expected, actual)
})

testing <- significant_barplots(combined = test_condbatch_combined)
test_that("significant_barplots() gave some plots?", {
  expect_equal(class(testing[["deseq"]]), c("gg", "ggplot"))
  expect_equal(class(testing[["limma"]]), c("gg", "ggplot"))
  expect_equal(class(testing[["edger"]]), c("gg", "ggplot"))
})

testing <- de_venn(test_condbatch_combined[["data"]][[1]])
test_that("de_venn() gave some plots?", {
  expect_equal(class(testing[["up_noweight"]]), "recordedplot")
  expect_equal(class(testing[["down_noweight"]]), "recordedplot")
})

testing <- plot_num_siggenes(test_condbatch_combined[["data"]][[1]])
expected <- c("gg", "ggplot")
test_that("plot_num_siggenes() gave some plots?", {
  expect_equal(class(testing[["up"]]), expected)
  expect_equal(class(testing[["down"]]), expected)
  expect_equal(class(testing[["pup"]]), expected)
  expect_equal(class(testing[["pdown"]]), expected)
})

testing <- extract_abundant_genes(test_sva, excel = NULL)
test_that("extract_abundant_genes() gave some stuff?", {
  expect_equal(100, length(testing[["abundances"]][["deseq"]][["high"]][["mut0"]]))
})

testing <- write_de_table(data = test_sva, type = "deseq")
test_that("Did write_de_table() write something?", {
  expect_equal(testing, 1)
})

## Test that we can plot the similarities/difference between
## experiments with Steve Christensen's function.
compare <- rank_order_scatter(test_condbatch, test_cond)
test_that("Did we compare the two de results with a rank order plot?", {
  expect_equal(class(compare[["plot"]])[1], "gg")
  expect_gt(as.numeric(compare[["correlation"]][["estimate"]]), 0.99)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 065de_all.R in ", elapsed,  " seconds.")
