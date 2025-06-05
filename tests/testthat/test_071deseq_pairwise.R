start <- as.POSIXlt(Sys.time())
context("071de_deseq.R")

## All of these functions will depend on an expt to play with:
pombe_se <- make_pombe_se(annotation = FALSE)
pombe_subset <- subset_se(
  pombe_se,
  subset = "minute == 0 | minute == 15 | minute == 30") %>%
  set_expt_batches(fact = "replicate")

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
