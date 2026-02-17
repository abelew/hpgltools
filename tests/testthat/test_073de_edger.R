start <- as.POSIXlt(Sys.time())
context("073de_edger.R")

pombe_se <- make_pombe_se(annotation = FALSE)
pombe_subset <- subset_se(
  pombe_se,
  subset = "minute == 0 | minute == 15 | minute == 30") %>%
  set_batches(fact = "replicate")

## edger_pairwise()
testing <- edger_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("edgeR performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt0_vs_mut0"]]
actual <- sum(test[["logFC"]] > 2)
expected <- 66
test_that("edgeR got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

actual <- sum(as.numeric(test[["PValue"]]) < 0.1)
expected <- 314
test_that("edgeR got some expected results (adjp)?", {
  expect_equal(expected, actual)
})

test <- write_edger(testing, excel = "test_edger_pairwise.xlsx")
test_that("write_edger() did something?", {
  expect_true(file.exists("test_edger_pairwise.xlsx"))
})
