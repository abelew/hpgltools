start <- as.POSIXlt(Sys.time())
context("075de_noiseq.R")

## All of these functions will depend on an expt to play with:
pombe_se <- make_pombe_se(annotation = FALSE)
pombe_subset <- subset_se(
  pombe_se,
  subset = "minute == 0 | minute == 15 | minute == 30") %>%
  set_expt_batches(fact = "replicate")

testing <- noiseq_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("noiseq performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt0_vs_mut0"]]
## Random thought, why did I choose 4 fold here?  I might drop this to 2 fold.
actual <- sum(test[["logFC"]] >= 1.0)
expected <- 21
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
