start <- as.POSIXlt(Sys.time())
context("070de_basic.R")

## All of these functions will depend on an expt to play with:
pombe_se <- make_pombe_se(annotation = FALSE)
pombe_subset <- subset_se(
  pombe_se,
  subset = "minute == 0 | minute == 15 | minute == 30") %>%
  set_expt_batches(fact = "replicate")

## Well, in the previous test, we created pombe_se, so let us use it.
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
