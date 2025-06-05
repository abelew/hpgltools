start <- as.POSIXlt(Sys.time())
context("072de_dream.R")

## All of these functions will depend on an expt to play with:
pombe_se <- make_pombe_se(annotation = FALSE)
pombe_subset <- subset_se(
  pombe_se,
  subset = "minute == 0 | minute == 15 | minute == 30") %>%
  set_expt_batches(fact = "replicate")

testing <- dream_pairwise(pombe_subset)
actual <- length(testing[["contrasts_performed"]])
expected <- 15
test_that("Dream performed the expected number of contrasts?", {
  expect_equal(expected, actual)
})

test <- testing[["all_tables"]][["wt0_vs_mut0"]]
actual <- sum(test[["logFC"]] > 2)
expected <- 2
test_that("Dream got some expected results (logFC)?", {
  expect_equal(expected, actual)
})

## Add a little fudge-factor to some of these tests.
actual <- sum(as.numeric(test[["P.Value"]]) < 0.1)
expected <- 292
test_that("Dream got some expected results (p)?", {
  expect_equal(expected, actual, tolerance = 3)
})

test <- write_dream(testing, excel = "test_dream_pairwise.xlsx")
test_that("write_dream() did something?", {
  expect_true(file.exists("test_basic_pairwise.xlsx"))
})
