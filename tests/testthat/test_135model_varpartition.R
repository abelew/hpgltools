start <- as.POSIXlt(Sys.time())
context("135model_varpartition.R")
## 2017-12, exported functions in model_varpartition:
## replot_varpart_percent(), varpart(), varpart_summaries()

pombe_expt <- make_pombe_expt(annotation = TRUE)

pombe_varpart <- simple_varpart(expt = pombe_expt)

## I decided to move away from the mixed models.
##expected <- "(1 | condition) + (1 | batch)"
expected <- "condition + batch"
actual <- as.character(pombe_varpart[["model_used"]])[2]
test_that("Do we get the assumed model?", {
  expect_equal(expected, actual)
})
expected <- c(6266, 3)
actual <- dim(pombe_varpart[["fitted_df"]])
test_that("Did we get an expected table of post-fitting percentages?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

expected <- "gg"
actual <- class(pombe_varpart[["percent_plot"]])[1]
test_that("Does the percent plot get generated?", {
  expect_equal(expected, actual)
})

actual <- class(pombe_varpart[["partition_plot"]])[1]
test_that("Does the partition plot get generated?", {
  expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 135model_varpartition.R in ", elapsed,  " seconds.")
