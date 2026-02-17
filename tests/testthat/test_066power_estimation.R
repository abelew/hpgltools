start <- as.POSIXlt(Sys.time())
context("066power_estimation.R")

pombe_se <- make_pombe_se()
apr <- new.env()
tt <- load(file = "test_065_de.rda", envir = apr)
apr <- apr[["cb_de"]]
small_combined <- new.env()
tt <- load(file = "test_065_combined.rda", envir = small_combined)
small_combined <- small_combined[["test_065_combined"]]
test_proper <- simple_proper(small_combined, apr = apr, mtrx = pombe_se,
                             reps = c(3, 5), nsims = 10)
expected <- 6
actual <- nrow(test_proper[[1]][["power_table"]])
test_that("Minimal check for proper functionality:", {
  expect_equal(expected, actual)
})

expected <- "recordedplot"
actual <- class(test_proper[[1]][["powerfd_plot"]])
test_that("Minimal check for proper plotting:", {
  expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 066power_estimation.R in ", elapsed,  " seconds.")
