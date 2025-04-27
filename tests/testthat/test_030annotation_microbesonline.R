start <- as.POSIXlt(Sys.time())
context("030annotation_microbesonline.R")

spy_annot <- load_microbesonline_annotations(id = 160490)
actual <- dim(spy_annot)
expected <- c(1871, 18)
test_that("Do we get the expected amount of pyogenes data?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

spy_go <- load_microbesonline_go(id = 160490)
actual <- dim(spy_go)
expected <- c(4594, 2)
test_that("Do we get the expected amount of ecoli GO data?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 030nnotation_microbesonline.R in ", elapsed,  " seconds.")
