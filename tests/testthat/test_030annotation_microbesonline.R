start <- as.POSIXlt(Sys.time())
context("030annotation_microbesonline.R")

## Load Streptococcus pyogenes annotations by the genome ID number
spy_annot <- load_microbesonline_annotations(id = 160490)
actual <- dim(spy_annot)
expected <- c(1871, 18)
test_that("Do we get the expected amount of pyogenes data?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## Use the same information to extract the GO data.
spy_go <- load_microbesonline_go(id = 160490)
actual <- dim(spy_go)
expected <- c(4594, 2)
test_that("Do we get the expected amount of ecoli GO data?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

## Use a portion of the species name to extract the ID number, then load the annotations.
## This time we will use Streptococcus agalactiae A909
a909_annot <- load_microbesonline_annotations(species = "A909")
actual <- dim(a909_annot)
expected <- c(2229, 18)
test_that("Do we get the expected amount of agalactiae data?", {
  expect_equal(expected[1], actual[1])
  expect_equal(expected[2], actual[2])
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 030nnotation_microbesonline.R in ", elapsed,  " seconds.")
