start <- as.POSIXlt(Sys.time())
context("015annotation_genbank.R")

## I changed my genbank loaders to use restez and so they are pretty much completely different.
## The good news: they do something again
## I need to decide what I want to test now, though...
testing <- suppressWarnings(load_genbank_annotations())
test_that("Do we get some data?", {
  expect_gt(nrow(testing), 1800)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 015annotation_genbank.R in ", elapsed,  " seconds.")
