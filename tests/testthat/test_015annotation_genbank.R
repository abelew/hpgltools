start <- as.POSIXlt(Sys.time())
context("015annotation_genbank.R")

## I changed my genbank loaders to use restez and so they are pretty much completely different.
## The good news: they do something again
## I need to decide what I want to test now, though...
testing <- suppressWarnings(load_genbank_annotations())
test_that("Do we get some data?", {
  expect_gt(nrow(testing), 1800)
})

actual_loci <- sort(head(testing[["locus_tag"]]), n = 10)
expected_loci <- c("spyM18_0001", "spyM18_0002", "spyM18_0004",
                   "spyM18_0005", "spyM18_0007", "spyM18_0008")
test_that("Do we get expected Spy IDs?", {
  expect_equal(expected_loci, actual_loci)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 015annotation_genbank.R in ", elapsed,  " seconds.")
