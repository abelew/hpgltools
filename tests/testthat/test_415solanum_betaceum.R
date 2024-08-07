start <- as.POSIXlt(Sys.time())
context("415solanum_betaceum.R: Test some Pseudomonas data.")

sb_annot <- load_trinotate_annotations(
  trinotate = system.file("share/sb/trinotate_head.csv.xz", package = "hpgldata"))

expected <- c(49999, 31)
actual <- dim(sb_annot)
test_that("Do we get the expected trinotate annotation data?", {
  expect_equal(expected, actual)
})

## Try loading some data generated with Sandra and make sure it still works.
untarred <- utils::untar(tarfile = system.file("share/sb/preprocessing.tar.xz",
                                               package = "hpgldata"))

## A little context: Sandra is working with an organism with no genome and no love.
## We therefore passed all the reads to trinity and trinotate and used the resulting,
## entirely bad, putative transcriptome as an input for kallisto and rsem.
## Given that as input, we think we can hunt down some meaningful genes.
## Let us see if that is true.

sb_expt <- create_expt(metadata = "preprocessing/kept_samples.xlsx", gene_info = sb_annot)
actual <- dim(exprs(sb_expt))
expected <- c(234330, 8)
test_that("Do we get an expected matrix from create_expt()?", {
  expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 065de_all.R in ", elapsed,  " seconds.")
