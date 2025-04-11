start <- as.POSIXlt(Sys.time())
context("020annotation_gff.R")
## I moved get_gff_gene_lengths() to get_genelengths() and made it less stupid.

pa_gff <- system.file("share/paeruginosa_pa14.gff", package = "hpgldata")
pa_fasta <- system.file("share/paeruginosa_pa14.fasta", package = "hpgldata")

pa_irange <- gff2irange(pa_gff)
test_that("Do we get suitable irange data?", {
  expect_equal("GRanges object with 11946 ranges and 11 metadata columns",
               summary(pa_irange))
})

pa_annot <- load_gff_annotations(pa_gff)
test_that("Do we get some gff data for Pseudomonas?", {
  expect_equal(11946, nrow(pa_annot))
  expect_equal(16, ncol(pa_annot))
})

pa_tas <- pattern_count_genome(pa_fasta, gff = pa_gff)
expected <- c(26, 16, 20, 39, 14, 14)
actual <- head(pa_tas[["number"]])
test_that("Do we get sensible numbers of TAs in the pseudomonas genome?", {
  expect_equal(expected, actual)
})

pa_attribs_genes <- sequence_attributes(pa_fasta, gff = pa_gff)
expected <- c(0.62589, 0.37411, 0.4757282, 0.5242718)
actual <- as.numeric(pa_attribs_genes["gene1650835", ])
test_that("Do we get sensible gene attributes by gene?", {
  expect_equal(expected, actual, tolerance = 0.001)
})

pa_attribs_genome <- sequence_attributes(pa_fasta)
expected <- c(0.6629220, 0.3370763, 0.4998674, 0.5001309)
actual <- as.numeric(pa_attribs_genome)
test_that("Do we get sensible gene attributes by genome?", {
  expect_equal(expected, actual, tolerance = 0.001)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 020annotation_gff.R in ", elapsed,  " seconds.")
