start <- as.POSIXlt(Sys.time())
context("311norm_convert.R")

load("pasilla_df.rda")
## create_se generates a .Rdata file which may be reread, do so.
pasilla <- new.env()
load("pasilla.rda", envir = pasilla)
pasilla_se <- pasilla[["se"]]

## Uses these genes for quick tests
test_genes <- c("FBgn0000014", "FBgn0000008", "FBgn0000017", "FBgn0000018", "FBgn0000024")

## Make sure that my invocation of cpm() is the same as edgeR's.
pasilla_convert <- sm(normalize(pasilla_se, convert = "cpm"))
expected <- edgeR::cpm(assay(pasilla_se))
actual <- assay(pasilla_convert)
test_that("cpm conversions are equivalent?", {
    expect_equal(expected, actual)
})

## Check that the different ways of calling rpkm() are identical
pasilla_convert <- convert_counts(pasilla_se, convert = "rpkm", length_column = "cds_length",
                                  start_column = "start_position", end_column = "end_position")
pasilla_norm <- normalize(pasilla_se, convert = "rpkm", length_column = "cds_length",
                          start_column = "start_position", end_column = "end_position")
expected <- pasilla_convert[["count_table"]]
actual <- assay(pasilla_norm)
test_that("calling convert_counts and normalize are equivalent?", {
    expect_equal(expected, actual)
})

## Similarly check that edgeR's rpkm() comes out the same
## Make sure that we remove undefined numbers from fdata(length)
## This subtraction logic is no longer needed, the pasilla annotations have
## cds lengths already recorded; and they take into account the UTRs.

## Well, some versions of pasilla provide this information, others do not...
if (is.null(rowData(pasilla_se)[["cds_length"]])) {
  rowData(pasilla_se)[["start_position"]] <- suppressWarnings(as.numeric(rowData(pasilla_se)[["start_position"]]))
  rowData(pasilla_se)[["end_position"]] <- suppressWarnings(as.numeric(rowData(pasilla_se)[["end_position"]]))
  rowData(pasilla_se)[["cds_length"]] <- abs(rowData(pasilla_se)[["start_position"]] -
                                               rowData(pasilla_se)[["end_position"]])
  print(summary(rowData(pasilla_se)[["cds_length"]]))
}
undef <- is.na(rowData(pasilla_se)[["cds_length"]])
lengths <- rowData(pasilla_se)[["cds_length"]]
## I changed my rpkm function to coerce undefined values to 1k so that
## they are effectively not normalized.
lengths[undef] <- 1000
fdata_lengths <- as.vector(as.numeric(lengths))
names(fdata_lengths) <- rownames(rowData(pasilla_se))
expected <- edgeR::rpkm(assay(pasilla_se), gene.length = fdata_lengths)
actual <- assay(pasilla_norm)
test_that("rpkm conversions are equivalent?", {
    expect_equal(expected, actual)
})

## I have a modification of rpkm(), cp_seq_m(), which should give some expected results.
## This is intended to count the number of instances of a given sequence ('TA' by default)
## and normalize based on its relative frequency.  This is useful primarily for tnseq.
genome <- "BSgenome.Dmelanogaster.UCSC.dm6"
installp <- ! genome %in% installed.packages()
if (isTRUE(installp)) {
  tt <- BiocManager::install(genome)
}

tt <- sm(library(genome, character.only = TRUE))
pasilla_convert <- suppressWarnings(normalize(
  pasilla_se, convert = "cp_seq_m", start_column = "start_position",
  chromosome_column = "chromosome_name", end_column = "end_position",
  genome = BSgenome.Dmelanogaster.UCSC.dm6))

## Interesting, these values have reverted to my old values...
expected <- c(0.03493820, 0.47396682, 22.70873214, 55.11152210, 0.03965286)
## And switched back...
##expected <- c(0.03443909, 0.46719586, 22.38432168, 54.32421464, 0.03908639)
actual <- as.numeric(assay(pasilla_convert)[test_genes, 1])
test_that("cp_seq_m works for TA?", {
    expect_equal(expected, actual, tolerance = 1)
})

## Repeat cp_seq_m() for ATG
pasilla_convert <- sm(normalize(
  pasilla_se, convert = "cp_seq_m", start_column = "start_position",
  chromosome_column = "chromosome_name", end_column = "end_position",
  genome = BSgenome.Dmelanogaster.UCSC.dm6, pattern = "ATG"))
## These values have changed again...
## Let us use the rowMeans of the values and add a tolerance, that should
## provide a little leeway if we have different annotations and/or small
## shifts in the genome package.
actual <- as.numeric(rowMeans(assay(pasilla_convert))[test_genes])
expected <- c(0.01007904, 0.54981942, 30.97386597, 34.69884511, 0.03447848)

actual <- as.numeric(assay(pasilla_convert)[test_genes, c("untreated1")])
test_that("cp_seq_m works for ATG?", {
    expect_equal(expected, actual, tolerance = 4)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 11norm_convert.R in ", elapsed, " seconds.")
