start <- as.POSIXlt(Sys.time())
context("005annotation_biomart.R")

## This pretty consistently times out when run non-interactively, which is weird because
## I use the same query pretty often in my work and never have a problem.
## This I think I will change the query to a simpler, yeast dataset.
##testing <- load_biomart_annotations(species = "hsapiens",
##                                    year = "2022", month = "apr")
sad <- FALSE
testing <- try(load_biomart_annotations(
  species = "spombe", host = "nov2020-fungi.ensembl.org",
  trymart = "fungi_mart", trydataset = "spombe_eg_gene",
  gene_requests = c("pombase_transcript", "ensembl_gene_id", "ensembl_transcript_id",
                    "hgnc_symbol", "description", "gene_biotype"), overwrite = TRUE))
if ("try-error" %in% class(testing)) {
  sad <- TRUE
}

testthat::skip_if(sad)
annotations <- testing[["annotation"]]
gene_ids <- head(sort(annotations[["ensembl_transcript_id"]]))
expected <- c(
  "ENSRNA049622790-T1", "ENSRNA049622863-T1", "ENSRNA049622903-T1", "ENSRNA049622995-T1",
  "ENSRNA049623050-T1", "ENSRNA049623137-T1")
test_that("Do we get expected gene IDs?", {
  expect_equal(expected, gene_ids)
})

expected <- 7260
actual <- nrow(annotations)
test_that("Do we receive expected output from load_biomart_annotations()?", {
  expect_gt(actual, expected)
})

testing <- load_biomart_go(
  species = "spombe", host = "nov2020-fungi.ensembl.org",
  trymart = "fungi_mart", trydataset = "spombe_eg_gene", overwrite = TRUE)
data <- testing[["go"]]
expected <- 59800
actual <- nrow(data)
test_that("Do we receive expected output from load_biomart_go()?", {
  expect_gt(actual, expected)
})

## It appears the feb 2021 archives are not responding right now.
## I think I should spend some time and make a mirror
testing <- load_biomart_orthologs(gene_ids = gene_ids, first_species = "hsapiens",
                                  second_species = "mmusculus", year = 2021,
                                  month = "feb")
data <- testing[["all_linked_genes"]]
actual <- nrow(data)
expected <- 23000
test_that("Do we get expected orthologs from load_biomart_orthologs()?", {
  expect_gt(actual, expected)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 005annotation_biomart.R in ", elapsed,  " seconds.")
