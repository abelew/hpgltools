start <- as.POSIXlt(Sys.time())
context("145nmer.R")

tt <- sm(library("BSgenome.Dmelanogaster.UCSC.dm6"))
num_starts <- count_nmer(BSgenome.Dmelanogaster.UCSC.dm6)
expected <- 2440952
actual <- sum(num_starts)
test_that("Can we count ATGs via count_nmer()?", {
    expect_equal(expected, actual, tolerance = 10)
})

load("pasilla_df.rda")
## create_se generates a .Rdata file which may be reread, do so.
pasilla <- new.env()
load("pasilla.rda", envir = pasilla)
pasilla_se <- pasilla[["se"]]

## There is a discrepency in how pasilla names chromosomes and how
## the bsgenome does...
annot_df <- rowData(pasilla_se)
annot_df[["short_chromosome_name"]] <- gsub(x = annot_df[["chromosome_name"]],
                                            pattern = "^chr", replacement = "")

## Note, sometimes ensembl works, sometimes it doesn't; depending on
## that we will need to use different columns from the annotation
## data.
print(colnames(annot_df))
if (is.null(annot_df[["TXSTRAND"]])) {
  message("The annotations came from an orgdb, ensembl must have been down.")
  annot_df[["chromosome_name"]] <- paste0("chr", annot_df[["chromosome_name"]])
  padding <- gather_utrs_padding(
    BSgenome.Dmelanogaster.UCSC.dm6, annot_df = annot_df, name_column = "ensembl_gene_id",
    start_column = "start_position", chr_column = "chromosome_name", end_column = "end_position",
    strand_column = "strand", gene_type = "protein_coding", type_column = "gene_biotype")
  expected <- 4506
  actual <- nrow(padding[["fiveprime_plus_table"]])
  test_that("Collect padding of CDS from the melanogaster genome?", {
    expect_equal(expected, actual, tolerance = 2)
  })
} else {
  message("The annotations came from ensembl.")
  padding <- gather_utrs_padding(
    BSgenome.Dmelanogaster.UCSC.dm6, annot_df = annot_df, name_column = "ensembl",
    start_column = "start_position", chr_column = "chromosome_name", end_column = "end_position",
    strand_column = "TXSTRAND", gene_type = "protein-coding", type_column = "genetype")
  expected <- 4508
  actual <- nrow(padding[["fiveprime_plus_table"]])
  test_that("Collect padding of CDS from the melanogaster genome?", {
    expect_equal(expected, actual, tolerance = 2)
  })
}

attributes <- sequence_attributes(BSgenome.Dmelanogaster.UCSC.dm6)
actual <- attributes["chr2L", "gc"]
expected <- 0.4178
test_that("Collect padding of CDS from the melanogaster genome?", {
    expect_equal(expected, actual, tolerance = 0.1)
})

fasta_file <- system.file("share", "paeruginosa_pa14.fasta", package = "hpgldata")
bs <- make_bsgenome_from_fasta(
  fasta = fasta_file, pkgname = "BSgenome.pseudomonas.aeruginosa.pa14",
  organism = "Pseudomonas aeruginosa strain PA14",
  title = "Test pseudomonas", common_name = "PA14", provider = "NCBI", installp = FALSE)
test_that("Able to create a BSgenome from fasta file?", {
  expect_equal(bs[["contigs"]], 1)
})

unlink("build", recursive = TRUE)

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 145nmer.R in ", elapsed,  " seconds.")
