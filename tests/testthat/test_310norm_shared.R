start <- as.POSIXlt(Sys.time())
context("310norm_shared.R")

## Note to self: Some recent changed to the creation of my expressionsets lead to changes in the order of the resulting data frames.
## This is intended to make it easier for me to keep track of what is happening to the data by forcing it into a consistent order.
## Sadly, this means that most of these tests fail because they assume the previous generic order of genes.  Thus I am now adding
## the gene names to the tests.

load("pasilla_df.rda")
## create_se generates a .Rdata file which may be reread, do so.
pasilla <- new.env()
load("pasilla.rda", envir = pasilla)
pasilla_se <- pasilla[["se"]]

test_genes <- c("FBgn0000014","FBgn0000008","FBgn0000017","FBgn0000018", "FBgn0000024")

## First make sure the pasilla_se still has the stuff we expect
expected <- "This is a summarized experiment."
actual <- metadata(pasilla_se)[["title"]]
test_that("Pasilla title?", {
    expect_equal(expected, actual)
})

## Ensure that the beginning count table library sizes are identical.
expected <- colSums(counts)
actual <- libsize(pasilla_se)
names(expected) <- c("untreated1", "untreated2", "untreated3", "untreated4",
                     "treated1", "treated2", "treated3")
test_that("Pasilla libsize?", {
    expect_equal(expected, actual)
})

## Check a few arbitrary counts to make sure they are maintained.
expected <- counts[test_genes, "untreated1"]
testing_counts <- assay(pasilla_se)
actual <- as.numeric(testing_counts[test_genes, "untreated1"])
test_that("Pasilla count tables? (untreated1)", {
    expect_equal(expected, actual)
})

## Check that all samples agree for 1 gene.
test_gene <- "FBgn0062565"
expected <- as.numeric(counts[test_gene, ])
actual <- as.numeric(assay(pasilla_se)[test_gene, ])
expected <- c(4, 7, 3, 3, 9, 10, 9)
test_that("Pasilla count tables? (gene FBgn0063565)", {
    expect_equal(expected, actual)
})

## Ensure that normalize does not mess up the data when called without arguments (this wasn't true once)
unmolested <- sm(normalize(pasilla_se))
expected <- as.matrix(assay(pasilla_se))  ## I accidently changed this to potentially return a data.frame
actual <- assay(unmolested)
test_that("Pasilla (un)normalized counts?", {
    expect_equal(expected, actual)
})

## Make sure that the colData information is maintained through normalization
expected <- colData(pasilla_se)
actual <- colData(unmolested)
test_that("Pasilla (un)normalized pdata?", {
    expect_equal(expected, actual)
})

## Also ensure that the library sizes (which are very important for limma) are not messed up.
expected <- pasilla_se[["libsize"]]
actual <- unmolested[["libsize"]]
test_that("Pasilla (un)normalized libsize?", {
    expect_equal(expected, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 10norm_shared.R in ", elapsed, " seconds.")
