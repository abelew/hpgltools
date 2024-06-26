start <- as.POSIXlt(Sys.time())
context("315norm_batch.R: Are normalizations consistent over time (Batch estimation/correction)?")

## Note to self: Some recent changed to the creation of my expressionsets lead to changes in the order of the resulting data frames.
## This is intended to make it easier for me to keep track of what is happening to the data by forcing it into a consistent order.
## Sadly, this means that most of these tests fail because they assume the previous generic order of genes.  Thus I am now adding
## the gene names to the tests.
## Uses these genes for quick tests
test_genes <- c("FBgn0000014", "FBgn0000008", "FBgn0000017", "FBgn0000018", "FBgn0000024")
## create_expt generates a .Rdata file which may be reread, do so.

load("pasilla_df.rda")
## create_expt generates a .Rdata file which may be reread, do so.
pasilla <- new.env()
load("pasilla.rda", envir = pasilla)
pasilla_expt <- pasilla[["expt"]]

## Test batch
## The following is the previous result, it seems to have changed
## expected <- c(3.333333, 64.500000, 3040.166667, 383.916667, 7.083333)
expected <- c(2.032443, 70.820173, 3357.734214, 379.051162, 6.500123)
names(expected) <- test_genes
pasilla_batch <- sm(normalize_expt(pasilla_expt, batch = "limma"))
actual_df <- exprs(pasilla_batch)
actual <- actual_df[test_genes, c("untreated1")]
test_that("limma batch gives expected values?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

expected <- c(2.03245815, 0.01467331, 0.00000000, 0.11385040, 0.41124803)
names(expected) <- test_genes
pasilla_batch <- sm(normalize_expt(pasilla_expt, batch = "limmaresid"))
actual_df <- exprs(pasilla_batch)
actual <- actual_df[test_genes, c("untreated1")]
test_that("limma-residuals batch gives expected values?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

expected <- c(3.095141, 60.542865, 3032.546240, 355.354483, 6.666536 )
names(expected) <- test_genes
pasilla_batch <- normalize_expt(pasilla_expt, batch = "combatmod", scale = FALSE)
actual_df <- exprs(pasilla_batch)
actual <- actual_df[test_genes, c("untreated1")]
test_that("combatmod from cbcbSEQ batch gives expected values?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

##expected <- c(0.1139956, 7.7406815, 384.8292656, 34.1636051, 0.4937972)
## 20180404: Something here has changed, my initial assumption was that
## yesterday's code which more carefully checks that I have sane options for
## num_surrogates/surrogate_method is incorrect.  I just added some extra text
## to see if that is true, and it appears to be correct.  I am therefore (for
## the moment) assuming that the previous expected values were due to an
## incorrect assumption generated by the previous, less intelligent code.

## Another change: 20181129: I decided to attempt to make the code intelligent
## with respect to the assumptions of the surrogate estimator, for example:
## sva assumes microarray input data (log2ish).  However, when I finished
## collecting the surrogate estimates, I would apply them to the base10 data.
## That seems a bit dumb.  So I added a check for the input state and desired output
## before applying those estimates.  If they do not match, I now have it convert back
## to the appropriate state.
## The good news I suppose, is that the new estimates are almost identical to the old,
## less intelligent way of doing it -- which serves as a reminder that I truly do not
## understand matrix math, I thought they would be a bit more different (not hugely, but
## more than 0.1% - 0.5% different as it appears to be).
## In any event, the new estimates follow, so you can see that they are very similar.
##expected <- c(0.2702588, 6.9005801, 333.1821527, 38.3304667, 0.6674336)
##expected <- c(0.2562567, 6.8468318, 332.5213682, 37.8301610, 0.6589789)
expected <- c(4.047263, 83.959051, 4180.724939, 509.452888, 8.783783)
names(expected) <- test_genes

pasilla_batch <- normalize_expt(pasilla_expt, batch = "sva", convert = "raw")
actual_df <- exprs(pasilla_batch)
actual <- actual_df[test_genes, c("untreated1")]
test_that("sva batch gives expected values?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

new_test_genes <- c("FBgn0000008", "FBgn0000017", "FBgn0000018",
                    "FBgn0000032", "FBgn0000042", "FBgn0000043")
pasilla_batch <- sm(normalize_expt(pasilla_expt, batch = "combat"))
actual_df <- exprs(pasilla_batch)
actual <- actual_df[new_test_genes, c("untreated1")]
## expected <- c(58.61422, 2743.68439, 414.96514, 1082.94955, 74271.77394, 20532.79848)
expected <- c(58.51442, 2738.49321, 415.25486, 1084.05676, 74233.52188, 20476.63402)
names(expected) <- new_test_genes
test_that("combat_noscale gives expected values?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

pasilla_batch <- sm(normalize_expt(pasilla_expt, batch = "combat_scale", filter = TRUE))
## Adding the filter drops FBgn0000014 and 24.
##expected <- c(70.92609, 3436.42054, 411.06522, 1035.16745, 75487.65204, 24292.55511)
expected <- c(70.89932, 3434.87591, 411.16049, 1035.55248, 75502.32575, 24298.59947)
names(expected) <- new_test_genes
actual_df <- exprs(pasilla_batch)
actual <- actual_df[new_test_genes, c("untreated1")]
test_that("combat_scale gives expected values?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

## pasilla_batch <- normalize_expt(pasilla_expt, batch = "combat_noprior_scale") ## takes forever
## The previous result
##expected <- c(0.1139956, 7.7406815, 384.8292656, 34.1636051, 0.4937972)
##expected <- c(4.610009, 82.109047, 4099.039062, 519.407500, 9.116170)
expected <- c(4.047263, 83.959051, 4180.724940, 509.452889, 8.783783)
names(expected) <- test_genes
pasilla_batch <- normalize_expt(pasilla_expt, batch = "svaseq")
actual_df <- exprs(pasilla_batch)
actual <- actual_df[test_genes, c("untreated1")]
test_that("svaseq gives expected values?", {
    expect_equal(expected, actual, tolerance = 0.0001)
})

expected <- c(4, 83, 4091, 496, 9)
names(expected) <- test_genes
pasilla_batch <- normalize_expt(pasilla_expt, batch = "ruvg")
actual_df <- exprs(pasilla_batch)
actual <- actual_df[test_genes, c("untreated1")]
test_that("ruvg gives expected values?", {
    expect_equal(expected, actual, tolerance = 0.5)
})

## The following tests take too much memory

if (FALSE) {
  expected <- c(4.610009, 82.109047, 4099.039071, 519.407501, 9.116170)
  names(expected) <- test_genes
  pasilla_batch <- normalize_expt(pasilla_expt, batch = "varpart")
  actual_df <- exprs(pasilla_batch)
  actual <- actual_df[test_genes, c("untreated1")]
  test_that("variancePartition gives expected values?", {
    expect_equal(expected, actual, tolerance = 0.0001)
  })

  pasilla_batch <- normalize_expt(pasilla_expt, batch = "combat_noprior", filter = TRUE)
}

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 15norm_batch.R in ", elapsed,  " seconds.")
