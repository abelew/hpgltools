start <- as.POSIXlt(Sys.time())
context("302graph_metrics.R: Is it possible to graph the various metrics with hpgltools?")

pasilla <- new.env()
load("pasilla.rda", envir = pasilla)
pasilla_se <- pasilla[["se"]]
## Uses these genes for quick tests
test_genes <- c("FBgn0000014", "FBgn0000008", "FBgn0000017", "FBgn0000018", "FBgn0000024")

## What graphs can we make!?
libsize_plot <- plot_libsize(pasilla_se)
actual <- libsize_plot[["table"]][["sum"]]
expected <- c(13971670, 21909886, 8357876, 9840745, 18668667, 9571213, 10343219)
test_that("The libsize plot is as expected?", {
    expect_equal(expected, actual)
})

nonzero_plot <- plot_nonzero(pasilla_se)
actual <- nonzero_plot[["table"]][["nonzero_genes"]]
expected <- c(9863, 10074, 9730, 9786, 10087, 9798, 9797)
test_that("The non-zero genes is as expected?", {
    expect_equal(expected, actual)
})

## These tests have also been affected by the changed order of expressionsets.
density <- sm(plot_density(pasilla_se))
density_plot <- density[["plot"]]
density_table <- density[["table"]]
expected <- c(92, 5, 4664, 583, 10, 1446)
actual <- head(density_table[["counts"]])
test_that("Density plot data is as expected?", {
    expect_equal(expected, actual)
})

hist_plot <- sm(plot_histogram(data.frame(exprs(pasilla_se))))
actual <- head(hist_plot[["data"]][["values"]])
## The values of expected have not changed
test_that("Histogram data is as expected?", {
    expect_equal(expected, actual)
})

box_plot <- sm(plot_boxplot(pasilla_se))
actual <- head(box_plot[["plot"]][["data"]][["reads"]])
test_that("Box plot data is as expected?", {
    expect_equal(expected, actual, tolerance = 1)
})

## Ahh yes I changed the cbcb_filter options to match those from the cbcbSEQ vignette.
## Note that the filtering has changed slightly, and this affects the results.
norm <- sm(normalize(pasilla_se, transform = "log2", convert = "cbcbcpm",
                     norm = "quant", filter = TRUE))
expected <- "recordedplot"  ## for all the heatmaps

corheat_plot <- plot_corheat(norm)
actual <- class(corheat_plot[["plot"]])
test_that("corheat is a recorded plot?", {
    expect_equal(expected, actual)
})

disheat_plot <- plot_disheat(norm)
actual <- class(disheat_plot[["plot"]])
test_that("disheat is a recorded plot?", {
    expect_equal(expected, actual)
})

sampleheat_plot <- plot_sample_heatmap(norm)
actual <- class(sampleheat_plot)
test_that("sampleheat is a recorded plot?", {
    expect_equal(expected, actual)
})

smc_plot <- sm(plot_sm(norm, method = "pearson"))
actual <- head(smc_plot[["plot"]][["data"]][["sm"]])
expected <- c(0.9759981, 0.9824316, 0.9759981, 0.9821373, 0.9784851, 0.9786376)
test_that("Is the normalized smc data expected?", {
    expect_equal(expected, actual, tolerance = 0.004)
})

smd_plot <- sm(plot_sm(norm, method = "euclidean"))
actual <- head(smd_plot[["plot"]][["data"]][["sm"]])
## 201812 Changed due to peculiarities in normalization methods.
## 201907 I changed the normalization back, so the values need to return.
## expected <- c(42.43941, 36.43116, 42.43941, 36.60569, 40.01228, 40.04465)
expected <- c(38.65624, 33.27842, 38.65624, 32.89692, 38.66213, 37.13604)
## 10
test_that("Is the normalized smd data expected?", {
    expect_equal(expected, actual, tolerance = 0.01)
})

pca_stuff <- plot_pca(norm)
pca_plot <- pca_stuff[["plot"]]
pca_pca <- head(pca_stuff[["pca"]])

actual <- pca_plot[["data"]][["PC1"]]
##expected <- c(-0.3595704, -0.4034877, -0.2721046, -0.2432659, 0.2850077, 0.4995150, 0.4939059)
expected <- c(-0.3497219, -0.3822099, -0.2997402, -0.2642127,  0.3419079,  0.4802172,  0.4737595)
test_that("Is the pca data as expected for PC1?", {
    expect_equal(expected, actual, tolerance = 0.01)
})

actual <- as.numeric(head(pca_stuff[["result"]][["v"]][, 1]))
##expected <- c(-0.3595704, -0.4034877, -0.2721046, -0.2432659, 0.2850077, 0.4995150)
expected <- c(-0.3497219, -0.3822099, -0.2997402, -0.2642127, 0.3419079, 0.4802172)
test_that("Is the SVD 'v' element expected?", {
    expect_equal(expected, actual, tolerance = 0.01)
})

actual <- pca_stuff[["residual_df"]][[1]]
##expected <- c(42.44, 31.26, 13.09, 5.84, 4.13, 3.24)
expected <- c(43.78, 28.86, 14.30, 5.97, 3.99, 3.11)
test_that("Is the pca residual table as expected?", {
    expect_equal(expected, actual, tolerance = 0.01)
})

actual <- pca_stuff[["prop_var"]]
##expected <- c(42.44, 31.26, 13.09, 5.84, 4.13, 3.24)
expected <- c(43.78, 28.86, 14.30, 5.97, 3.99, 3.11)
test_that("Is the variance list as expected?", {
    expect_equal(expected, actual, tolerance = 0.01)
})

actual <- pca_stuff[["table"]][["PC2"]]
##expected <- c(0.3014798, 0.2752759, -0.4578880, -0.3894024, 0.6349581, -0.1471152, -0.2173081)
expected <- c(0.3093946, 0.2812113, -0.4154854, -0.3527661, 0.6445740, -0.1908238, -0.2761046)
test_that("Is the PCA PC2 as expected?", {
    expect_equal(expected, actual, tolerance = 0.01)
})

tsne_stuff <- plot_tsne(norm, seed = 1)
tsne_stuff$plot
actual <- tsne_stuff[["table"]][["Factor1"]]
##expected <- c(117.7296, 117.7064, -287.8868, -287.8774, 105.9297, 114.4243, 119.9742)
expected <- c(176.7191, 156.0065, 165.2553, 165.8340, -225.6521, -220.7042, -217.4586)
## These values seem to have changed in the new version of Rtsne.
test_that("Is the tsne data as expected for Comp1?", {
    expect_equal(expected, actual, tolerance = 0.1)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 03graph_metrics.R in ", elapsed, " seconds.")
