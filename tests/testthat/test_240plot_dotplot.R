start <- as.POSIXlt(Sys.time())
context("240plot_dotplot.R")

pombe_se <- make_pombe_se(annotation = FALSE)
pombe_subset <- subset_se(
  pombe_se,
  subset = "minute == 0 | minute == 15 | minute == 30")
pombe_norm <- normalize(pombe_se, transform = "log2", convert = "cpm", filter = TRUE)

pombe_estimates <- all_adjusters(pombe_norm, estimate_type = "sva")
adjustments <- pombe_estimates[["model_adjust"]]

test <- plot_svfactor(pombe_norm, adjustments)
gg_class <- "ggplot2::ggplot"
actual <- class(test)[1]
test_that("Can we plot surrogate variables by factor?", {
  expect_equal(gg_class, actual)
})

test <- plot_batchsv(pombe_norm, adjustments)
actual1 <- class(test[[1]])[1]
actual2 <- class(test[[2]])[1]
actual3 <- class(test[[3]])[1]
test_that("Can we plot surrogates by batch?", {
  expect_equal(gg_class, actual1)
  expect_equal(gg_class, actual2)
  expect_equal(gg_class, actual3)
})

pombe_pca <- plot_pca(pombe_norm)
pombe_pcs <- pombe_pca[["table"]]
test <- plot_pcfactor(pombe_pcs, pombe_norm)
actual <- class(test)[1]
test_that("Can we plot principle components by factor?", {
  expect_equal(gg_class, actual)
})

test <- plot_sm(pombe_norm)
actual <- class(test[["plot"]])[1]
test_that("Can we plot standard medians?", {
  expect_equal(gg_class, actual)
})

end <- as.POSIXlt(Sys.time())
elapsed <- round(x = as.numeric(end - start))
message("\nFinished 240plot_dotplot.R in ", elapsed,  " seconds.")
