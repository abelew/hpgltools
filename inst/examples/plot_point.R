if (is.null(pombe_se)) {
  source(system.file("examples", "create_pombe.R", package = "hpgltools"))
}

## Biological Coefficient of variation plots are most helpful for raw data.
bcv <- plot_bcv(pombe_se)
bcv[["plot"]]

if (is.null(pombe_norm)) {
  source(system.file("examples", "normalize_pombe.R", package = "hpgltools"))
}

## Since I have the data loaded, I will plot two samples
## of my pombe data against each other...
test_df <- as.data.frame(assay(pombe_norm[, c(1, 2)]))

test_scatter <- plot_dist_scatter(test_df)
test_scatter

boring_scatter <- plot_scatter(test_df)
boring_scatter

test_scatter <- plot_linear_scatter(test_df)
summary(test_scatter)
test_scatter[["scatter"]]
## Hey, check it, they are very similar.

## See how many nonzero genes there are per sample vs. coverage
nonzero <- plot_nonzero(pombe_se)
nonzero

## MA samples against each other
## In this case I am just comparing the first two samples
## If I did not do this subset, it would do all pairwise comparisons.
sample_ma <- plot_pairwise_ma(pombe_se[, c(1, 2)])
sample_ma
