if (is.null(pombe_se)) {
  source(system.file("examples", "create_se.R", package = "hpgltools"))
}

## Most distribution plots are particularly informative with the raw data.

## Make a boxplot of the SE
a_boxplot <- plot_boxplot(pombe_se)
a_boxplot

## The analagous density plot of the samples
a_density <- plot_density(pombe_se)
a_density

pombe_subset <- pombe_se[, c(1, 2, 3, 4)]
## QQ plots comparing samples to each other
qq <- plot_qq_all(pombe_subset)
qq[["logs"]]

## Show the percentage of all reads comprised by the top-n genes.
## This might be smarter if moved to the log scale?
topn_genes <- plot_topn(pombe_subset)
topn_genes

## Create a violin of the coefficient of variance on a per-condition basis.
varcoef <- plot_variance_coefficients(pombe_se)
varcoef
