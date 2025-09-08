#' Extract groups of samples from a SE and make a boxplot of their expression
#'
#' @param exp Input dataset
#' @param conditions Pair of conditions to compare (hopefully more than a pair soon)
#' @param genes Vector of genes to plot
#' @param norm Normalization to apply
#' @param convert Conversion to apply
#' @param filter Filter to apply
#' @param transform Transformation
#' @param batch Batch method to apply
#' @param name_column Gene information column from which to get gene names.
#' @param condition_column Metadata column containing the conditions to compare.
#' @param test Pairwise test to invoke.
#' @importFrom dplyr arrange
#' @export
ggsignif_paired_genes <- function(exp, conditions = NULL, genes = NULL, norm = "quant",
                                  convert = "cpm", filter = TRUE, transform = "log2",
                                  batch = "svaseq", name_column = "hgnc_symbol",
                                  condition_column = "condition", test = "wilcox.test") {
  normed <- normalize(exp, filter = filter, convert = convert, transform = transform,
                      norm = norm, batch = batch)
  starting_meta <- colData(exp)
  subset_idx <- starting_meta[[condition_column]] %in% conditions
  exp_subset <- normed[, subset_idx]
  annot <- as.data.frame(rowData(exp_subset))
  meta <- as.data.frame(colData(exp_subset))

  row_idx <- NULL
  wanted_annotations <- data.frame()
  wanted_expression <- data.frame()
  if (is.null(genes)) {
    genes <- head(rowData(normed)[[name_column]], n = 10)
    row_idx <- seq_len(10)
  } else {
    all_ids <- annot[[name_column]]
    row_idx <- all_ids %in% genes
  }
  wanted_annotations <- as.data.frame(annot[row_idx, ]) %>%
    dplyr::arrange(factor(!!sym(name_column), levels = genes))
  wanted_rows <- rownames(wanted_annotations)
  wanted_expression <- as.data.frame(assay(normed)[wanted_rows, ])

  start_df <- merge(wanted_expression, wanted_annotations, by = "row.names")
  plot_df <- start_df %>%
    reshape2::melt() %>%
    as.data.frame()
  plot_df <- merge(plot_df, meta, by.x = "variable", by.y = "row.names")
  plot_df[["pair"]] <- paste0(plot_df[[name_column]], "_", plot_df[[condition_column]])

  comparison_list <- list()
  axis_labels <- c()
  level_order <- c()
  for (i in genes) {
    element <-  c(glue("{i}_{conditions[1]}"), glue("{i}_{conditions[2]}"))
    comparison_list[[i]] <- element
    axis_labels <- c(axis_labels, i, "")
    level_order <- c(level_order, element)
  }
  pair <- value <- condition <- NULL
  merged[["pair"]] <- factor(merged[["pair"]], levels = level_order)
  plot <- ggplot(data = merged, aes(x = pair, y = value, fill = condition)) +
    ggplot2::geom_boxplot() +
    ggsignif::geom_signif(comparisons = comparison_list, step_increase = 0.01, test = test) +
    ggplot2::scale_x_discrete(labels = axis_labels) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
  return(plot)
}
