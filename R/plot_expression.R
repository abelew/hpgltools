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
ggsignif_paired_genes <- function(exp, conditions = NULL, genes = NULL, norm = "raw",
                                  convert = "cpm", filter = TRUE, transform = "log2",
                                  batch = "svaseq", name_column = "hgnc_symbol",
                                  condition_column = "condition", test = "wilcox.test",
                                  excel = NULL) {
  normed <- normalize(exp, filter = filter, convert = convert, transform = transform,
                      norm = norm, batch = batch)
  starting_meta <- colData(exp)
  subset_idx <- starting_meta[[condition_column]] %in% conditions
  exp_subset <- normed[, subset_idx]
  plot_colors <- get_colors_by_condition(exp_subset)
  start_subset <- assay(exp)[, subset_idx]
  annot <- as.data.frame(rowData(exp_subset))
  meta <- as.data.frame(colData(exp_subset))

  row_idx <- NULL
  wanted_annotations <- data.frame()
  start_expression <- data.frame()
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
  start_expression <- as.data.frame(start_subset[wanted_rows, ])
  start_df <- merge(start_expression, wanted_annotations, by = "row.names")
  norm_df <- merge(wanted_expression, wanted_annotations, by = "row.names")
  plot_df <- norm_df %>%
    reshape2::melt() %>%
    as.data.frame()
  plot_df <- merge(plot_df, meta, by.x = "variable", by.y = "row.names")
  plot_df[["pair"]] <- paste0(plot_df[[name_column]], "_", plot_df[[condition_column]])
  raw_df <- start_df %>%
    reshape2::melt() %>%
    as.data.frame()
  raw_df <- merge(raw_df, meta, by.x = "variable", by.y = "row.names")
  raw_df[["pair"]] <- paste0(raw_df[[name_column]], "_", raw_df[[condition_column]])

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
  plot_df[["pair"]] <- factor(plot_df[["pair"]], levels = level_order)
  raw_df[["pair"]] <- factor(raw_df[["pair"]], levels = level_order)
  plot <- ggplot(data = plot_df, aes(x = pair, y = value, fill = condition)) +
    ggplot2::geom_boxplot() +
    ggsignif::geom_signif(comparisons = comparison_list, step_increase = 0.01, test = test) +
    ggplot2::scale_x_discrete(labels = axis_labels) +
    ggplot2::scale_fill_manual(name = "Condition", values = plot_colors) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
  raw_plot <- ggplot(data = raw_df, aes(x = pair, y = value, fill = condition)) +
    ggplot2::geom_boxplot() +
    ggsignif::geom_signif(comparisons = comparison_list, step_increase = 0.01, test = test) +
    ggplot2::scale_x_discrete(labels = axis_labels) +
    ggplot2::scale_fill_manual(name = "Condition", values = plot_colors) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))

  if (!is.null(excel)) {
    xlsx <- init_xlsx(excel)
    wb <- xlsx[["wb"]]
    excel_basename <- xlsx[["basename"]]
    new_row <- 1
    new_col <- 1
    ## Set up a vector of images to clean up when finished.
    image_files <- c()
    ## Get the number of columns in the plot to define width:
    plot_width <- (length(plot_colors) * length(genes)) / 2
    plot_height <- 6
    xls_result <- write_xlsx(data = start_df, wb = wb, start_row = new_row,
                             rownames = FALSE, sheet = "raw", start_col = 1,
                             title = "Raw reads.")
    xls_result <- write_xlsx(data = norm_df, wb = wb, start_row = new_row,
                             rownames = FALSE, sheet = "norm", start_col = 1,
                             title = "Normalized reads.")
    try_result <- xlsx_insert_png(plot, wb = wb, sheet = "norm_plot",
                                  width = plot_width, height = plot_height,
                                  start_col = 1, start_row = 1,
                                  plotname = "01_norm", savedir = excel_basename,
                                  fancy_type = "svg")
    if ("try-error" %in% class(try_result)) {
      warning("Failed to add the normalized plot.")
    } else {
      image_files <- c(image_files, try_result[["filename"]])
    }
    try_result <- xlsx_insert_png(raw_plot, wb = wb, sheet = "raw_plot",
                                  width = plot_width, height = plot_height,
                                  start_col = 1, start_row = 1,
                                  plotname = "01_norm", savedir = excel_basename,
                                  fancy_type = "svg")
    if ("try-error" %in% class(try_result)) {
      warning("Failed to add the raw plot.")
    } else {
      image_files <- c(image_files, try_result[["filename"]])
    }
    save_result <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
    message("Saving to ", excel)
    for (img in image_files) {
      removed <- try(suppressWarnings(file.remove(img)), silent = TRUE)
    }
  }
  retlist <- list(
    "raw_df" = start_df,
    "norm_df" = norm_df,
    "raw_plot" = raw_plot,
    "plot" = plot)
  class(retlist) <- "paired_expression_plot"
  return(retlist)
}

## EOF
