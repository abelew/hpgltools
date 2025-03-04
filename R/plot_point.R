## plot_scatter.r: Various scatter plots

#' Steal edgeR's plotBCV() and make it a ggplot2.
#'
#' This was written primarily to understand what that function is doing in edgeR.
#'
#' @param data Dataframe/expt/exprs with count data
#' @return Plot of the BCV a la ggplot2.
#' @seealso [edgeR::plotBCV()] [ggplot2]
#' @examples
#' \dontrun{
#'  bcv <- plot_bcv(expt)
#'  summary(bcv$data)
#'  bcv$plot
#' }
#' @export
plot_bcv <- function(data) {
  data_class <- class(data)[1]
  if (data_class == "expt" || data_class == "SummarizedExperiment") {
    data <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    data <- exprs(data)
  } else if (data_class == "matrix" || data_class == "data.frame") {
    data <- as.data.frame(data)
    ## some functions prefer matrix, so I am keeping this explicit for the moment
  } else {
    stop("This function only understands types: expt, ExpressionSet, data.frame, and matrix.")
  }
  data <- edgeR::DGEList(counts = data)
  edisp <- edgeR::estimateDisp(data)
  avg_log_cpm <- edisp[["AveLogCPM"]]
  if (is.null(avg_log_cpm)) {
    avg_log_cpm <- edgeR::aveLogCPM(edisp[["counts"]], offset = edgeR::getOffset(edisp))
  }
  disper <- edgeR::getDispersion(edisp)
  if (is.null(disper)) {
    stop("No dispersions to plot")
  }
  if (attr(disper, "type") == "common") {
    disper <- rep(disper, length = length(avg_log_cpm))
  }
  disp_df <- data.frame("A" = avg_log_cpm,
                        "disp" = sqrt(disper))
  fitted_disp <- gplots::lowess(disp_df[["A"]], disp_df[["disp"]], f = 0.5)
  f <- stats::approxfun(fitted_disp, rule = 2)
  disp_df[["label"]] <- rownames(disp_df)
  disp_plot <- ggplot(
    disp_df, aes(x = .data[["A"]], y = .data[["disp"]],
                 label = .data[["label"]])) +
    ggplot2::geom_point() +
    ggplot2::xlab("Average log(CPM)") +
    ggplot2::ylab("Dispersion of Biological Variance") +
    ggplot2::stat_density2d(geom = "tile", aes(fill = ggplot2::after_stat(density ^ 0.25)),
                            contour = FALSE, show.legend = FALSE) +
    ggplot2::scale_fill_gradientn(
      colours = grDevices::colorRampPalette(c("white", "black"))(256)) +
    ggplot2::geom_smooth(method = "loess") +
    ggplot2::stat_function(fun = f, colour = "red") +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.position = "none",
                   axis.text = ggplot2::element_text(size = base_size, colour = "black"))
  ret <- list("data" = disp_df, "plot" = disp_plot)
  return(ret)
}

#' Make a scatter plot between two sets of numbers with a cheesy distance metric
#' and some statistics of the two sets.
#'
#' The distance metric should be codified and made more intelligent.
#' Currently it creates a dataframe of distances which are absolute
#' distances from each axis, multiplied by each other, summed by axis,
#' then normalized against the maximum.
#'
#' @param df Dataframe likely containing two columns.
#' @param size Size of the dots.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
#' @return Ggplot2 scatter plot.  This plot provides a "bird's eye"
#'  view of two data sets.  This plot assumes the two data structures
#'  are not correlated, and so it calculates the median/mad of each
#'  axis and uses these to calculate a stupid, home-grown distance
#'  metric away from both medians.  This distance metric is used to
#'  color dots which are presumed the therefore be interesting because
#'  they are far from 'normal.'  This will make a fun clicky googleVis
#'  graph if requested.
#' @seealso [ggplot2::geom_point()] [plot_linear_scatter()]
#' @examples
#' \dontrun{
#'  dist_scatter(lotsofnumbers_intwo_columns)
#' }
#' @export
plot_dist_scatter <- function(df, size = 2, xlab = NULL, ylab = NULL) {
  df <- data.frame(df[, c(1, 2)])
  df <- df[complete.cases(df), ]
  df_columns <- colnames(df)
  df_x_axis <- df_columns[1]
  df_y_axis <- df_columns[2]
  if (is.null(xlab)) {
    xlab <- glue("Expression of {df_x_axis}")
  }
  if (is.null(ylab)) {
    ylab <- glue("Expression of {df_y_axis}")
  }
  colnames(df) <- c("first", "second")
  first_median <- summary(df[, 1])["Median"]
  second_median <- summary(df[, 2])["Median"]
  first_mad <- stats::mad(df[, 1])
  second_mad <- stats::mad(df[, 2])
  mydist <- sillydist(df[, 1], df[, 2], first_median, second_median)
  mydist[["x"]] <- abs((mydist[, 1] - first_median) / abs(first_median))
  mydist[["y"]] <- abs((mydist[, 2] - second_median) / abs(second_median))
  mydist[["x"]] <- mydist[["x"]] / max(mydist[["x"]])
  mydist[["y"]] <- mydist[["y"]] / max(mydist[["y"]])
  mydist[["dist"]] <- mydist[["x"]] * mydist[["y"]]
  mydist[["dist"]] <- mydist[["dist"]] / max(mydist[["dist"]])
  line_size <- size / 2
  df[["label"]] <- rownames(df)
  first_vs_second <- ggplot(
    df, aes(x = .data[["first"]], y = .data[["second"]],
            label = .data[["label"]])) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::geom_vline(
      color = "grey", xintercept = (first_median - first_mad), size = line_size) +
    ggplot2::geom_vline(
      color = "grey", xintercept = (first_median + first_mad), size = line_size) +
    ggplot2::geom_vline(
      color = "darkgrey", xintercept = first_median, size = line_size) +
    ggplot2::geom_hline(
      color = "grey", yintercept = (second_median - second_mad), size = line_size) +
    ggplot2::geom_hline(
      color = "grey", yintercept = (second_median + second_mad), size = line_size) +
    ggplot2::geom_hline(color = "darkgrey", yintercept = second_median, size = line_size) +
    ggplot2::geom_point(
      colour = grDevices::hsv(mydist[["dist"]], 1, mydist[["dist"]]),
      alpha = 0.6, size = size) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.position = "none",
                   axis.text = ggplot2::element_text(size = base_size, colour = "black"))
  return(first_vs_second)
}

#' Make a scatter plot between two groups with a linear model superimposed and
#' some supporting statistics.
#'
#' @param df Dataframe likely containing two columns.
#' @param cormethod What type of correlation to check?
#' @param size Size of the dots on the plot.
#' @param loess Add a loess estimation?
#' @param xcol Column name of x-values
#' @param ycol Column name of y-values#'
#' @param text_col Column containing text annotations.
#' @param logfc Point out genes with a specific logfc.
#' @param identity Add the identity line?
#' @param z Use this z-score cutoff.
#' @param z_lines  Include lines defining the z-score boundaries.
#' @param first First column to plot.
#' @param second Second column to plot.
#' @param base_url Base url to add to the plot.
#' @param pretty_colors Colors!
#' @param xlab Alternate x-axis label.
#' @param ylab Alternate x-axis label.
#' @param color_high Chosen color for points significantly above the mean.
#' @param color_low Chosen color for points significantly below the mean.
#' @param alpha Choose an alpha channel to define how see-through the dots are.
#' @param ... Extra args likely used for choosing significant genes.
#' @return List including a ggplot2 scatter plot and some histograms.  This plot
#'  provides a "bird's eye" view of two data sets.  This plot assumes a
#'  (potential) linear correlation between the data, so it calculates the
#'  correlation between them.  It then calculates and plots a robust linear
#'  model of the data using an 'SMDM' estimator (which I don't remember how to
#'  describe, just that the document I was reading said it is good).  The
#'  median/mad of each axis is calculated and plotted as well.  The distance
#'  from the linear model is finally used to color the dots on the plot.
#'  Histograms of each axis are plotted separately and then together under a
#'  single cdf to allow tests of distribution similarity.  This will make a fun
#'  clicky googleVis graph if requested.
#' @seealso [robust] [stats] [ggplot2] [robust::lmRob] [stats::weights] [plot_histogram()]
#' @examples
#' \dontrun{
#'  plot_linear_scatter(lotsofnumbers_intwo_columns)
#' }
#' @export
plot_linear_scatter <- function(df, cormethod = "pearson", size = 2, loess = FALSE,
                                xcol = NULL, ycol = NULL, text_col = NULL, logfc = 2.0,
                                identity = FALSE, z = 1.5, z_lines = FALSE,
                                first = NULL, second = NULL, base_url = NULL,
                                color_weights = TRUE, xlab = NULL, ylab = NULL,
                                model_type = "robust", add_equation = TRUE, add_rsq = TRUE,
                                add_cor = TRUE, label_prefix = "Expression of",
                                color_high = NULL, color_low = NULL, alpha = 0.4, ...) {
  ## At this time, one might expect arglist to contain
  ## z, p, fc, n and these will therefore be passed to get_sig_genes()
  arglist <- list(...)
  if (isTRUE(color_high)) {
    color_high <- "#FF0000"
  }
  if (isTRUE(color_low)) {
    color_low <- "#7B9F35"
  }
  if (is.null(xcol)) {
    xcol <- colnames(df)[1]
    ycol <- colnames(df)[2]
  }

  correlation <- try(cor.test(df[[xcol]], df[[ycol]], method = cormethod, exact = FALSE))
  cor_value <- correlation[["estimate"]]
  if (class(correlation)[1] == "try-error") {
    correlation <- NULL
    cor_value <- NULL
  }
  df_columns <- colnames(df)
  if (is.null(xlab)) {
    xlab <- glue("{label_prefix} {xcol}")
  }
  if (is.null(ylab)) {
    ylab <- glue("{label_prefix} {ycol}")
  }
  test_formula <- as.formula(glue("{ycol} ~ {xcol}"))
  linear_model <- NULL
  linear_model_summary <- NULL
  linear_model_rsq <- NULL
  linear_model_weights <- NULL
  linear_model_intercept <- NULL
  linear_model_slope <- NULL
  if (model_type == "robust") {
    model_test <- try(robustbase::lmrob(formula = test_formula,
                                        data = df, method = "SMDM"), silent = TRUE)
  } else if (model_type == "glm") {
    model_test <- try(glm(formula = test_formula, data = df), silent = TRUE)
  } else if (model_type == "lm") {
    model_test <- try(lm(formula = test_formula, data = df), silent = TRUE)
  } else {
    stop("I don't yet know this model type.")
  }
  if (class(model_test)[1] == "try-error") {
    warning("Model type ", model_type, " failed, falling back to a default lm.")
    model_test <- try(lm(formula = test_formula, data = df), silent = TRUE)
    if (class(model_test)[1] == "try-error") {
      message("Could not create a linear model of the data.")
      message("Going to perform a scatter plot without linear model.")
      plot <- plot_scatter(df)
      ret <- list(data = df, scatter = plot)
      return(ret)
    }
  }
  linear_model <- model_test
  linear_model_summary <- summary(linear_model)
  linear_model_rsq <- linear_model_summary[["r.squared"]]
  linear_model_weights <- stats::weights(linear_model, type = "robustness", na.action = NULL)
  if (is.null(linear_model_weights)) {
    color_weights = FALSE
  }
  linear_model_intercept <- stats::coef(linear_model_summary)[1]
  linear_model_slope <- stats::coef(linear_model_summary)[2]
  first_median <- summary(df[[xcol]])[["Median"]]
  second_median <- summary(df[[ycol]])[["Median"]]
  first_mad <- stats::mad(df[[xcol]], na.rm = TRUE)
  second_mad <- stats::mad(df[[ycol]], na.rm = TRUE)
  line_size <- size / 2
  df[["label"]] <- rownames(df)
  .data <- NULL ## aes figured this out via NSE shenanigans.
  if (is.null(text_col)) {
    aesthetics <- aes(x = .data[[xcol]], y = .data[[ycol]],
                      label = .data[["label"]])
  } else {
    aesthetics <- aes(x = .data[[xcol]], y = .data[[ycol]],
                      label = .data[["label"]], text = .data[[text_col]])
  }
  first_vs_second <- ggplot(df, mapping = aesthetics) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::geom_vline(
      color = "grey", xintercept = (first_median - first_mad), size = line_size) +
    ggplot2::geom_vline(
      color = "grey", xintercept = (first_median + first_mad), size = line_size) +
    ggplot2::geom_hline(
      color = "grey", yintercept = (second_median - second_mad), size = line_size) +
    ggplot2::geom_hline(
      color = "grey", yintercept = (second_median + second_mad), size = line_size) +
    ggplot2::geom_hline(
      color = "darkgrey", yintercept = second_median, size = line_size) +
    ggplot2::geom_vline(
      color = "darkgrey", xintercept = first_median, size = line_size) +
    ggplot2::geom_abline(
      colour = "grey", slope = linear_model_slope,
      intercept = linear_model_intercept, size = line_size)
  if (isTRUE(identity)) {
    first_vs_second <- first_vs_second +
      ggplot2::geom_abline(colour = "darkgreen", slope = 1, intercept = 0, size = 0.8)
  }

  ## The axes and guide-lines are set up, now add the points
  low_df <- high_df <- NULL
  if (!is.null(color_low) || !is.null(color_high)) {
    ## If you want to color the above or below identity line points, then you
    ## will need subsets to define them
    tmpdf <- df
    tmpdf[["ratio"]] <- tmpdf[[ycol]] - tmpdf[[xcol]]
    subset_points <- get_sig_genes(tmpdf, column = "ratio", lfc = logfc, z = z)
    high_subset <- subset_points[["up_genes"]]
    low_subset <- subset_points[["down_genes"]]
    original_df <- tmpdf
    high_index <- rownames(original_df) %in% rownames(high_subset)
    high_df <- original_df[high_index, ]
    low_index <- rownames(original_df) %in% rownames(low_subset)
    low_df <- original_df[low_index, ]
    first_vs_second <- first_vs_second +
      ggplot2::geom_point(colour = "black", size = size, alpha = alpha)
  }

  if (isTRUE(z_lines)) {
    first_vs_second <- first_vs_second +
      ggplot2::geom_abline(colour = "grey", slope = linear_model_slope,
                           intercept = linear_model_intercept + z, size = line_size / 3) +
      ggplot2::geom_abline(colour = "grey", slope = linear_model_slope,
                           intercept = linear_model_intercept - z, size = line_size / 3)
  }

  ## Add a color to the dots which are lower than the identity line by some amount
  if (!is.null(color_low)) {
    first_vs_second <- first_vs_second +
      ggplot2::geom_point(data = low_df, colour = color_low)
  }
  if (!is.null(color_high)) {
    first_vs_second <- first_vs_second +
      ggplot2::geom_point(data = high_df, colour = color_high)
  }

  if (isTRUE(color_weights)) {
    first_vs_second <- first_vs_second +
      ggplot2::geom_point(size = size, alpha = alpha,
                          colour = grDevices::hsv(linear_model_weights * (9 / 20),
                                                  linear_model_weights / 20 + (19 / 20),
                                                  (1.0 - linear_model_weights)))
  } else {
    first_vs_second <- first_vs_second +
      ggplot2::geom_point(colour = "black", size = size, alpha = alpha)
  }

  annot_string <- glue("
")
  if (isTRUE(add_equation)) {
    annot_string <- glue("{annot_string}
  Equation: {ymxb_print(linear_model)}", .trim = FALSE)
  }
  if (isTRUE(add_rsq)) {
    annot_string <- glue("{annot_string}
  R^2: {signif(x=linear_model_rsq, digits=3)}", .trim = FALSE)
  }
  if (isTRUE(add_cor)) {
    annot_string <- glue("{annot_string}
  {cormethod} correlation: {signif(x=cor_value, digits=3)}", .trim = FALSE)
  }
  if (!is.null(annot_string)) {
    first_vs_second <- first_vs_second +
      ggplot2::annotate("text", x = -Inf, y = Inf, label = annot_string, vjust = 1, hjust = 0)
  }

  if (isTRUE(loess)) {
    first_vs_second <- first_vs_second +
      ggplot2::geom_smooth(method = "loess")
  }

  first_vs_second <- first_vs_second +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.position = "none",
                   axis.text = ggplot2::element_text(size = base_size, colour = "black"))

  x_histogram <- plot_histogram(data.frame(df[[xcol]]), fillcolor = "lightblue", color = "blue")
  y_histogram <- plot_histogram(data.frame(df[[ycol]]), fillcolor = "pink", color = "red")
  both_histogram <- plot_multihistogram(df[ , c(xcol, ycol)])
  plots <- list(
    "data" = df,
    "scatter" = first_vs_second,
    "x_histogram" = x_histogram,
    "y_histogram" = y_histogram,
    "both_histogram" = both_histogram,
    "correlation" = correlation,
    "lm_model" = linear_model,
    "lm_summary" = linear_model_summary,
    "lm_weights" = linear_model_weights,
    "lm_rsq" = linear_model_rsq,
    "first_median" = first_median,
    "first_mad" = first_mad,
    "second_median" = second_median,
    "second_mad" = second_mad)
  class(plots) <- "linear_scatter"
  return(plots)
}

#' Quick point-recolorizer given an existing plot, df, list of rownames to
#' recolor, and a color.
#'
#' This function should make it easy to color a family of genes in any of the
#' point plots.
#'
#' @param plot Geom_point based plot
#' @param df Data frame used to create the plot
#' @param ids Set of ids which must be in the rownames of df to recolor
#' @param color Chosen color for the new points.
#' @param ... Extra arguments are passed to arglist.
#' @return prettier plot.
recolor_points <- function(plot, df, ids, color = "red", ...) {
  arglist <- list(...)
  alpha <- 0.3
  if (!is.null(arglist[["alpha"]])) {
    alpha <- arglist[["alpha"]]
  }

  point_index <- rownames(df) %in% ids
  newdf <- df[point_index, ]
  newplot <- plot +
    ggplot2::geom_point(data = newdf,  colour = color, fill = color, alpha = alpha)
  return(newplot)
}

#' Make a ggplot graph of the number of non-zero genes by sample.
#'
#' This puts the number of genes with > 0 hits on the y-axis and CPM on the
#' x-axis. Made by Ramzi Temanni <temanni at umd dot edu>.
#'
#' @param data Expt, expressionset, or dataframe.
#' @param design Eesign matrix.
#' @param colors Color scheme.
#' @param plot_labels How do you want to label the graph? 'fancy' will use
#'  directlabels() to try to match the labels with the positions without
#'  overlapping anything else will just stick them on a 45' offset next to the
#'  graphed point.
#' @param expt_names Column or character list of preferred sample names.
#' @param max_overlaps Permit this many labels to overlap before dropping some.
#' @param label_chars How many characters for sample names before abbreviation.
#' @param plot_legend Print a legend for this plot?
#' @param plot_title Add a title?
#' @param cutoff Minimum proportion (or number) of genes below which samples might be in trouble.
#' @param ... rawr!
#' @return a ggplot2 plot of the number of non-zero genes with respect to each
#'  library's CPM.
#' @seealso [ggplot2]
#' @examples
#' \dontrun{
#'  nonzero_plot <- plot_nonzero(expt = expt)
#' }
#' @export
plot_nonzero <- function(data, design = NULL, colors = NULL, plot_labels = "repel",
                         expt_names = NULL, max_overlaps = 5, label_chars = 10, plot_legend = FALSE,
                         plot_title = NULL, cutoff = 0.65, ...) {
  arglist <- list(...)

  condition <- design[["condition"]]
  batch <- design[["batch"]]
  if (!is.null(expt_names) && class(expt_names)[1] == "character") {
    if (length(expt_names) == 1) {
      colnames(data) <- make.names(design[[expt_names]], unique = TRUE)
    } else {
      colnames(data) <- expt_names
    }
  }
  if (!is.null(label_chars) && is.numeric(label_chars)) {
    colnames(data) <- abbreviate(colnames(data), minlength = label_chars)
  }
  nz_df <- data.frame(
    "id" = colnames(data),
    "nonzero_genes" = colSums(data > 0),
    "cpm" = colSums(data) * 1e-6,
    "condition" = condition,
    "batch" = batch,
    "color" = as.character(colors))

  ## Add a little logic to warn the user if samples have poor representation
  ## using a cutoff which may either be a proportion of the number of available
  ## rows, or an aribtrary cutoff
  sad_samples <- NULL
  if (!is.null(cutoff)) {
    if (cutoff < 1) { ## Then it is a proportion
      cutoff <- nrow(data) * cutoff
    }
    sad_idx <- nz_df[["nonzero_genes"]] <= cutoff
    sad_samples <- nz_df[sad_idx, "id"]
    if (length(sad_samples) > 0) {
      message("The following samples have less than ", cutoff, " genes.")
      print(sad_samples)
    }
  }

  color_listing <- nz_df[, c("condition", "color")]
  color_listing <- unique(color_listing)
  color_list <- as.character(color_listing[["color"]])
  names(color_list) <- as.character(color_listing[["condition"]])
  nz_df[["label"]] <- rownames(nz_df)
  num_batches <- length(unique(nz_df[["batch"]]))

  non_zero_plot <- ggplot(
    data = nz_df,
    aes(x = .data[["cpm"]], y = .data[["nonzero_genes"]], label = .data[["label"]]))
  if (num_batches <= 5) {
    non_zero_plot <- non_zero_plot +
      ggplot2::geom_point(size = 3,
                          aes(shape = .data[["batch"]],
                              colour = as.factor(.data[["condition"]]),
                              fill = as.factor(.data[["condition"]]))) +
      ggplot2::geom_point(size = 3, colour = "black", show.legend = FALSE,
                          aes(shape = .data[["batch"]],
                              fill = as.factor(.data[["condition"]]))) +
      ggplot2::scale_color_manual(name = "Condition",
                                  values = color_list) +
      ggplot2::scale_fill_manual(name = "Condition",
                                 values = color_list) +
      ggplot2::scale_shape_manual(
        name = "Batch",
        labels = levels(as.factor(nz_df[["batch"]])),
        values = 21:25)
  } else {
    non_zero_plot <- non_zero_plot +
      ggplot2::geom_point(size = 3, shape = 21,
                          aes(colour = as.factor(.data[["condition"]]),
                              fill = as.factor(.data[["condition"]]))) +
      ggplot2::geom_point(size = 3, shape = 21, colour = "black", show.legend = FALSE,
                          aes(fill = as.factor(.data[["condition"]])))
  }

  non_zero_plot <- non_zero_plot +
    ggplot2::scale_color_manual(name = "Condition",
                                guide = "legend",
                                values = color_list) +
    ggplot2::scale_fill_manual(name = "Condition",
                               guide = "legend",
                               values = color_list) +
    ggplot2::ylab("Number of non-zero genes observed") +
    ggplot2::xlab("Number of reads mapped (millions)") +
    ggplot2::theme_bw(base_size = base_size)

  if (isTRUE(plot_labels)) {
    plot_labels <- "repel"
  }
  if (is.null(plot_labels)) {
    plot_labels <- "repel"
  }
  if (plot_labels == FALSE) {
    message("Not putting labels on the plot.")
  } else if (plot_labels == "normal") {
    non_zero_plot <- non_zero_plot +
      ggplot2::geom_text(aes(x = .data[["cpm"]], y = .data[["nonzero_genes"]],
                             label = .data[["id"]], angle = 45, size = 4, vjust = 2))
  } else if (plot_labels == "oldrepel") {
    non_zero_plot <- non_zero_plot +
      ggrepel::geom_text_repel(ggplot2::aes(label = .data[["id"]]),
                               size = 5, box.padding = ggplot2::unit(0.5, "lines"),
                               point.padding = ggplot2::unit(1.6, "lines"),
                               arrow = ggplot2::arrow(length = ggplot2::unit(0.01, "npc")))
  } else if (plot_labels == "dlsmart") {
    non_zero_plot <- non_zero_plot +
      directlabels::geom_dl(aes(label = .data[["id"]]), method = "smart.grid")
  } else if (plot_labels == "repel") {
    non_zero_plot <- non_zero_plot +
      ggrepel::geom_text_repel(ggplot2::aes(label = .data[["id"]]), max.overlaps = max_overlaps)
  } else {
    non_zero_plot <- non_zero_plot +
      directlabels::geom_dl(ggplot2::aes(label = .data[["id"]]), method = "first.qp")
  }

  if (!is.null(plot_title)) {
    non_zero_plot <- non_zero_plot + ggplot2::ggtitle(plot_title)
  }
  non_zero_plot <- non_zero_plot +
    ggplot2::theme(axis.ticks = ggplot2::element_blank(),
                   axis.text = ggplot2::element_text(size = base_size, colour = "black"))
  if (isFALSE(plot_legend)) {
    non_zero_plot <- non_zero_plot +
      ggplot2::theme(legend.position = "none")
  }

  retlist <- list(
    "plot" = non_zero_plot,
    "table" = nz_df)
  class(retlist) <- "nonzero_plot"
  return(retlist)
}
setGeneric("plot_nonzero")

#' Plot all pairwise MA plots in an experiment.
#'
#' Use affy's ma.plot() on every pair of columns in a data set to help diagnose
#' problematic samples.
#'
#' @param data Expt expressionset or data frame.
#' @param log Is the data in log format?
#' @param ... Options are good and passed to arglist().
#' @return List of affy::maplots
#' @seealso [affy::ma.plot()]
#' @examples
#' \dontrun{
#'  ma_plots = plot_pairwise_ma(expt = some_expt)
#' }
#' @export
plot_pairwise_ma <- function(data, log = NULL, ...) {
  data_class <- class(data)[1]
  if (data_class == "expt" || data_class == "SummarizedExperiment") {
    design <- pData(data)
    colors <- data[["colors"]]
    data <- exprs(data)
  } else if (data_class == "ExpressionSet") {
    data <- exprs(data)
  } else if (data_class == "matrix" || data_class == "data.frame") {
    ## some functions prefer matrix, so I am keeping this explicit for the moment
    data <- as.data.frame(data)
  } else {
    stop("This function understands types: expt, ExpressionSet, data.frame, and matrix.")
  }
  plot_list <- list()
  for (c in seq(from = 1, to = length(colnames(data)) - 1)) {
    nextc <- c + 1
    for (d in seq(from = nextc, to = length(colnames(data)))) {
      first <- as.numeric(data[, c])
      second <- as.numeric(data[, d])
      if (max(first) > 1000) {
        if (is.null(log)) {
          message("I suspect you want to set log = TRUE for this.")
          message("In fact, I am so sure, I am doing it now.")
          message("If I am wrong, set log = FALSE, but I'm not.")
          log <- TRUE
        }
      } else if (max(first) < 80) {
        if (!is.null(log)) {
          message("I suspect you want to set log = FALSE for this.")
          message("In fact, I am so  sure, I am doing it now.")
          message("If I am wrong, set log = TRUE.")
          log <- FALSE
        }
      }
      firstname <- colnames(data)[c]
      secondname <- colnames(data)[d]
      name <- glue("{firstname}_{secondname}")
      if (isTRUE(log)) {
        first <- log2(first + 1.0)
        second <- log2(second + 1.0)
      }
      m <- first - second
      a <- (first + second) / 2

      tmp_file <- tmpmd5file(pattern = "ma", fileext = ".png")
      this_plot <- png(filename = tmp_file)
      controlled <- dev.control("enable")
      affy::ma.plot(A = a, M = m, plot.method = "smoothScatter",
                    show.statistics = TRUE, add.loess = TRUE)
      title(glue("MA of {firstname} vs {secondname}."))
      plot_list[[name]] <- grDevices::recordPlot()
      dev.off()
      removed <- suppressWarnings(file.remove(tmp_file))
      removed <- unlink(dirname(tmp_file))

    }
  }
  return(plot_list)
}

#' Make a pretty scatter plot between two sets of numbers.
#'
#' This function tries to supplement a normal scatterplot with some information
#' describing the relationship between the columns of data plotted.
#'
#' @param df Dataframe likely containing two columns.
#' @param size Size of the dots on the graph.
#' @param color Color of the dots on the graph.
#' @param xlab Alternate x-axis label.
#' @param ylab Alternate x-axis label.
#' @param alpha Define how see-through the dots are.
#' @return Ggplot2 scatter plot.
#' @seealso [plot_linear_scatter()] [all_pairwise()]
#' @examples
#' \dontrun{
#'  plot_scatter(lotsofnumbers_intwo_columns)
#' }
#' @export
plot_scatter <- function(df, color = "black", xlab = NULL, xcol = NULL, ycol = NULL,
                         ylab = NULL, alpha = 0.6, size = 2) {
  if (is.null(xcol)) {
    xcol <- 1
  }
  if (is.null(ycol)) {
    ycol <- 2
  }
  df <- data.frame(df[, c(xcol, ycol)])
  df <- df[complete.cases(df), ]
  df_columns <- colnames(df)
  df_x_axis <- df_columns[1]
  df_y_axis <- df_columns[2]
  if (is.null(xlab)) {
    xlab <- glue("Expression of {df_x_axis}")
  }
  if (is.null(ylab)) {
    ylab <- glue("Expression of {df_y_axis}")
  }
  colnames(df) <- c("first", "second")
  df[["label"]] <- rownames(df)
  first_vs_second <- ggplot(df, aes(x = .data[["first"]], y = .data[["second"]],
                                    label = .data[["label"]])) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::geom_point(colour = color, alpha = alpha, size = size) +
    ggplot2::theme(legend.position = "none",
                   axis.text = ggplot2::element_text(size = 10, colour = "black"))
  return(first_vs_second)
}

## EOF
