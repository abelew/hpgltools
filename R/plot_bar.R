#' Create a plot showing relative expression with respect to each chromosome/contig.
#'
#' @param expt Input expressionset.
#' @param chromosome_column Annotation column containing the chromosome ID.
#' @param scaffolds Include scaffolds in addition to the actual chromosomes.
#' @param min_genes The minimum number of genes which should be on the 'chromosome' before it
#'  is considered worth considering.
#' @export
plot_exprs_by_chromosome <- function(expt, chromosome_column = "chromosome", scaffolds = TRUE,
                                     min_genes = 10) {
  start <- data.frame(row.names = unique(fData(expt)[["chromosome"]]))
  start[["genes"]] <- 0
  start[["exprs_mean"]] <- 0
  start[["exprs_stdev"]] <- 0
  start[["exprs_var"]] <- 0
  start[["exprs_min"]] <- 0
  start[["exprs_qt1"]] <- 0
  start[["exprs_qt3"]] <- 0
  start[["exprs_median"]] <- 0
  start[["exprs_max"]] <- 0
  for (ch in rownames(start)) {
    gene_id_idx <- fData(expt)[["chromosome"]] == ch
    gene_ids <- rownames(fData(expt))[gene_id_idx]
    start[ch, "genes"] <- length(gene_ids)
    subset_exprs <- exprs(expt)[gene_ids, ]
    if (length(gene_ids) == 1) {
      start[ch, "exprs_mean"] <- mean(subset_exprs)
      start[ch, "exprs_stdev"] <- stats::sd(subset_exprs)
      start[ch, "exprs_min"] <- min(subset_exprs)
      start[ch, "exprs_max"] <- max(subset_exprs)
      start[ch, "exprs_median"] <- median(subset_exprs)
      start[ch, "exprs_var"] <- stats::var(subset_exprs)
      start[ch, "exprs_qt1"] <- as.numeric(summary(subset_exprs))[2]
      start[ch, "exprs_qt3"] <- as.numeric(summary(subset_exprs))[5]
    } else {
      start[ch, "exprs_mean"] <- mean(rowMeans(subset_exprs))
      start[ch, "exprs_stdev"] <- stats::sd(rowMeans(subset_exprs))
      start[ch, "exprs_min"] <- min(rowMeans(subset_exprs))
      start[ch, "exprs_max"] <- max(rowMeans(subset_exprs))
      start[ch, "exprs_median"] <- median(rowMeans(subset_exprs))
      start[ch, "exprs_var"] <- stats::var(rowMeans(subset_exprs))
      start[ch, "exprs_qt1"] <- as.numeric(summary(rowMeans(subset_exprs)))[2]
      start[ch, "exprs_qt3"] <- as.numeric(summary(rowMeans(subset_exprs)))[5]
    }

    min_idx <- start[["genes"]] >= min_genes
    start <- start[min_idx, ]

    plt <- ggplot(start, aes(y = exprs_mean, x = genes)) +
      ggplot2::geom_point() +
      ggplot2::scale_y_log10()
    plt <- ggplot(start, aes(y = exprs_var, x = genes)) +
      ggplot2::geom_point()
  }
  retlist <- list(
    "plot" = plt,
    "info" = start)
  return(retlist)
}

#' Make a ggplot graph of library sizes.
#'
#' It is often useful to have a quick view of which samples have more/fewer
#' reads.  This does that and maintains one's favorite color scheme and tries to
#' make it pretty!
#'
#' @param data Expt, dataframe, or expressionset of samples.
#' @param condition Vector of sample condition names.
#' @param colors Color scheme if the data is not an expt.
#' @param text Add the numeric values inside the top of the bars of the plot?
#' @param order Explicitly set the order of samples in the plot?
#' @param plot_title Title for the plot.
#' @param yscale Whether or not to log10 the y-axis.
#' @param expt_names Design column or manually selected names for printing sample names.
#' @param label_chars Maximum number of characters before abbreviating sample names.
#' @param ... More parameters for your good time!
#' @return a ggplot2 bar plot of every sample's size
#' @seealso [ggplot2] [prettyNum] [plot_sample_bars()]
#' @examples
#' \dontrun{
#'  libsize_plot <- plot_libsize(expt = expt)
#'  libsize_plot  ## ooo pretty bargraph
#' }
#' @export
plot_libsize <- function(data, condition = NULL, colors = NULL,
                         text = TRUE, order = NULL, plot_title = NULL, yscale = NULL,
                         expt_names = NULL, label_chars = 10,
                         ...) {
  arglist <- list(...)
  if (is.null(text)) {
    text <- TRUE
  }

  ## In response to Keith's recent comment when there are more than 8 factors
  chosen_palette <- "Dark2"
  if (!is.null(arglist[["palette"]])) {
    chosen_palette <- arglist[["palette"]]
  }

  mtrx <- as.matrix(data)
  if (is.null(colors)) {
    pal <- RColorBrewer::brewer.pal(ncol(mtrx), chosen_palette)
    colors <- grDevices::colorRampPalette(pal)(ncol(mtrx))
  }

  ## Get conditions
  if (is.null(condition)) {
    stop("Missing condition label vector.")
  }
  values <- as.numeric(mtrx)
  integerp <- all.equal(values, as.integer(values))

  colors <- as.character(colors)
  sum <- NULL

  if (!is.null(expt_names) && class(expt_names) == "character") {
    if (length(expt_names) == 1) {
      colnames(mtrx) <- make.names(design[[expt_names]], unique = TRUE)
    } else {
      colnames(mtrx) <- expt_names
    }
  }
  if (!is.null(label_chars) && is.numeric(label_chars)) {
    colnames(mtrx) <- abbreviate(colnames(mtrx), minlength = label_chars)
  }

  libsize_df <- data.frame("id" = as.character(colnames(mtrx)),
                           "sum" = colSums(mtrx),
                           "condition" = condition,
                           "colors" = as.character(colors))
  summary_df <- data.table::setDT(libsize_df)[, list(
    "min" = min(sum),
    "1st" = quantile(x = sum, probs = 0.25),
    "median" = median(x = sum),
    "mean" = mean(sum),
    "3rd" = quantile(x = sum, probs = 0.75),
    "max" = max(sum)),
    by = "condition"]
  libsize_plot <- plot_sample_bars(
    libsize_df, condition = condition, colors = colors,
    text = text, order = order, plot_title = plot_title, integerp = integerp,
    yscale = yscale, ...)
  retlist <- list(
    "plot" = libsize_plot,
    "table" = libsize_df,
    "summary" = summary_df)
  class(retlist) <- "libsize_plot"
  return(retlist)
}
setGeneric("plot_libsize")

#' Visualize genes observed before/after filtering.
#'
#' Thanks to Sandra Correia for this!  This function attempts to represent the
#' change in the number of genes which are well/poorly represented in the data
#' before and after performing a low-count filter.
#'
#' @param expt Input expressionset.
#' @param low_limit Threshold to define 'low-representation.'
#' @param filter Method used to low-count filter the data.
#' @param ... Extra arbitrary arguments to pass to normalize_expt()
#' @return Bar plot showing the number of genes below the low_limit before and
#'  after filtering the data.
#' @seealso [plot_libsize()] [filter_counts()]
#' @export
plot_libsize_prepost <- function(expt, low_limit = 2, filter = TRUE,
                                 num_color = "black", num_size = 4, ...) {
  start <- plot_libsize(expt, text = FALSE)
  norm <- sm(normalize_expt(expt, filter = filter,
                            ...))
  end <- plot_libsize(norm)

  ## Gather the number of genes which are <= the low limit, before and after filtering
  lt_min_start <- colSums(exprs(expt) <= low_limit)
  lt_min_end <- colSums(exprs(norm) <= low_limit)
  start_tab <- as.data.frame(start[["table"]])
  end_tab <- as.data.frame(end[["table"]])
  ## Get the number of counts in each samples.
  start_tab[["sum"]] <- as.numeric(start_tab[["sum"]])
  start_tab[["colors"]] <- as.character(start_tab[["colors"]])
  start_tab[["alpha"]] <- ggplot2::alpha(start_tab[["colors"]], 0.75)
  ## Get the number of low genes in each sample.
  start_tab[["low"]] <- lt_min_start
  ##start_tab[["sub_low"]] <- 0
  start_tab[["subtraction"]] <- 0
  start_tab[["subtraction_string"]] <- ""

  ## Get the number of counts after filtering.
  end_tab[["sum"]] <- as.numeric(end_tab[["sum"]])
  end_tab[["colors"]] <- as.character(end_tab[["colors"]])
  end_tab[["alpha"]] <- ggplot2::alpha(end_tab[["colors"]], 1.0)
  ## Find how many counts were lost from filtering.
  subtract_count_sums <- start_tab[["sum"]] - end_tab[["sum"]]
  end_tab[["subtraction"]] <- subtract_count_sums
  ## Get the number of genes after filtering.
  end_tab[["low"]] <- lt_min_end
  ## end_tab[["sub_low"]] <- 0
  end_tab[["subtraction_string"]] <- ""

  ## Get the number of genes lost from filtering.
  subtract_gene_sums <- start_tab[["low"]] - end_tab[["low"]]
  undef <- is.na(subtract_gene_sums)
  ## When this comes up NA, just use the original number.
  subtract_gene_sums[undef] <- start_tab[undef, "low"]
  ## start_tab[["sub_low"]] <- subtract_gene_sums
  start_tab[["subtraction"]] <- subtract_gene_sums
  start_tab[["subtraction_string"]] <- paste0(subtract_count_sums, " counts from ",
                                              subtract_gene_sums, " genes.")
  all_tab <- rbind(start_tab, end_tab)

  count_title <- glue("Counts remaining after filtering less than {low_limit} reads,
labeled by counts/genes removed.")
  ## Suppressing 'Using alpha for a discrete variable is not advised.'
  count_columns <- ggplot(all_tab,
                          aes(x = .data[["id"]], y = .data[["sum"]])) +
    ggplot2::geom_col(position = "identity", color = "black",
                      aes(fill = .data[["colors"]])) +
    ggplot2::scale_fill_manual(values = c(levels(as.factor(all_tab[["colors"]])))) +
    ggplot2::geom_text(parse = FALSE, angle = 90, size = num_size,
                       color = num_color, hjust = 1.2,
                       aes(x = .data[["id"]],
                           ## label='as.character(all_tab$subtraction)')) +
                           label = .data[["subtraction_string"]])) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 10, colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
                   legend.position = "none") +
    ggplot2::scale_y_continuous(trans = "log10", labels = scales::scientific) +
    ggplot2::ggtitle(count_title)
  low_title <- glue("Genes with less than {low_limit} reads, labeled by delta.")
  low_columns <- ggplot(all_tab, aes(x = .data[["id"]], y = .data[["low"]])) +
    ggplot2::geom_col(position = "identity", color = "black",
                      aes(alpha = .data[["alpha"]], fill = .data[["colors"]])) +
    ggplot2::scale_fill_manual(values = c(levels(as.factor(all_tab[["colors"]])))) +
    ggplot2::geom_text(parse = FALSE, angle = 90, size = num_size,
                       color = num_color, hjust = 0,
                       aes(x = .data[["id"]], label = .data[["low"]])) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = 10, colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
                   legend.position = "none") +
    ggplot2::ggtitle(low_title)
  retlist <- list(
    "start" = start,
    "end" = end,
    "table" = all_tab,
    "count_plot" = count_columns,
    "lowgene_plot" = low_columns)
  class(retlist) <- "prepost_filter"
  return(retlist)
}

#' Make a ggplot graph of the percentage/number of reads kept/removed.
#'
#' The function expt_exclude_genes() removes some portion of the original reads.
#' This function will make it possible to see what is left.
#'
#' @param data Dataframe of the material remaining, usually expt$summary_table
#' @param row Row name to plot.
#' @param condition Vector of sample condition names.
#' @param colors Color scheme if the data is not an expt.
#' @param names Alternate names for the x-axis.
#' @param text Add the numeric values inside the top of the bars of the plot?
#' @param plot_title Title for the plot.
#' @param yscale Whether or not to log10 the y-axis.
#' @param ... More parameters for your good time!
#' @return a ggplot2 bar plot of every sample's size
#' @seealso [plot_sample_bars()]
#' @examples
#' \dontrun{
#'  kept_plot <- plot_pct_kept(expt_removed)
#'  kept_plot  ## ooo pretty bargraph
#' }
#' @export
plot_pct_kept <- function(data, row = "pct_kept", condition = NULL, colors = NULL,
                          names = NULL, text = TRUE, plot_title = NULL, yscale = NULL, ...) {
  arglist <- list(...)
  table <- data
  if (class(data) == "expt" || class(data) == "SummarizedExperiment") {
    table <- data[["summary_table"]]
  }
  ## In response to Keith's recent comment when there are more than 8 factors
  chosen_palette <- "Dark2"
  if (!is.null(arglist[["palette"]])) {
    chosen_palette <- arglist[["palette"]]
  }

  data_class <- class(data)[1]
  if (data_class == "expt" || data_class == "SummarizedExperiment") {
    condition <- pData(data)[["condition"]]
    colors <- data[["colors"]]
    names <- data[["names"]]
    data <- exprs(data)  ## Why does this need the simplifying
    ## method of extracting an element? (eg. data['expressionset'] does not work)
    ## that is _really_ weird!
  } else if (data_class == "ExpressionSet") {
    data <- exprs(data)
  } else if (data_class == "matrix" || data_class == "data.frame") {
    data <- as.data.frame(data)
  } else {
    stop("This function understands types: expt, ExpressionSet, data.frame, and matrix.")
  }

  if (is.null(colors)) {
    colors <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(ncol(data),
                               chosen_palette))(ncol(data))
  }
  if (is.null(condition)) {
    stop("Missing condition label vector.")
  }
  colors <- as.character(colors)
  kept_df <- data.frame("id" = colnames(data),
                        "sum" = table[row, ],
                        "condition" = condition,
                        "colors" = as.character(colors))
  kept_plot <- plot_sample_bars(kept_df, condition = condition, colors = colors,
                                names = names, text = text,
                                plot_title = plot_title, yscale = yscale, ...)
  return(kept_plot)
}

#' The actual library size plotter.
#'
#' This makes a ggplot2 plot of library sizes.
#'
#' @param sample_df Expt, dataframe, or expressionset of samples.
#' @param condition Vector of sample condition names.
#' @param colors Color scheme if the data is not an expt.
#' @param integerp Is this comprised of integer values?
#' @param order Explicitly set the order of samples in the plot?
#' @param text Add the numeric values inside the top of the bars of the plot?
#' @param plot_title Title for the plot.
#' @param yscale Whether or not to log10 the y-axis.
#' @param ... Used to catch random arguments which are unused here.
plot_sample_bars <- function(sample_df, condition = NULL, colors = NULL,
                             integerp = FALSE, order = NULL,
                             text = TRUE, plot_title = NULL, yscale = NULL, ...) {
  arglist <- list(...)

  y_label <- "Library size in pseudocounts."
  if (isTRUE(integerp)) {
    y_label <- "Library size in counts."
  }
  if (!is.null(arglist[["y_label"]])) {
    y_label <- arglist[["y_label"]]
  }

  sample_df[["order"]] <- factor(sample_df[["id"]], as.character(sample_df[["id"]]))
  if (is.null(order)) {
    ## Order it by sample names lexically
    lexical <- order(levels(sample_df[["order"]]))
    sample_df <- sample_df[lexical, ]
    ## The next two lines are an attempt to ensure that the factor
    ## levels make sense and are in the same order as the samples.
    new_order <- as.character(sample_df[["order"]])
    sample_df[["order"]] <- as.factor(new_order)
  } else {
    new_df <- data.frame()
    for (o in order) {
      matches <- grep(pattern = o, x = sample_df[["order"]])
      adders <- sample_df[matches, ]
      new_df <- rbind(new_df, adders)
    }
    sample_df <- new_df
    sample_df[["order"]] <- factor(sample_df[["order"]], as.character(sample_df[["order"]]))
  }

  color_listing <- sample_df[, c("condition", "colors")]
  color_listing <- unique(color_listing)
  color_list <- as.character(color_listing[["colors"]])
  names(color_list) <- as.character(color_listing[["condition"]])

  ## Set a cutoff of the sum of the rgbs of the bar colors, if the color is > than that cutoff
  ## Set the text color to black.
  ## I reversed the logic for these colors!
  ## sample_df[["text_color"]] <- "#ffffff"
  sample_df[["text_color"]] <- "#FFFFFF"
  text_lst <- color_int(sample_df[["colors"]])
  for (r in seq_len(nrow(sample_df))) {
    ## Some recommendations for color cutoffs from:
    ## https://stackoverflow.com/questions/3942878/how-to-decide-font-color-in-white-or-black-depending-on-background-color
    ## Thanks to Theresa
    color_sum <- (text_lst[["red"]][[r]] * 0.299) +
      (text_lst[["green"]][[r]] * 0.587) +
      (text_lst[["blue"]][[r]] * 0.114)
    if (color_sum >= 186) {
      ## Hahah I am a doofus and reversed the logic!
      sample_df[r, "text_color"] <- "#000000"
    }
  }
  text_colors <- sample_df[["text_color"]]

  sample_plot <- ggplot(data = sample_df,
                        colour = colors,
                        aes(x = .data[["order"]], y = .data[["sum"]])) +
    ggplot2::geom_bar(stat = "identity",
                      colour = "black",
                      fill = sample_df[["colors"]],
                      aes(x = .data[["order"]])) +
    ggplot2::xlab("Sample ID") +
    ggplot2::ylab(y_label) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = base_size, colour = "black"),
                   axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))

  if (isTRUE(text)) {
    if (!isTRUE(integerp)) {
      sample_df[["sum"]] <- sprintf("%.2f", round(as.numeric(sample_df[["sum"]]), 2))
    }

    sample_plot <- sample_plot +
      ggplot2::geom_text(
        parse = FALSE, angle = 90, size = 4, hjust = 1.2, colour = text_colors,
        aes(x = .data[["order"]],
            label = prettyNum(as.character(sample_df[["sum"]]), big.mark = ",")))
  }

  if (!is.null(plot_title)) {
    sample_plot <- sample_plot + ggplot2::ggtitle(plot_title)
  }
  if (is.null(yscale)) {
    scale_difference <- max(as.numeric(sample_df[["sum"]])) /
      min(as.numeric(sample_df[["sum"]]))
    if (scale_difference > 10.0) {
      mesg("The scale difference between the smallest and largest
libraries is > 10. Assuming a log10 scale is better, set scale = FALSE if not.")
      yscale <- TRUE
    } else {
      yscale <- FALSE
    }
  }
  if (isTRUE(yscale)) {
    sample_plot <- sample_plot +
      ggplot2::scale_y_continuous(trans = "log10",
                                  labels = scales::scientific)
  } else if (yscale == "log10") {
    sample_plot <- sample_plot +
      ggplot2::scale_y_continuous(trans = "log10",
                                  labels = scales::scientific)
  } else if (yscale == "log2") {
    sample_plot <- sample_plot +
      ggplot2::scale_y_continuous(trans = "log2",
                                  labels = scales::scientific)
  } else {
    sample_plot <- sample_plot +
      ggplot2::scale_y_continuous(labels = scales::scientific)
  }
  return(sample_plot)
}

#' Make relatively pretty bar plots of coverage in a genome.
#'
#' This was written for ribosome profiling coverage / gene.
#' It should however, work for any data with little or no modification, it was
#' also written when I was first learning R and when I look at it now I see a
#' few obvious places which can use improvement.
#'
#' @param input Coverage / position filename.
#' @param workdir Where to put the resulting images.
#' @param output Output image filename.
#' @param name Gene name to print at the bottom of the plot.
#' @param start Relative to 0, where is the gene's start codon.
#' @param end Relative to 0, where is the gene's stop codon.
#' @param strand Is this on the + or - strand? (+1/-1)
#' @param padding How much space to provide on the sides?
#' @return coverage plot surrounging the ORF of interest
#' @seealso [ggplot2]
plot_rpm <- function(input, workdir = "images", output = "01.svg", name = "LmjF.01.0010",
                     start = 1000, end = 2000, strand = 1, padding = 100) {

  mychr <- gsub(pattern = "\\.\\d+$", replacement = "", x = name, perl = TRUE)
  plotted_start <- start - padding
  plotted_end <- end + padding
  my_start <- start
  my_end <- end
  region_idx <- input[["chromosome"]] == mychr &
    input[["position"]] >= plotted_start &
    input[["position"]] <= plotted_end
  rpm_region <- rpm_region[region_idx, ]
  rpm_region <- rpm_region[, -1]
  rpm_region[["log"]] <- log2(rpm_region[["rpm"]] + 1)

  ## pre_start = subset(rpm_region, position < my_start)
  start_idx <- rpm_region[["position"]] < my_start
  pre_start <- rpm_region[start_idx, ]
  ## post_stop = subset(rpm_region, position > my_end)
  stop_idx <- rpm_region[["position"]] > my_end
  post_stop <- rpm_region[stop_idx, ]
  ## cds = subset(rpm_region, position >= my_start & position < my_end)
  cds_idx <- rpm_region[["position"]] >= my_start & rpm_region[["position"]] < my_end
  cds <- rpm_region[cds_idx, ]

  eval(substitute(
    expr = {
      stupid <- aes(y = 0, yend = 0, x = my_start, xend = my_end)
    },
    env <- list(my_start = my_start, my_end = my_end)))

  if (strand == "+") {
    gene_arrow <- grid::arrow(type = "closed", ends = "last")
  } else {
    gene_arrow <- grid::arrow(type = "closed", ends = "first")
  }
  xlabel_string <- glue("{name}: {my_start} to {my_end}")
  my_plot <- ggplot(rpm_region, aes(x = .data[["position"]], y = .data[["log"]])) +
    ggplot2::xlab(xlabel_string) +
    ggplot2::ylab("Log2(RPM) reads") +
    ggplot2::geom_bar(data = rpm_region, stat = "identity", fill = "black", colour = "black") +
    ggplot2::geom_bar(data = pre_start, stat = "identity", fill = "red", colour = "red") +
    ggplot2::geom_bar(data = post_stop, stat = "identity", fill = "red", colour = "red") +
    ggplot2::geom_segment(data = rpm_region, mapping = stupid,
                          arrow = gene_arrow, size = 2, color = "blue") +
    ggplot2::theme_bw(base_size = base_size)

  return(my_plot)
}

#' Plot significant genes by contrast with different colors for significance levels.
#'
#' This is my attempt to recapitulate some plots made in Laura and Najib's mbio
#' paper.  The goal of the plot is to show a few ranges of significance as
#' differently colored and stacked bars.  The colors are nice because Najib and
#' Laura chose them.
#'
#' @param ups Set of up-regulated genes.
#' @param downs Set of down-regulated genes.
#' @param maximum Maximum/minimum number of genes to display.
#' @param text Add text at the ends of the bars describing the number of genes >/< 0 fc.
#' @param color_list Set of colors to use for the bars.
#' @param color_names Categories associated with aforementioned colors.
#' @return weird significance bar plots
#' @seealso [ggplot2] [extract_significant_genes()]
#' @export
plot_significant_bar <- function(ups, downs, maximum = NULL, text = TRUE,
                                 color_list = c("lightcyan", "lightskyblue", "dodgerblue",
                                                "plum1", "orchid", "purple4"),
                                 color_names = c("a_up_inner", "b_up_middle", "c_up_outer",
                                                 "a_down_inner", "b_down_middle", "c_down_outer")) {
  choose_max <- function(u, d) {
    ## m is the maximum found in the ups/downs
    m <- 0
    ## which is extracted from ups and downs
    um <- max(as.numeric(u))
    dm <- max(as.numeric(d))
    ## choose the maximum by which is biggest!
    if (um >= dm) {
      m <- um
    } else {
      m <- dm
    }
    ## Figure out the number of digits in the number
    digits <- nchar(as.character(m))
    ## And the number of zeroes in it.
    num_zeroes <- digits - 1.0
    ## Add 1 to the first digit
    first_digit <- as.numeric(strsplit(x = as.character(m), split = "")[[1]][[1]]) + 1.0
    ## And set maximum to that number * 10 to the number of zeroes.
    maximum <- first_digit * (10 ^ num_zeroes)
    return(maximum)
  }

  up_sums <- list()
  down_sums <- list()
  comp_names <- ups[ups[["variable"]] == "a_up_inner", ][["comparisons"]]
  for (comp in seq_along(comp_names)) {
    comp_name <- comp_names[[comp]]
    idx <- ups[["comparisons"]] == comp_name
    up_sums[[comp_name]] <- sum(as.numeric(ups[idx, ][["value"]]))
    idx <- downs[["comparisons"]] == comp_name
    down_sums[[comp_name]] <- sum(as.numeric(downs[idx, ][["value"]])) * -1.0
  }

  if (is.null(maximum)) {
    maximum <- choose_max(up_sums, down_sums)
  }

  ## Try to ensure that ggplot orders my colors and bars in the specific order I want.
  ## holy ass crackers this is annoying and difficult to get correct,  as the
  ## ordering is (to my eyes) arbitrary.
  names(color_list) <- color_names
  levels(ups[["variable"]]) <- c("c_up_outer", "b_up_middle", "a_up_inner")
  levels(downs[["variable"]]) <- c("c_down_outer", "b_down_middle", "a_down_inner")
  sigbar_plot <- ggplot() +
    ggplot2::geom_col(
      data = ups,
      aes(x = .data[["comparisons"]], y = .data[["value"]],
          fill = .data[["variable"]])) +
    ggplot2::geom_col(
      data = downs,
      aes(x = .data[["comparisons"]], y = .data[["value"]],
          fill = .data[["variable"]])) +
    ggplot2::scale_fill_manual(values = color_list) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "none")

  if (isTRUE(text)) {
    for (comp in seq_along(comp_names)) {
      comp_name <- comp_names[[comp]]
      upstring <- as.character(up_sums[[comp_name]])
      downstring <- as.character(down_sums[[comp_name]])
      sigbar_plot <- sigbar_plot +
        ggplot2::annotate("text", x = comp, y = maximum, label = upstring, angle = -90) +
        ggplot2::annotate("text", x = comp, y = maximum * -1, label = downstring, angle = 90)
    }
  }
  return(sigbar_plot)
}

## EOF
