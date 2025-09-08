#' Try plotting a chromosome (region)
#'
#' genoplotr is cool, I don't yet understand it though
#'
#' @param accession An accession to plot, this will download it.
#' @param start First segment to plot (doesn't quite work yet).
#' @param end Final segment to plot (doesn't quite work yet).
#' @param plot_title Put a title on the resulting plot.
#' @return Hopefully a pretty plot of a genome
#' @seealso [genoPlotR]
#' @export
genoplot_chromosome <- function(accession = "AE009949", start = NULL, end = NULL,
                                plot_title = "Genome plot") {
  tt <- download_gbk(accession)
  segments <- try(genoPlotR::read_dna_seg_from_file(glue("{accession}.gb")))
  if (is.null(start)) {
    start <- 1
  }
  if (is.null(end)) {
    end <- nrow(segments)
  }

  mid_pos <- genoPlotR::middle(segments)
  xlims <- list(c(Inf, -Inf), c(-Inf, Inf), c(start, end))
  genoPlotR::plot_gene_map(dna_segs = list(segments),
                           main = plot_title,
                           gene_type = "side_blocks",
                           dna_seg_scale = TRUE,
                           scale = FALSE)
}

#' Plot coverage from a grange using ggcoverage
#'
#' ggcoverage makes a bunch of assumptions about how the input gff is formatted.
#'
#' @param gr input granges
#' @param from starting point for the plot
#' @param to ending point
#' @param id_column mcols() portion of the granges containing the gene names
#' @param coverage_files Set of bigwig file used to define the coverage.
#' @param padding not currently used, but add some padding to the plot
#' @param span degree of smoothing
#' @param smoothed Smooth the coverage plot?
#' @param ribbon fill in the smoothed data?
#' @param cov_color Color for the coverage region
#' @param y_max Explicitly set the maximum of the y axis.
#' @importFrom stats predict loess
#' @export
plot_ggcoverage <- function(gr, from, to, id_column, coverage_files, padding = 100,
                            span = 0.03, smoothed = TRUE, ribbon = FALSE, cov_color = "black",
                            y_max = NULL) {
  g1 <- gr_from_to(gr, from, to, id_column)
  dup_idx <- duplicated(mcols(g1))
  g1 <- g1[!dup_idx, ]
  region_string <- paste0(seqnames(range(g1)[1]), ":",
                          min(start(g1)) - padding, "-", max(start(g1)) + padding)
  g1_coverage <- ggcoverage::LoadTrackFile(track.file = coverage_files,
                                           gtf.gr = g1,
                                           format = "bw",
                                           region = region_string)
  if (isTRUE(smoothed)) {
    cov_plot <- ggcoverage::ggcoverage(data = g1_coverage, color = "white") +
      ggplot2::geom_smooth(data = g1_coverage, se = FALSE, method = "loess", color = cov_color,
                  method.args = list(family = "symmetric"), span = span,
                  aes(x = start, y = score, ymin = 0))
  } else {
    cov_plot <- ggcoverage::ggcoverage(data = g1_coverage, color = cov_color)
  }
  if (isTRUE(ribbon)) {
    cov_plot <- cov_plot +
      ggplot2::geom_ribbon(data = g1_coverage, outline.type = "full", show.legend = FALSE,
                  aes(x = start, y = score, ymin = 0,
                      ymax = stats::predict(stats::loess(score ~ start, span = span))),
                  fill = cov_color)
  }

  if (is.null(y_max)) {
    cov_plot <- cov_plot +
      ggplot2::scale_y_continuous(minor_breaks = scales::breaks_width(1),
                                  position = "left", expand = c(0, 0))
  } else {
    cov_plot <- cov_plot +
      ggplot2::scale_y_continuous(minor_breaks = scales::breaks_width(1), limits = c(0, y_max),
                                  position = "left", exapnd = c(0, 0))
  }

  cov_plot <- cov_plot +
    ggplot2::theme(axis.ticks.length = grid::unit(10, "pt"),
                   axis.minor.ticks.length = ggplot2::rel(0.5),
                   axis.text = ggplot2::element_text(color = "black", size = 12),
                   axis.text.y = ggplot2::element_text(color = "black", size = 10, hjust = 1),
                   axis.ticks.y = ggplot2::element_line(colour = "darkblue")) +
    ggcoverage::geom_gene(gtf.gr = g1, gene.size = 4.0, label.size = 4.0,
                          overlap.style = "tight", arrow.type = "closed",
                          arrow.gap = NULL, arrow.num = NULL)
  retlist <- list(
    "region" = gr,
    "coverage_data" = g1_coverage,
    "plot" = cov_plot)
  class(retlist) <- "coverage_plot"
  return(retlist)
}

#' Print a gg coverage plot
#'
#' FIXME: Add something useful to this
#' @param x List containing some information and the plot of interest.
#' @param ... Other args for the generic.
#' @export
print.coverage_plot <- function(x, ...) {
  plot(x[["plot"]])
  return(invisible(x))
}

#' Given a SE, plot a coverage map with some helpful defaults and workarounds.
#'
#' ggcoverage is really cool, but has some hard-coded values which are
#' not useful. This function seeks to overcome some of these
#' limitations and fill in some default values.
#'
#' @param se Input summarized experiment, it should have a full genome
#'  GRanges and colors.
#' @param from Starting gene region's name or number.
#' @param to Ending gene region's name or number.
#' @param id_column Column of the rowData containing the gene IDs.
#' @param coverage_column Column of the colData containing the
#'  location of the coverage files in bigwig format.
#' @param convert Perform a cpmish conversion of the data?
#' @param norm Normalize the data?
#' @param transform log transform the data?
#' @param padding Add this amount of sequence before/after the
#'  granges.
#' @param span When smoothing, adjust the loess transformation.
#' @param smoothed Smooth the coverage data? (currently only loess)
#' @param ribbon Add a ribbon to the smoothing. (currently disabled)
#' @param cov_color Specify colors for the coverage lines.
#' @param y_max Specify a y-maximum.
#' @param bin_size Specify the coverage bin size.
#' @param hide_bars Hide the coverage bars when smoothing the data?
#' @param feature_type_column Use this rowData column to extract the
#'  feature types.
#' @param feature_type Use this type for creating the GR.
#' @param meta_type_column Use this colData column to get the
#'  ggcoverage Type column.
#' @param feature_name_column Use this rowData column to fill in the
#'  gene names below the gene arrows.
#' @param gene_name Choose a gene?
#' @param meta_group_column Use this colData column to group samples.
#' @param overlap_type Define the behavior of the arrows/row.
#' @param overlap_gap Modify the arrow overlaps with this.
#' @param arrow_type Either closed or open.
#' @param arrow_gap Not sure what this does.
#' @param label_column FIXME: perhaps redundant with feature_name_column?
#' @param bg_color Define the plot background color.
#' @param gene_color_by Play with the gene arrow colors (I want to
#'  modify ggcoverage to allow me to color by FC status).
#' @param gene_colors specify factor of gene colors.
#' @param facet Use a specifc factor to facet the data.
#' @param y_scale Set a specifc y-axis scale.
plot_ggcoverage_se <- function(se, from = 1 , to = 10, id_column = "gene_id",
                               coverage_column = "deeptools_coverage", convert = "cpm",
                               norm = NULL, transform = NULL,
                               padding = c(1, 2), span = 0.03, smoothed = TRUE, ribbon = FALSE,
                               cov_color = NULL, y_max = NULL, bin_size = 10, hide_bars = FALSE,
                               feature_type_column = "type",
                               feature_type = "gene", meta_type_column = "media",
                               feature_name_column = "Parent",
                               gene_name = NULL, meta_group_column = "condition",
                               overlap_type = "tight", overlap_gap = 0.1,
                               arrow_type = "closed", arrow_gap = NULL, label_column = NULL,
                               bg_color = "white", gene_color_by = "strand", gene_colors = NULL,
                               facet = TRUE, y_scale = NULL) {
  se <- set_conditions(se, fact = meta_group_column)
  se_meta <- colData(se)
  se_info <- rowData(se)
  se_exprs <- assay(se)
  gr <- metadata(se)[["grange"]]
  colors_by_cond <- get_colors_by_condition(se)

  if (is.null(gene_colors)) {
    gene_colors <- c("-" = "cornflowerblue", "+" = "darkolivegreen3")
  }
  bs <- metadata(se)[["genome"]]
  padding_start <- 0
  padding_end <- 0
  if (length(padding) == 1) {
    padding_start <- padding
    padding_end <- padding
  } else if (length(padding) == 2) {
    padding_start <- padding[1]
    padding_end <- padding[2]
  }
  coverage_files <- se_meta[[coverage_column]]
  message("Getting grange subset with padding, start: ", padding_start, " end: ", padding_end, ".")
  gr_subset <- gr_from_to(gr, from, to, column = id_column, type_column = feature_type_column,
                          type = feature_type, padding_start = padding_start,
                          padding_end = padding_end)
  gr_start <- min(start(gr_subset))
  gr_end <- max(end(gr_subset))
  mesg("The region from ", from, " to ", to, " starts at ", gr_start, " and ends at ", gr_end, ".")
  dup_idx <- duplicated(mcols(gr_subset))
  if (sum(dup_idx) > 0) {
    mesg("There were ", sum(dup_idx), " duplicates in this region.")
    gr_subset <- gr_subset[!dup_idx, ]
  }
  region_string <- paste0(seqnames(range(gr_subset)[1]), ":",
                          min(start(gr_subset)), "-", max(start(gr_subset)))
  ## Match the expectations of ggcoverage
  convert = tolower(convert)
  coverage_info = load_se_tracks(se, track_column = coverage_column,
                                 region_string = region_string, name_column = feature_name_column,
                                 gene_name = gene_name, padding = padding,
                                 bin_size = bin_size, convert = convert, transform = transform,
                                 norm = norm,
                                 type_column = meta_type_column,
                                 group_column = meta_group_column, cores = cores)
  coverage_info <- as.data.frame(coverage_info)
  coverage_info[["shared"]] <- "all"

  ## algG is from 1586560 to 1588191
  test_position <- 1587000
  region <- coverage_info[["end"]] == test_position
  coverage_info[region, ][, c("sampleid.y", "condition", "score")]
  ## We should see larger numbers for the IPTG samples than not.
  ## This is true.

  ## Define the rows in the plot
  facet_key <- "shared"
  if (!is.null(facet)) {
    if (isTRUE(facet)) {
      facet_key = "Type"
    } else if ("character" %in% class(facet)) {
      facet_key <- facet
    }
  }
  message("Using ", facet_key, " as the facet.")

  max_coverage <- max(coverage_info[["score"]])
  mesg("The maximum observed coverage in this region is: ", max_coverage, ".")
  if (isTRUE(smoothed)) {
    hidden_colors <- rep(bg_color, length(colors_by_cond))
    smoothed_info <- stats::loess(score ~ start, data = coverage_info, family = "symmetric", span = span)
    loess_max <- max(smoothed_info[["fitted"]])
    mesg("The maximum observed coverage in the smoothed region is: ", loess_max, ".")
    max_coverage <- loess_max
    Type <- NULL
    cov_plot <- ggcoverage::ggcoverage(data = coverage_info, color = hidden_colors,
                                       facet.key = facet_key, range.position = "out") +
      ggplot2::geom_smooth(data = coverage_info, se = FALSE, method = "loess",
                           method.args = list(family = "symmetric"), span = span,
                           show.legend = FALSE,
                           aes(x = start, y = score, ymin = 0, color = Type))
  } else {
    cov_plot <- ggcoverage::ggcoverage(data = coverage_info, color = cov_color)
  }

  cov_plot <- cov_plot +
    ggplot2::scale_color_manual(values = colors_by_cond, guide = "none")
  if (!isTRUE(hide_bars)) {
    cov_plot <- cov_plot +
      ggplot2::scale_fill_manual(values = colors_by_cond, guide = "none")
  }

  ## Play with the strip colors via ggh4x?
  #if (facet_key != "shared") {
  #  strip_ggh4x <- ggh4x::strip_themed(
  #    background_y = ggh4x::elem_list_rect(
  #      fill = colors_by_cond,
  #      color = colors_by_cond))
  #  strip_form <- as.formula(glue("~ {meta_group_column}"))
  #  cov_plot <- cov_plot +
  #    ggh4x::facet_wrap2(strip_form,
  #                       nrow = length(colors_by_cond), ncol = 1,
  #                       strip.position = "right",
  #                       strip = strip_ggh4x)
  #}
  if (is.null(y_max)) {
    cov_plot <- cov_plot +
      ggplot2::scale_y_continuous(minor_breaks = scales::breaks_width(1), limits = c(0, max_coverage),
                                  position = "right", expand = c(0, 0))
  } else {
    cov_plot <- cov_plot +
      ggplot2::scale_y_continuous(minor_breaks = scales::breaks_width(1), limits = c(0, y_max),
                                  position = "left", expand = c(0, 0))
  }
  if (!is.null(y_scale) && y_scale == "log2") {
    cov_plot <- cov_plot +
      ggplot2::scale_y_continuous(trans = "log2")
  }

  ## ggcoverage hard-codes the label column to 'gene'
  if (!is.null(label_column)) {
    message("Changing the label column from gene_name to ", label_column, ".")
    ## Note, mcols may have the label column as some weird list thing like 'CompressedCharacterList'
    mcols(gr_subset)[["gene_name"]] <- as.character(mcols(gr_subset)[[label_column]])
  }

  cov_plot <- cov_plot +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bg_color, colour = NA),
      axis.ticks.length = grid::unit(10, "pt"),
      axis.minor.ticks.length = ggplot2::rel(0.5),
      axis.text = ggplot2::element_text(color = "black", size = 12),
      axis.text.y = ggplot2::element_text(color = "black", size = 10, hjust = 1),
      axis.ticks.y = ggplot2::element_line(colour = "darkblue")) +
    ggcoverage::geom_gene(gtf.gr = gr_subset, gene.size = 4.0, label.size = 4.0,
                          overlap.style = overlap_type, arrow.type = arrow_type,
                          overlap.gene.gap = overlap_gap, arrow.gap = arrow_gap,
                          arrow.num = NULL, fill.color = gene_colors,
                          color.by = gene_color_by)

  region_start <- min(start(gr_subset), na.rm = TRUE)
  region_end <- max(end(gr_subset), na.rm = TRUE)
  region_width <- region_end - region_start
  observed_genes <- mcols(gr_subset)[[feature_name_column]]
  retlist <- list(
    "region" = gr_subset,
    "region_start" = region_start,
    "region_end" = region_end,
    "region_width" = region_width,
    "coverage_data" = coverage_info,
    "observed_genes" = observed_genes,
    "plot" = cov_plot)
  class(retlist) <- "ggcoverage_plot"
  return(retlist)
}

#' @export
print.ggcoverage_plot <- function(x, ...) {
  message("Something something useful.")
  print(x[["plot"]])
  return(invisible(x))
}

## EOF
