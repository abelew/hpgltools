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
#' @param padding not currently used, but add some padding to the plot
#' @param span degree of smoothing
#' @param smoothed Smooth the coverage plot?
#' @param ribbon fill in the smoothed data?
#' @param cov_color Color for the coverage region
#' @param y_max Explicitly set the maximum of the y axis.
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
    cov_plot <- ggcoverage(data = g1_coverage, color = "white") +
      geom_smooth(data = g1_coverage, se = FALSE, method = "loess", color = cov_color,
                  method.args = list(family = "symmetric"), span = span,
                  aes(x = start, y = score, ymin = 0))
  } else {
    cov_plot <- ggcoverage(data = g1_coverage, color = cov_color)
  }
  if (isTRUE(ribbon)) {
    cov_plot <- cov_plot +
      geom_ribbon(data = g1_coverage, outline.type = "full", show.legend = FALSE,
                  aes(x = start, y = score, ymin = 0,
                      ymax = predict(loess(score ~ start, span = span))),
                  fill = cov_color)
  }

  if (is.null(y_max)) {
    cov_plot <- cov_plot +
      scale_y_continuous(minor_breaks = scales::breaks_width(1),
                         position = "left", expand = c(0, 0))
  } else {
    cov_plot <- cov_plot +
      scale_y_continuous(minor_breaks = scales::breaks_width(1), limits = c(0, y_max),
                         position = "left", exapnd = c(0, 0))
  }

  cov_plot <- cov_plot +
    theme(axis.ticks.length = unit(10, "pt"),
          axis.minor.ticks.length = rel(0.5),
          axis.text = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 10, hjust = 1),
          axis.ticks.y = element_line(colour = "darkblue")) +
    geom_gene(gtf.gr = g1, gene.size = 4.0, label.size = 4.0,
              overlap.style = "tight", arrow.type = "closed",
              arrow.gap = NULL, arrow.num = NULL)
  retlist <- list(
    "region" = gr,
    "coverage_data" = g1_coverage,
    "plot" = cov_plot)
  return(retlist)
}

## EOF
