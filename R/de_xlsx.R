## de_xlsx.r: Attempt to export DE results to consistent xlsx reports. This is,
## in effect, a dumping ground for the various metrics which people have asked
## me to include when attempting to decide if a given DE result is
## 'significant'.  I therefore try in this file to find an appropriate place to
## dump each metric into an excel workbook.

## Make sure I pull the imports from 01_hpgltools.R
#' @include 01_hpgltools.R
NULL

#' Convert a vector of yes/no by DE method to a list.
#'
#' This compiles the set of possible methods to include in an
#' all_pairwise() from a series of booleans into a simpler list and
#' checks that the elements have some data that may be used.
#'
#' @param includes List of methods desired
#' @param apr Input pairwise result, not always needed.
#' @param methods In a few place I use a list called methods instead of includes,
#'  this should be removed.
#' @param ... Extra args, currently ignored.
#' @return List containing TRUE/FALSE for each method desired,
#'  depending on if we actually have the relevant data.
check_includes <- function(includes = NULL, apr = NULL, methods = NULL, ...) {
  arglist <- list(...)
  final <- list()
  if (is.null(includes)) {
    final <- list(
      "basic" = TRUE, "deseq" = TRUE, "dream" = TRUE,
      "ebseq" = TRUE, "edger" = TRUE, "limma" = TRUE,
      "noiseq" = TRUE)
  } else if ("character" %in% class(includes)) {
    for (char in includes) {
      final[[char]] <- TRUE
    }
  }
  if (!is.null(methods)) {
    if (!is.null(methods[["basic"]])) {
      final[["basic"]] <- methods[["basic"]]
    }
    if (!is.null(methods[["deseq"]])) {
      final[["deseq"]] <- methods[["deseq"]]
    }
    if (!is.null(methods[["ebseq"]])) {
      final[["ebseq"]] <- methods[["ebseq"]]
    }
    if (!is.null(methods[["edger"]])) {
      final[["edger"]] <- methods[["edger"]]
    }
    if (!is.null(methods[["limma"]])) {
      final[["limma"]] <- methods[["limma"]]
    }
    if (!is.null(methods[["noiseq"]])) {
      final[["noiseq"]] <- methods[["noiseq"]]
    }
    if (!is.null(methods[["dream"]])) {
      final[["dream"]] <- methods[["dream"]]
    }
  }
  if (!is.null(arglist[["do_basic"]])) {
    final[["basic"]] <- arglist[["do_basic"]]
  }
  if (!is.null(arglist[["do_deseq"]])) {
    final[["deseq"]] <- arglist[["do_deseq"]]
  }
  if (!is.null(arglist[["do_dream"]])) {
    final[["dream"]] <- arglist[["do_dream"]]
  }
  if (!is.null(arglist[["do_ebseq"]])) {
    final[["ebseq"]] <- arglist[["do_ebseq"]]
  }
  if (!is.null(arglist[["do_edger"]])) {
    final[["edger"]] <- arglist[["do_edger"]]
  }
  if (!is.null(arglist[["do_limma"]])) {
    final[["limma"]] <- arglist[["do_limma"]]
  }
  if (!is.null(arglist[["do_noiseq"]])) {
    final[["noiseq"]] <- arglist[["do_noiseq"]]
  }
  if (!is.null(apr)) {
    queries <- c("basic", "deseq", "ebseq", "edger", "dream", "limma", "noiseq")
    for (q in queries) {
      if ("try-error" %in% class(apr[[q]]) || is.null(apr[[q]])) {
        final[[q]] <- FALSE
      }
    }
  }

  ## EBSeq made an incompatible change in its most recent release.
  ## I unthinkingly changed my code to match it without considering the old
  ## bioconductor release.
  if (as.numeric(as.character(BiocManager::version())) < 3.18) {
    warning("I changed ebseq_pairwise for the new bioc release, it needs >= 3.18.")
    final[["ebseq"]] <- FALSE
  }
  return(final)
}

#' Combine portions of deseq/limma/edger table output.
#'
#' This hopefully makes it easy to compare the outputs from
#' limma/DESeq2/EdgeR on a table-by-table basis.
#'
#' @param apr Output from all_pairwise().
#' @param extra_annot Add some annotation information?
#' @param keepers List of reformatted table names to explicitly keep
#'  certain contrasts in specific orders and orientations.
#' @param excludes List of columns and patterns to use for excluding
#'  genes.
#' @param adjp Perhaps you do not want the adjusted p-values for
#'  plotting?
#' @param includes List comprised of the set of analyses to include.
#' @param add_plots Add plots to the end of the sheets with expression values?
#' @param loess Add time intensive loess estimation to plots?
#' @param plot_dim Number of inches squared for the plot if added.
#' @param compare_plots Add some plots comparing the results.
#' @param padj_type Add a consistent p adjustment of this type.
#' @param fancy Save a set of fancy plots along with the xlsx file?
#' @param lfc_cutoff In this context, only used for plotting volcano/MA plots.
#' @param p_cutoff In this context, used for volcano/MA plots.
#' @param excel_title Title for the excel sheet(s).  If it has the
#'  string 'YYY', that will be replaced by the contrast name.
#' @param increment_start When incrementing the table number for each contrast,
#'  look for this string and increment when it is found.  It should therefore
#'  be found in the excel_title.
#' @param start_worksheet_num Start writing data at this worksheet number.
#'  (in case you want to put other stuff in)
#' @param rda Write a rda file of the results.
#' @param rda_input Include the input all_pairwise() result in the rda?
#' @param label Label this number of top-n genes on the plots?
#' @param label_column Use this gene annotation column to pick up gene labels.
#' @param format_sig Use this many significant digits for printing
#'  wacky numbers.
#' @param excel Filename for the excel workbook, or null if not
#'  printed.
#' @param plot_columns A guesstimate of how wide plots are with
#'  respect to 'normally' sized columns in excel.
#' @param alpha Use the alpha channel with this transparency when
#'  plotting.
#' @param z Use this z-score for defining significant in coefficient
#'  plots.
#' @param z_lines Add z-score lines to coefficient plots?
#' @param wanted_genes Character vector of genes to explicitly keep.
#' @param scale_p Scale the p-value in order to allow direct comparisons to dream.
#' @param ... Extra arguments passed to plotters etc.
#' @return Table combining limma/edger/deseq outputs.
#' @seealso [all_pairwise()] [extract_significant_genes()]
#' @examples
#' \dontrun{
#'  expt <- create_expt(metadata="some_metadata.xlsx", gene_info=funkytown)
#'  big_result <- all_pairwise(expt, model_batch=FALSE)
#'  pretty <- combine_de_tables(big_result, table='t12_vs_t0')
#'  pretty <- combine_de_tables(big_result, table='t12_vs_t0',
#'                              keepers=list("avsb"=c("a","b")))
#'  pretty <- combine_de_tables(big_result, table='t12_vs_t0',
#'                              keepers=list("avsb"=c("a","b")),
#'                              excludes=list("description"=c("sno","rRNA")))
#' }
#' @export
combine_de_tables <- function(apr, extra_annot = NULL, keepers = "all", excludes = NULL,
                              adjp = TRUE, includes = NULL,
                              add_plots = TRUE, loess = FALSE, plot_dim = 6,
                              compare_plots = TRUE, padj_type = "fdr", fancy = FALSE,
                              lfc_cutoff = 1.0, p_cutoff = 0.05,
                              excel_title = "Table SXXX: Combined Differential Expression of YYY",
                              increment_start = "SXXX", start_worksheet_num = 2,
                              rda = NULL, rda_input = FALSE, label = 10, label_column = "hgnc_symbol",
                              format_sig = 4, excel = NULL, plot_columns = 10,
                              alpha = 0.4, z = 1.5, z_lines = TRUE, wanted_genes = NULL,
                              scale_p = FALSE, ...) {
  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]
  do_excel <- TRUE
  if (is.null(wb)) {
    do_excel <- FALSE
  }
  ## hey, id10t, keep this function right at the top, if you move it later
  ## badness almost certainly will ensue when you perform weird, non-standard
  ## contrasts like mixed models via lmer or interaction models.
  includes <- check_includes(apr, includes = includes)

  plot_colors <- get_colors_by_condition(apr[["input"]])
  ## Create a list of image files so that they may be properly cleaned up
  ## after writing the xlsx file.
  image_files <- c()
  possible_methods <- c("limma", "deseq", "edger", "ebseq", "basic", "noiseq", "dream")
  model_fstrings <- list()
  for (m in seq_along(possible_methods)) {
    method <- possible_methods[m]
    if (isTRUE(includes[[method]])) {
      if (!is.null(apr[[method]][["model_fstring"]])) {
        model_fstrings[[method]] <- apr[[method]][["model_fstring"]]
      }
    } else {
      apr[[method]] <- list()
    }
  }

  model_used <- NULL
  ## I want to be able to print the model
  ## When summarizing the result, I think in most cases
  ## checking for a deseq or limma model should suffice.
  if (!is.null(apr[["deseq"]][["model_fstring"]])) {
    model_used <- apr[["deseq"]][["model_fstring"]]
  } else if (!is.null(apr[["limma"]][["model_fstring"]])) {
    model_used <- apr[["limma"]][["model_fstring"]]
  }
  ## A common request is to have the annotation data added to the table.  Do that here.
  annot_df <- rowData(apr[["input"]])
  if (!is.null(extra_annot)) {
    annot_df <- merge(annot_df, extra_annot, by = "row.names", all.x = TRUE)
    rownames(annot_df) <- annot_df[["Row.names"]]
    annot_df[["Row.names"]] <- NULL
  }

  ## Write the legend.
  legend <- write_combined_legend(
    wb, excel_basename, plot_dim, apr,
    includes, padj_type, fancy = fancy)
  image_files <- c(legend[["image_files"]], image_files)
  table_names <- legend[["table_names"]]

  ## Because we sometimes use different names for tables in the keepers
  ## This is probably not the correct way to handle the list of coefficients.
  ## Instead, we should just grab the conditions list from deseq/edger/etc

  ## FIXME, I am just grabbing the tablenames from the first type observed
  ## This should be improved to get coefficients on a per-method basis
  ## at some point. But I am assuming for the moment that all coefficients
  ## are available in each method's result
  all_coefficients <- unlist(strsplit(x = table_names[[1]], split = "_vs_"))
  ## But I want to see what happend before I change this.
  ## all_coefficients <- apr[["deseq"]][["conditions"]]

  ## Extract and write the data tables.
  extracted <- list(
    "data" = list(),
    "table_names" = list(),
    "plots" = list(),
    "summaries" = data.frame(),
    "numerators" = c(),
    "denominators" = c())
  ## Here, we will look for only those elements in the keepers list.
  ## In addition, if someone wanted a_vs_b, but we did b_vs_a, then this will
  ## flip the logFCs.
  extracted <- extract_keepers(
    extracted, keepers, table_names, all_coefficients, apr,
    adjp = adjp, annot_df = annot_df, includes = includes,
    excludes = excludes, padj_type = padj_type, fancy = fancy, loess = loess,
    lfc_cutoff = lfc_cutoff, p_cutoff = p_cutoff,
    format_sig = format_sig,
    plot_colors = plot_colors, z = z, alpha = alpha, z_lines = z_lines,
    label = label, label_column = label_column, wanted_genes = wanted_genes,
    scale_p = scale_p)
  numerators <- extracted[["numerators"]]
  denominators <- extracted[["denominators"]]

  ## At this point, we have done everything we can to combine the requested
  ## tables. So lets dump the tables to the excel file and compare how the
  ## various tools performed with some venn diagrams, and finally dump the plots
  ## from above into the sheet.
  comp <- list()
  ## The following if() is too long and should be split into its own function.
  if (isTRUE(do_excel)) {
    ## Starting a new counter of sheets.
    ## I am considering adding some logic to collect the linear models
    ## Then check to see if the slopes/intercepts are duplicated across any
    ## of the contrasts, if this is true, then it is highly likely a mistake was made
    ## when setting up the contrasts such that something got duplicated.
    worksheet_number <- start_worksheet_num
    tnames <- names(extracted[["table_names"]])
    ## tsources <- as.character(extracted[["table_names"]])
    for (x in seq_along(tnames)) {
      tab <- tnames[x]
      numerator <- numerators[x]
      numerator_name <- as.character(numerator[[1]][1])
      denominator <- denominators[x]
      denominator_name <- as.character(denominator[[1]][1])
      written_table <- extracted[["data"]][[tab]]
      if (! tabularp(written_table)) {
        message("There is no data for ", tab, ", skipping it.")
        next
      }
      final_excel_title <- gsub(pattern = "YYY", replacement = tab, x = excel_title)
      final_excel_title <- paste0(final_excel_title, "; ", tab, ": Contrast numerator: ",
                                  numerator_name, ".", " Contrast denominator: ",
                                  denominator_name, ".")
      ## Replace the excel table name with the incrementing value
      sheet_increment_string <- glue("S{worksheet_number}")
      worksheet_number <- worksheet_number + 1
      final_excel_title <- gsub(pattern = increment_start,
                                replacement = sheet_increment_string,
                                x = final_excel_title)
      ## Dump each table to the appropriate excel sheet
      xls_result <- write_xlsx(data = written_table, wb = wb, sheet = tab,
                               title = final_excel_title, rownames = TRUE)
      ## The function write_xlsx has some logic in it to get around excel name
      ## limitations (30 characters), therefore set the sheetname to what was
      ## returned in case it had to change the sheet's name.
      sheetname <- xls_result[["sheet"]]
      current_row <- 1
      current_column <- xls_result[["end_col"]] + 2
      if (isTRUE(add_plots)) {
        ## Text on row 1, plots from 2-17 (15 rows)
        mesg("Adding venn plots for ", tnames[x], ".")
        venn_info <- write_venns_de_xlsx(
          written_table, tab, wb, sheetname, current_row, current_column, excel_basename,
          plot_dim, image_files, lfc_cutoff = lfc_cutoff, p_cutoff = p_cutoff,
          includes = includes, plot_columns = plot_columns)
        image_files <- venn_info[["image_files"]]
        wb <- venn_info[["wb"]]
        if (venn_info[["current_row"]]) {
          current_row <- venn_info[["current_row"]]
          current_column <- venn_info[["current_column"]]
        }
        de_plots_written <- write_plots_de_xlsx(includes, extracted, sheetname,
                                                current_row, current_column, tab,
                                                xls_result, wb, plot_dim,
                                                excel_basename, image_files)
        image_files <- de_plots_written[["image_files"]]
        wb <- de_plots_written[["wb"]]
        current_row <- de_plots_written[["current_row"]]
        current_column <- de_plots_written[["current_column"]]
      }
    }
    xls_result <- NULL
    if (nrow(extracted[["data"]][[1]]) < 50) {
      message("The result table is too small for meaningful comparisons.")
      message("The first table has only: ", nrow(extracted[["data"]][[1]]), ".")
    } else {
      ## Now add some summary data and some plots comparing the tools.
      xls_result <- write_combined_summary(wb, excel_basename, apr, extracted, compare_plots,
                                         lfc_cutoff = lfc_cutoff, p_cutoff = p_cutoff)
      image_files <- c(image_files, xls_result[["image_files"]])
    }

    ## Finish up.
    if (!is.null(apr[["original_pvalues"]])) {
      mesg("Appending a data frame of the original pvalues before sva messed with them.")
      xls_result <- write_xlsx(
        wb, data = apr[["original_pvalues"]], sheet = "original_pvalues",
        title = "Original pvalues for all contrasts before sva adjustment.",
        start_row = 1, rownames = TRUE)
    }

    design_result <- write_sample_design(wb, apr)

    mesg("Writing ", excel, ".")
    save_result <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
    if (class(save_result)[1] == "try-error") {
      warning("Saving xlsx failed.")
    }
  } ## End if !is.null(excel)

  ## Finished!  Dump the important stuff into a return list.
  ret <- list(
    "input" = apr,
    "data" = extracted[["data"]],
    "table_names" = extracted[["table_names"]],
    "plots" = extracted[["plots"]],
    "comp_plot" = comp,
    "keepers" = keepers,
    ## Kept is currently broken.
    "kept" = extracted[["kept"]],
    "model_used" = model_used,
    "model_fstrings" = model_fstrings,
    "de_summary" = extracted[["summaries"]])
  class(ret) <- c("combined_de", "list")

  if (!is.null(rda)) {
    varname <- gsub(x = basename(rda), pattern = "\\.rda$", replacement = "")
    tmp <- ret
    if (!isTRUE(rda_input)) {
      mesg("Removing the input and plots from the rda result.")
      tmp[["plots"]] <- NULL
      tmp[["input"]] <- NULL
    }
    assigned <- assign(varname, tmp)
    removed <- rm(list = "tmp")
    ## When saving a rda of the combined tables, it is less likely that one wants a copy of the
    ## entire differential expression analysis produced by all_pairwise().
    mesg("Saving de result as ", varname, " to ", rda, ".")
    saved <- save(list = varname, file = rda, compress = "xz")
    removed <- rm(varname)
  }
  ## Cleanup the saved image files.
  image_dir <- ""
  for (img in image_files) {
    removed <- try(suppressWarnings(file.remove(img)), silent = TRUE)
    image_dir <- dirname(img)
  }
  nodir <- try(unlink(image_dir), silent = TRUE)

  return(ret)
}

#' Print a combined differential expression analysis.
#'
#' @param x List containing the dataframes for each contrast, the
#'  various plots, the set of wanted contrasts, models used, and
#'  summaries of the data.
#' @param ... Other args to match the generic.
#' @export
print.combined_de <- function(x, ...) {
  message("A set of combined differential expression results.")
  summary_table <- x[["de_summary"]][, c("table", "deseq_sigup", "deseq_sigdown",
                                         "edger_sigup", "edger_sigdown",
                                         "limma_sigup", "limma_sigdown")]
  print(summary_table)
  print(upsetr_combined_de(x))
  return(invisible(x))
}

#' Gather data required to make MA/Volcano plots for pairwise comparisons.
#'
#' It should be possible to significantly simplify the arguments passed to this
#' function, but I have thus far focused only on getting it to work with the
#' newly split-apart combine_de_tables() functions.
#'
#' @param name Name of the table to plot.
#' @param combined Modified pairwise result, containing the various DE methods.
#' @param denominator Name of the denominator coefficient.
#' @param numerator Name of the numerator coefficient.
#' @param plot_inputs The individual outputs from limma etc.
#' @param loess Add a loess estimation?
#' @param lfc_cutoff For Volcano/MA plot lines.
#' @param p_cutoff For Volcano/MA plot lines.
#' @param found_table Table name when there are multiple possible choices.
#' @param p_type Use this/these methods' p-value for determining
#'  significance.
#' @param plot_colors Use these colors for numerators/denominators.
#' @param fancy Include fancy pdf/svg versions of plots for publication?
#' @param adjp Use adjusted p-values?
#' @param do_inverse Flip the numerator/denominator?
#' @param invert_colors Conversely, keep the values the same, but flip
#'  the colors.  I think these invert parameters are not needed anymore.
#' @param z Use a z-score cutoff for coefficient plots.
#' @param alpha Add some transparency to the plots.
#' @param z_lines Add lines for zscore cutoffs?
#' @param label Label this number of the top genes.
#' @param label_column Label the top genes with this column.
combine_extracted_plots <- function(name, combined, denominator, numerator, plot_inputs,
                                    loess = FALSE, lfc_cutoff = 1, p_cutoff = 0.05, found_table = NULL,
                                    p_type = "all", plot_colors = NULL, fancy = FALSE,
                                    adjp = TRUE, do_inverse = FALSE, invert_colors = FALSE,
                                    z = 1.5, alpha = 0.4, z_lines = FALSE,
                                    label = 10, label_column = "hgnc_symbol") {
  combined_data <- combined[["data"]]
  plots <- list()
  types <- c()
  for (i in seq_along(plot_inputs)) {
    type <- names(plot_inputs)[i]
    scatter_name <- glue("{type}_scatter_plots")
    ma_name <- glue("{type}_ma_plots")
    vol_name <- glue("{type}_vol_plots")
    p_name <- glue("{type}_p_plots")
    adjp_name <- glue("{type}_adjp_plots")
    plots[[scatter_name]] <- list()
    plots[[ma_name]] <- list()
    plots[[vol_name]] <- list()
    plots[[p_name]] <- list()
    plots[[adjp_name]] <- list()
    types <- c(type, types)
    ## I think I should plot the coefficient plot with the
    ## MA/volcano.  I am stealing the following 6 lines from the
    ## volcano/MA plotter to make pretty colors
    ## A quick short circuit for extra contrasts
    ## They will not have colors defined in the conditions of the input data and so will be NA
    ## here, so drop out now.
    color_high <- plot_colors[numerator]
    color_low <- plot_colors[denominator]
    if (is.na(color_low) || is.na(color_high)) {
      warning("I think this is an extra contrast table, the plots may be weird.")
      color_high <- "darkred"
      color_low <- "darkblue"
    }

    ## The invert parameter only flips the volcano x-axis
    ma_vol_cof <- list()
    ma_vol_coef <- try(extract_de_plots(
      plot_inputs, combined, type = type,
      invert = FALSE, invert_colors = invert_colors,
      numerator = numerator, denominator = denominator, alpha = alpha, z = z,
      lfc_cutoff = lfc_cutoff, p_cutoff = p_cutoff, adjp = adjp, found_table = found_table,
      p_type = p_type, color_low = color_low, color_high = color_high,
      z_lines = z_lines, label = label, label_column = label_column))
    ## FIXME: There appears to be an error in how I do the t-statistic for basic.
    ## Figure that out before re-enabling basic volcano plots.
    if (type == "basic") {
      plots[[vol_name]] <- NULL
    }
    if ("try-error" %in% class(ma_vol_coef)) {
      plots[[scatter_name]] <- NULL
      plots[[vol_name]] <- NULL
      plots[[ma_name]] <- NULL
    } else {
      if ("try-error" %in% class(ma_vol_coef[["coef"]])) {
        plots[[scatter_name]] <- NULL
      } else {
        plots[[scatter_name]] <- ma_vol_coef[["coef"]]
      }
      if ("try-error" %in% class(ma_vol_coef[["volcano"]])) {
        plots[[vol_name]] <- NULL
      } else {
        plots[[vol_name]] <- ma_vol_coef[["volcano"]][["plot"]]
      }
      if ("try-error" %in% class(ma_vol_coef[["ma"]])) {
        plots[[ma_name]] <- NULL
      } else {
        plots[[ma_name]] <- ma_vol_coef[["ma"]][["plot"]]
      }
    }

    if (is.null(p_name)) {
      mesg("Skipping p-value plot for ", t, ".")
    } else {
      ## If one sets the padj_type, then one will need to pull the correct columns
      ## from the data at this point.  In my vignette, I set padj_type to 'BH'
      ## and as a result I have a series of columns: 'limma_adj_bh' etc.
      ## Therefore we need to get that information to this function call.
      pval_plot <- try(plot_de_pvals(combined[["data"]], type = type, p_type = "raw"),
                       silent = TRUE)
      adjpval_plot <- try(plot_de_pvals(combined[["data"]], type = type, p_type = "adj"),
                       silent = TRUE)
      if ("try-error" %in% class(pval_plot)) {
        plots[[p_name]] <- NULL
        plots[[adjp_name]] <- NULL
      } else {
        plots[[p_name]] <- pval_plot
        plots[[adjp_name]] <- adjpval_plot
      }
    }
  }
  return(plots)
}

get_num_den <- function(string, split = "_vs_") {
  num_den_names <- strsplit(x = string, split = split)[[1]]
  num_name <- num_den_names[[1]]
  den_name <- num_den_names[[2]]
  ret <- list(
    "numerator" = num_name,
    "denominator" = den_name)
  return(ret)
}

check_single_de_table <- function(pairwise, table_name, wanted_numerator,
                                  wanted_denominator, extra = FALSE, type = "deseq") {
  num_den_names <- get_num_den(table_name)
  num_name <- num_den_names[["numerator"]]
  den_name <- num_den_names[["denominator"]]
  inverse_name <- paste0(den_name, "_vs_", num_name)
  forward_test_table <- pairwise[["all_tables"]][[table_name]]
  reverse_test_table <- pairwise[["all_tables"]][[inverse_name]]

  ret <- list(
    "table" = NULL,
    "include" = TRUE,
    "data_orientation" = "forward",
    "tablename_orientation" = "forward")
  if (is.null(pairwise) || class(pairwise)[1] == "try-error") {
    ## mesg("The ", type, " table for ", table_name, " is null.")
    ret[["include"]] <- FALSE
    return(ret)
  }

  if (wanted_numerator == num_name && wanted_denominator == den_name) {
    ret[["data_orientation"]] <- "forward"
  } else if (wanted_numerator == den_name && wanted_denominator == num_name) {
    ret[["data_orientation"]] <- "reverse"
  } else {
    stop("Something is terribly wrong when checking the orientation of tables.")
  }

  if (!is.null(forward_test_table)) {
    ret[["table"]] <- forward_test_table
  } else if (!is.null(reverse_test_table)) {
    ret[["table"]] <- reverse_test_table
    ret[["tablename_orientation"]] <- "reverse"
  }

  return(ret)
}

#' Combine data taken from map_keepers() into a single large table.
#'
#' This is part of an ongoing attempt to simplify and clean up the
#' combine_de_tables() function.  I am hoping that map_keepers and
#' this will be able to take over all the logic currently held in the
#' various extract_keepers_xxx() functions.
#'
#' @param entry Single entry from map_keepers() which provides
#'  orientation information about the table from all_pairwise(), along
#'  with the actual data.
#' @param includes List of methods to include.
#' @param adjp Used adjusted pvalues when defining 'significant.?
#' @param padj_type Perform this type of pvalue adjustment.
#' @param annot_df Include these annotations in the result tables.
#' @param excludes When provided as a list, remove any rows with values in the column defined by the list names, otherwise exclude rownames.
#' @param lfc_cutoff Use this value for a log2FC significance cutoff.
#' @param p_cutoff Use this value for a(n adjusted) pvalue
#'  significance cutoff.
#' @param format_sig Use this many significant digits for some of the
#'  unwieldy numbers.
#' @param sheet_count Start with these sheet number and increment for excel.
#' @param keep_underscore Sanitize underscores?
#' @param wanted_genes Vector of explicitly wanted genes.
#' @param scale_p Scale p-values to provide explicit comparison to dream.
combine_mapped_table <- function(entry, includes, adjp = TRUE, padj_type = "fdr",
                                 annot_df = NULL, excludes = NULL, lfc_cutoff = 1,
                                 p_cutoff = 0.05, format_sig = 4, sheet_count = 0,
                                 keep_underscore = TRUE, wanted_genes = NULL,
                                 scale_p = FALSE) {
  if (padj_type[1] != "ihw" && !(padj_type %in% p.adjust.methods)) {
    warning("The p adjustment ", padj_type, " is not in the set of p.adjust.methods.
Defaulting to fdr.")
    padj_type <- "fdr"
  }

  badf <- data.frame("basic_num" = 0, "basic_den" = 0, "basic_numvar" = 0,
                     "basic_denvar" = 0, "basic_logfc" = 0, "basic_t" = 0,
                     "basic_p" = 0, "basic_adjp" = 0)
  rownames(badf) <- "dropme"

  dedf <- data.frame("deseq_basemean" = 0, "deseq_logfc" = 0, "deseq_lfcse" = 0,
                     "deseq_stat" = 0, "deseq_p" = 0, "deseq_adjp" = 0,
                     "deseq_num" = 0, "deseq_den" = 0)
  rownames(dedf) <- "dropme"

  drdf <- data.frame("dream_logfc" = 0, "dream_ave" = 0, "dream_t" = 0,
                     "dream_p" = 0, "dream_adjp" = 0, "dream_b" = 0)
  rownames(drdf) <- "dropme"
  if (!is.null(entry[["dream"]][["data"]][["z.std"]])) {
    drdf[["dream_zstd"]] <- 0
  }

  ebdf <- data.frame("ebseq_logfc" = 0, "ebseq_c1mean" = 0,
                     "ebseq_c2mean" = 0, "ebseq_mean" = 0, "ebseq_var" = 0,
                     "ebseq_postfc" = 0, "ebseq_ppde" = 0, "ebseq_adjp" = 0)
  rownames(ebdf) <- "dropme"

  eddf <- data.frame("edger_logfc" = 0, "edger_logcpm" = 0, "edger_lr" = 0,
                     "edger_p" = 0, "edger_adjp" = 0)
  rownames(eddf) <- "dropme"

  lidf <- data.frame("limma_logfc" = 0, "limma_ave" = 0, "limma_t" = 0,
                     "limma_p" = 0, "limma_adjp" = 0, "limma_b" = 0)
  rownames(lidf) <- "dropme"

  nodf <- data.frame("noiseq_num" = 0, "noiseq_den" = 0, "noiseq_theta" = 0,
                     "noiseq_prob" = 0, "noiseq_logfc" = 0, "noiseq_p" = 0, "noiseq_adjp" = 0)
  rownames(nodf) <- "dropme"


  ## I am changing the logic of this function so that the inversion of the values
  ## is no longer connected to the inversion of colors
  invert_colors <- c()
  inverts <- list(
    "basic" = FALSE,
    "deseq" = FALSE,
    "dream" = FALSE,
    "ebseq" = FALSE,
    "edger" = FALSE,
    "limma" = FALSE,
    "noiseq" = FALSE)
  for (query in names(includes)) {
    if (!isTRUE(includes[[query]])) {
      next
    }
    if (is.null(entry[[query]])) {
      includes[[query]] <- FALSE
    } else {
      if (entry[[query]][["data_orientation"]] == "reverse" &&
        entry[[query]][["tablename_orientation"]] == "expected") {
        inverts[[query]] <- TRUE
      } else if (entry[[query]][["data_orientation"]] == "forward" &&
        entry[[query]][["tablename_orientation"]] == "unexpected") {
        inverts[[query]] <- TRUE
      }

      ## These little stanzas are a little redundant because I specify the column names above,
      ## but it does provide me a chance to play with the names in case I am using different
      ## metrics (like mean vs median).
      if (query == "basic") {
        if (is.null(entry[[query]][["data"]])) {
          mesg("Method ", query, " does not have this contrast.")
          badf <- data.frame()
          ba_stats <- data.frame()
          ba_lfc_adjp <- data.frame()
        } else {
          badf <- entry[[query]][["data"]]
          ## I recently changed basic to optionally do means or medians.  I need to take that into
          ## account when working with these tables.  For the moment, I think I will simply rename
          ## the column to _median to avoid confusion.
          colnames(badf) <- c("basic_num", "basic_den", "basic_numvar",
                              "basic_denvar", "basic_t", "basic_p",
                              "basic_logfc", "basic_adjp")
          ba_stats <- badf[, c("basic_num", "basic_den", "basic_numvar", "basic_denvar",
                               "basic_t", "basic_p")]
          ba_lfc_adjp <- badf[, c("basic_logfc", "basic_adjp")]
        }
      }
      if (query == "deseq") {
        if (is.null(entry[[query]][["data"]])) {
          mesg("Method ", query, " does not have this contrast.")
          dedf <- data.frame()
          de_stats <- data.frame()
          de_lfc_adjp <- data.frame()
        } else {
          dedf <- entry[[query]][["data"]]
          colnames(dedf) <- c("deseq_basemean", "deseq_logfc", "deseq_lfcse",
            "deseq_stat", "deseq_p", "deseq_adjp",
            "deseq_num", "deseq_den")
          de_stats <- dedf[, c("deseq_basemean", "deseq_lfcse", "deseq_stat",
            "deseq_p", "deseq_num", "deseq_den")]
          de_lfc_adjp <- dedf[, c("deseq_logfc", "deseq_adjp")]
        }
      }
      if (query == "dream") {
        if (is.null(entry[[query]][["data"]])) {
          mesg("dream does not have this entry.")
          drdf <- data.frame()
          dr_stats <- data.frame()
          dr_lfc_adjp <- data.frame
        } else {
          drdf <- entry[[query]][["data"]]
          if (ncol(drdf) == 6) {
            colnames(drdf) <- c("dream_logfc", "dream_ave", "dream_t", "dream_p",
                                "dream_adjp", "dream_b")
            dr_stats <- drdf[, c("dream_ave", "dream_t", "dream_p", "dream_b")]
          } else {
            colnames(drdf) <- c("dream_logfc", "dream_ave", "dream_t", "dream_p",
                                "dream_adjp", "dream_b", "dream_zstd")
            dr_stats <- drdf[, c("dream_ave", "dream_t", "dream_p", "dream_p", "dream_zstd")]

          }
          dr_lfc_adjp <- drdf[, c("dream_logfc", "dream_adjp")]
        }
      }
      if (query == "ebseq") {
        if (is.null(entry[[query]][["data"]])) {
          mesg("EBSeq does not have this entry.")
          ebdf <- data.frame()
          eb_stats <- data.frame()
          eb_lfc_adjp <- data.frame()
        } else {
          ebdf <- entry[[query]][["data"]]
          colnames(ebdf) <- c("ebseq_logfc", "ebseq_c1mean",
                              "ebseq_c2mean", "ebseq_mean", "ebseq_var",
                              "ebseq_postfc", "ebseq_ppde", "ebseq_adjp")
          eb_stats <- ebdf[, c("ebseq_c1mean", "ebseq_c2mean",
                               "ebseq_mean", "ebseq_postfc",
                               "ebseq_ppde")]
          eb_lfc_adjp <- ebdf[, c("ebseq_logfc", "ebseq_adjp")]
        }
      }
      if (query == "edger") {
        if (is.null(entry[[query]][["data"]])) {
          mesg("Edger does not have this entry.")
          eddf <- data.frame()
          ed_stats <- data.frame()
          ed_lfc_adjp <- data.frame()
        } else {
          eddf <- entry[[query]][["data"]]
          colnames(eddf) <- c("edger_logfc", "edger_logcpm", "edger_lr", "edger_p", "edger_adjp")
          ed_stats <- eddf[, c("edger_logcpm", "edger_lr", "edger_p")]
          ed_lfc_adjp <- eddf[, c("edger_logfc", "edger_adjp")]
        }
      }
      if (query == "limma") {
        if (is.null(entry[[query]][["data"]])) {
          mesg("limma does not have this entry.")
          lidf <- data.frame()
          li_stats <- data.frame()
          li_lfc_adjp <- data.frame
        } else {
          lidf <- entry[[query]][["data"]]
          colnames(lidf) <- c("limma_logfc", "limma_ave", "limma_t", "limma_p",
            "limma_adjp", "limma_b")
          li_stats <- lidf[, c("limma_ave", "limma_t", "limma_p", "limma_b")]
          li_lfc_adjp <- lidf[, c("limma_logfc", "limma_adjp")]
        }
      }
      if (query == "noiseq") {
        if (is.null(entry[[query]][["data"]])) {
          mesg("noiseq does not have this entry.")
          nodf <- data.frame()
          no_stats <- data.frame()
          no_lfc_adjp <- data.frame()
        } else {
          nodf <- entry[[query]][["data"]]
          colnames(nodf) <- c("noiseq_num", "noiseq_den", "noiseq_theta", "noiseq_prob",
                              "noiseq_logfc", "noiseq_p", "noiseq_adjp")
          nodf[["noiseq_mean"]] <- (nodf[["noiseq_num"]] + nodf[["noiseq_den"]]) / 2
          no_stats <- nodf[, c("noiseq_num", "noiseq_den", "noiseq_mean",
                               "noiseq_theta", "noiseq_prob", "noiseq_p")]
          no_lfc_adjp <- nodf[, c("noiseq_logfc", "noiseq_adjp")]
        }
      }
    }

    if (isTRUE(inverts[[query]])) {
      invert_colors <- c(invert_colors, query)
    }
  }

  datalst <- list()
  statslst <- list()
  if (isTRUE(includes[["basic"]])) {
    datalst[["basic"]] <- as.data.table(ba_lfc_adjp)
    datalst[["basic"]][["rownames"]] <- rownames(ba_lfc_adjp)
    statslst[["basic"]] <- as.data.table(ba_stats)
    statslst[["basic"]][["rownames"]] <- rownames(ba_stats)
  }
  if (isTRUE(includes[["deseq"]])) {
    datalst[["deseq"]] <- as.data.table(de_lfc_adjp)
    datalst[["deseq"]][["rownames"]] <- rownames(de_lfc_adjp)
    statslst[["deseq"]] <- as.data.table(de_stats)
    statslst[["deseq"]][["rownames"]] <- rownames(de_stats)
  }
  if (isTRUE(includes[["dream"]])) {
    datalst[["dream"]] <- as.data.table(dr_lfc_adjp)
    datalst[["dream"]][["rownames"]] <- rownames(dr_lfc_adjp)
    statslst[["dream"]] <- as.data.table(dr_stats)
    statslst[["dream"]][["rownames"]] <- rownames(dr_stats)
  }
  if (isTRUE(includes[["ebseq"]])) {
    datalst[["ebseq"]] <- as.data.table(eb_lfc_adjp)
    datalst[["ebseq"]][["rownames"]] <- rownames(eb_lfc_adjp)
    statslst[["ebseq"]] <- as.data.table(eb_stats)
    statslst[["ebseq"]][["rownames"]] <- rownames(eb_stats)
  }
  if (isTRUE(includes[["edger"]])) {
    datalst[["edger"]] <- as.data.table(ed_lfc_adjp)
    datalst[["edger"]][["rownames"]] <- rownames(ed_lfc_adjp)
    statslst[["edger"]] <- as.data.table(ed_stats)
    statslst[["edger"]][["rownames"]] <- rownames(ed_stats)
  }
  if (isTRUE(includes[["limma"]])) {
    datalst[["limma"]] <- as.data.table(li_lfc_adjp)
    datalst[["limma"]][["rownames"]] <- rownames(li_lfc_adjp)
    statslst[["limma"]] <- as.data.table(li_stats)
    statslst[["limma"]][["rownames"]] <- rownames(li_stats)
  }
  if (isTRUE(includes[["noiseq"]])) {
    datalst[["noiseq"]] <- as.data.table(no_lfc_adjp)
    datalst[["noiseq"]][["rownames"]] <- rownames(no_lfc_adjp)
    statslst[["noiseq"]] <- as.data.table(no_stats)
    statslst[["noiseq"]][["rownames"]] <- rownames(no_stats)
  }

  ## Make the initial data structure
  wanted <- names(datalst)[[1]]
  num_stats <- length(statslst)
  num_data <- length(datalst)
  if (num_data == 1) {
    if (num_stats == 1) {
      ## Then this is a chunk of data and associated stats.
      comb <- merge(datalst[[1]], statslst[[1]], by = "rownames", all = TRUE)
    } else {
      ## Then there must only be a chunk of data.
      comb <- datalst[[1]]
    }
  } else {
    ## There is more than one set of data to merge.
    comb <- datalst[[1]]
    for (i in seq(from = 2, to = length(datalst))) {
      comb <- merge(comb, datalst[[i]], by = "rownames", all = TRUE)
    }
    if (length(statslst) > 0) {
      for (j in seq_along(statslst)) {
        comb <- merge(comb, statslst[[j]], by = "rownames", all = TRUE)
      }
    }
  }
  ## Doing the merge in the way above will lead to a single row which is essentially blank
  ## It is probably the first row.

  ## The next lines are intended to drop that blank row.
  comb <- as.data.frame(comb)
  rownames(comb) <- comb[["rownames"]]
  dropme <- rownames(comb) == "dropme"
  comb <- comb[!dropme, ]

  keepers <- colnames(comb) != "rownames"
  comb <- comb[, keepers, drop = FALSE]
  comb[is.na(comb)] <- 0

  ## Invert logFC when needed.
  for (query in names(inverts)) {
    fc_column <- paste0(query, "_logfc")
    if (isTRUE(inverts[[query]])) {
      mesg("Inverting column: ", fc_column, " for ", query, ".")
      comb[[fc_column]] <- comb[[fc_column]] * -1.0
      if (query == "deseq") {
        comb[["deseq_stat"]] <- comb[["deseq_stat"]] * -1.0
        tmp <- comb[["deseq_num"]]
        comb[["deseq_num"]] <- comb[["deseq_den"]]
        comb[["deseq_den"]] <- tmp
      } else if (query == "noiseq") {
        tmp <- comb[["noiseq_num"]]
        comb[["noiseq_num"]] <- comb[["noiseq_den"]]
        comb[["noiseq_den"]] <- tmp
      }
    }
  }

  ## Add one final p-adjustment to ensure a consistent and user defined value.
  mesg(" Performing ", padj_type, " p-adjustment across all methods.")
  extra_adjust_columns <- c()
  if (!is.null(comb[["limma_p"]])) {
    colname <- glue("limma_adjp_{padj_type}")
    extra_adjust_columns <- c(extra_adjust_columns, colname)
    comb[[colname]] <- hpgl_padjust(comb, pvalue_column = "limma_p",
                                    mean_column = "limma_ave",
                                    method = padj_type, significance = p_cutoff)
    if (isTRUE(scale_p)) {
      comb[["limma_p_zstd"]] <- scale(comb[["limma_p"]], center = TRUE, scale = TRUE)
      extra_adjust_columns <- c(extra_adjust_columns, "limma_p_zstd")
    }
  }
  if (!is.null(comb[["dream_p"]])) {
    colname <- glue("dream_adjp_{padj_type}")
    extra_adjust_columns <- c(extra_adjust_columns, colname)
    comb[[colname]] <- hpgl_padjust(comb, pvalue_column = "dream_p",
                                    mean_column = "dream_ave",
                                    method = padj_type, significance = p_cutoff)
    if (isTRUE(scale_p)) {
      comb[["dream_p_zstd"]] <- scale(comb[["dream_p"]], center = TRUE, scale = TRUE)
      extra_adjust_columns <- c(extra_adjust_columns, "dream_p_zstd")
    }
  }
  if (!is.null(comb[["deseq_p"]])) {
    colname <- glue("deseq_adjp_{padj_type}")
    extra_adjust_columns <- c(extra_adjust_columns, colname)
    comb[[colname]] <- hpgl_padjust(comb, pvalue_column = "deseq_p",
                                    mean_column = "deseq_basemean",
                                    method = padj_type, significance = p_cutoff)
    if (isTRUE(scale_p)) {
      comb[["deseq_p_zstd"]] <- scale(comb[["limma_p"]], center = TRUE, scale = TRUE)
      extra_adjust_columns <- c(extra_adjust_columns, "deseq_p_zstd")
    }
  }
  if (!is.null(comb[["edger_p"]])) {
    colname <- glue("edger_adjp_{padj_type}")
    extra_adjust_columns <- c(extra_adjust_columns, colname)
    comb[[colname]] <- hpgl_padjust(comb, pvalue_column = "edger_p",
                                    mean_column = "edger_logcpm",
                                    method = padj_type, significance = p_cutoff)
    if (isTRUE(scale_p)) {
      comb[["edger_p_zstd"]] <- scale(comb[["limma_p"]], center = TRUE, scale = TRUE)
      extra_adjust_columns <- c(extra_adjust_columns, "edger_p_zstd")
    }
  }
  if (!is.null(comb[["ebseq_ppde"]])) {
    colname <- glue("ebseq_adjp_{padj_type}")
    extra_adjust_columns <- c(extra_adjust_columns, colname)
    comb[[colname]] <- hpgl_padjust(comb, pvalue_column = "ebseq_ppde",
                                    mean_column = "ebseq_mean",
                                    method = padj_type, significance = p_cutoff)
    if (isTRUE(scale_p)) {
      comb[["ebseq_p_zstd"]] <- scale(comb[["ebseq_ppde"]], center = TRUE, scale = TRUE)
      extra_adjust_columns <- c(extra_adjust_columns, "ebseq_p_zstd")
    }
  }
  if (!is.null(comb[["basic_p"]])) {
    colname <- glue("basic_adjp_{padj_type}")
    extra_adjust_columns <- c(extra_adjust_columns, colname)
    comb[[colname]] <- hpgl_padjust(comb, pvalue_column = "basic_p",
                                    mean_column = "basic_num",
                                    method = padj_type, significance = p_cutoff)
    if (isTRUE(scale_p)) {
      comb[["basic_p_zstd"]] <- scale(comb[["basic_p"]], center = TRUE, scale = TRUE)
      extra_adjust_columns <- c(extra_adjust_columns, "basic_p_zstd")
    }
  }
  if (!is.null(comb[["noiseq_p"]])) {
    colname <- glue("noiseq_adjp_{padj_type}")
    extra_adjust_columns <- c(extra_adjust_columns, colname)
    comb[[colname]] <- hpgl_padjust(comb, pvalue_column = "noiseq_p",
                                    mean_column = "noiseq_num",
                                    method = padj_type, significance = p_cutoff)
    if (isTRUE(scale_p)) {
      comb[["noiseq_p_zstd"]] <- scale(comb[["noiseq_p"]], center = TRUE, scale = TRUE)
      extra_adjust_columns <- c(extra_adjust_columns, "noiseq_p_zstd")
    }
  }

  ## Set a consistent numeric format for the adjusted p-values
  if (is.numeric(format_sig)) {
    for (colname in extra_adjust_columns) {
      comb[[colname]] <- as.numeric(format(x = comb[[colname]], digits = format_sig,
        scientific = TRUE, trim = TRUE))
    }
  }

  ## I made an odd choice in a moment to normalize.quantiles the combined fold changes
  ## This should be reevaluated
  temp_fc <- data.frame()
  if (isTRUE(includes[["limma"]]) &&
        isTRUE(includes[["deseq"]]) &&
        isTRUE(includes[["edger"]])) {
    ## Why did I rename this to temp_helpfc?
    ## I think this should be removed.
    temp_helpfc <- cbind(as.numeric(comb[["limma_logfc"]]),
                     as.numeric(comb[["edger_logfc"]]),
                     as.numeric(comb[["deseq_logfc"]]))
    temp_fc <- try(preprocessCore::normalize.quantiles(as.matrix(temp_helpfc)), silent = TRUE)
    ## I just got a strange error: 'vector types do not match in copyVector'
    if (! "try-error" %in% class(temp_fc)) {
      comb[["lfc_meta"]] <- rowMeans(temp_fc, na.rm = TRUE)
      comb[["lfc_var"]] <- genefilter::rowVars(temp_fc, na.rm = TRUE)
      comb[["lfc_varbymed"]] <- comb[["lfc_var"]] / comb[["lfc_meta"]]
      temp_p <- cbind(as.numeric(comb[["limma_p"]]),
                      as.numeric(comb[["edger_p"]]),
                      as.numeric(comb[["deseq_p"]]))
      comb[["p_meta"]] <- rowMeans(temp_p, na.rm = TRUE)
      comb[["p_var"]] <- genefilter::rowVars(temp_p, na.rm = TRUE)
      if (is.numeric(format_sig)) {
        comb[["lfc_meta"]] <- signif(x = comb[["lfc_meta"]], digits = format_sig)
        comb[["lfc_var"]] <- format(x = comb[["lfc_var"]], digits = format_sig,
                                    scientific = TRUE, trim = TRUE)
        comb[["lfc_varbymed"]] <- format(x = comb[["lfc_varbymed"]], digits = format_sig,
                                         scientific = TRUE, trim = TRUE)
        comb[["p_var"]] <- format(x = comb[["p_var"]], digits = format_sig,
                                  scientific = TRUE, trim = TRUE)
        comb[["p_meta"]] <- format(x = comb[["p_meta"]], digits = format_sig,
                                   scientific = TRUE, trim = TRUE)
      }
    }
  }
  if (!is.null(annot_df)) {
    if (isTRUE(keep_underscore)) {
      colnames(annot_df) <- gsub(pattern = "[^_[:^punct:]]",
                                 replacement = "", x = colnames(annot_df), perl = TRUE)
    } else {
      colnames(annot_df) <- gsub(pattern = "[[:punct:]]",
                                 replacement = "", x = colnames(annot_df), perl = TRUE)
    }
    comb <- merge(annot_df, comb, by = "row.names", all.y = TRUE)
    rownames(comb) <- comb[["Row.names"]]
    comb[["Row.names"]] <- NULL
    colnames(comb) <- make.names(tolower(colnames(comb)), unique = TRUE)
  }

  ## Exclude rows based on a list of unwanted columns/strings
  if ("list" %in% class(excludes)) {
    comb <- exclude_genes_from_combined_list(comb, excludes)
  } else if ("character" %in% class(excludes)) {
    keepers <- ! rownames(comb) %in% excludes
    mesg("Keeping ", sum(keepers), " genes out of ", length(rownames(comb)), ".")
    comb <- comb[keepers, ]
  }

  if (!is.null(wanted_genes)) {
    wanted_idx <- rownames(comb) %in% wanted_genes
    mesg("Keeping ", sum(wanted_idx), " rows found in the set of ",
         length(wanted_genes), " wanted genes.")
    comb <- comb[wanted_idx, ]
  }

  up_fc <- lfc_cutoff
  down_fc <- -1.0 * lfc_cutoff
  summary_table_name <- entry[[1]][["string"]]
  if (entry[[1]][["data_orientation"]] == "reverse") {
    summary_table_name <- glue("{summary_table_name}-inverted")
  }
  limma_p_column <- "limma_adjp"
  deseq_p_column <- "deseq_adjp"
  edger_p_column <- "edger_adjp"
  if (!isTRUE(adjp)) {
    limma_p_column <- "limma_p"
    deseq_p_column <- "deseq_p"
    edger_p_column <- "edger_p"
  }
  summary_lst <- summarize_combined(comb, up_fc, down_fc, p_cutoff)
  summary_lst[["table"]] <- summary_table_name
  ret <- list(
    "data" = comb,
    "includes" = includes,
    "inverts" = inverts,
    "summary" = summary_lst)
  class(ret) <- c("combined_table", "list")
  return(ret)
}

#' Print a single combined DE result.
#'
#' @param x Data table of combined differential expression results.
#' @param ... Other args to match the generic.
#' @export
print.combined_table <- function(x, ...) {
  message("A combined differential expression table.")
  message("Comprising ", nrow(x[["comb"]]), " genes and including ", x[["includes"]], ".")
  return(invisible(x))
}

exclude_genes_from_combined_list <- function(comb, excludes) {
  for (colnum in seq_along(excludes)) {
    col <- names(excludes)[colnum]
    for (exclude_num in seq_along(excludes[[col]])) {
      exclude <- excludes[[col]][exclude_num]
      remove_column <- comb[[col]]
      remove_idx <- grep(pattern = exclude, x = remove_column, perl = TRUE, invert = TRUE)
      removed_num <- sum(as.numeric(remove_idx))
      mesg("Removed ", removed_num, " genes using ",
           exclude, " as a string against column ", col, ".")
      comb <- comb[remove_idx, ]
    }  ## End iterating through every string to exclude
  }  ## End iterating through every element of the exclude list
  return(comb)
}

#' Find the correct tables given a set of definitions of desired tables, numerators/denominators.
#'
#' This is responsible for hunting down tables which correspond to the various ways one may
#' represent them.
#'
#' @param keepers List/scalar representation of desired tables.
#' @param table_names The actual list of results produced by the various methods employed.
#' @param datum The full dataset.
map_keepers <- function(keepers, table_names, datum) {
  keeper_table_map <- list()
  ## I changed table_names to be keyed off method because we cannot assume
  ## each method has the same tables at this time
  ## So: keeper_table_map->{method}->{table_name}
  for (m in names(table_names)) {
    keeper_table_map[[m]] <- list()
  }
  numerators <- denominators <- c()

  contrast_list <- names(keepers)
  for (a in seq_along(names(keepers))) {
    name <- names(keepers)[a]
    ## Initially, set same_string to the name of the table, then if there is a
    ## table with that, just use it directly rather than hunt for numerator/denominator.
    same_string <- name
    ## The numerators and denominators will be used to check that we are a_vs_b or b_vs_a
    numerator <- keepers[[name]][1]
    if (is.null(numerator)) {
      ## Then this should be the name of a special table, e.g. from extra_contrasts.
      mesg("Assuming ", name, " is a table from extra_contrasts.")
      inverse_string <- same_string
    } else {
      numerators[name] <- numerator
      denominator <- keepers[[name]][2]
      denominators[name] <- denominator
      same_string <- numerator
      inverse_string <- numerator
      if (!is.na(denominator)) {
        same_string <- glue("{numerator}_vs_{denominator}")
        inverse_string <- glue("{denominator}_vs_{numerator}")
      }
    }

    method_list <- list()
    for (m in seq_along(names(table_names))) {
      method <- names(table_names)[m]
      method_names <- table_names[[method]]
      if (same_string %in% method_names && inverse_string %in% method_names) {
        keeper_table_map[[method]][[name]] <- list(
          "string" = same_string,
          "reverse_string" = inverse_string,
          "data_orientation" = "both")
      } else if (same_string %in% method_names) {
        keeper_table_map[[method]][[name]] <- list(
          "string" = same_string,
          "reverse_string" = inverse_string,
          "data_orientation" = "forward")
      } else if (inverse_string %in% method_names) {
        keeper_table_map[[method]][[name]] <- list(
          "string" = inverse_string,
          "reverse_string" = inverse_string,
          "data_orientation" = "reverse")
      } else {
        keeper_table_map[[method]][[name]] <- list(
          "string" = FALSE,
          "data_orientation" = FALSE)
      }
      keeper_table_map[[method]][[name]][["wanted_numerator"]] <- numerator
      keeper_table_map[[method]][[name]][["wanted_denominator"]] <- denominator
      mesg("Keeper ", name, " is ", numerator, "/", denominator,
           ": ", keeper_table_map[[method]][[name]][["data_orientation"]], " ",
           keeper_table_map[[method]][[name]][["string"]], ".")
      ## Let us check that the individual pairwise contrasts have the same tables in the same order.
      individual_tables <- list()
      if (is.null(datum[[method]])) {
        message("The ", method, " slot of datum is missing.")
      } else {
        contrast_string <- keeper_table_map[[method]][[name]][["string"]]
        mesg("Checking ", method, " for name ", name, ":", contrast_string)
        available_contrasts <- names(datum[[method]][["all_tables"]])
        position <- which(available_contrasts %in% contrast_string)
        test_name <- names(datum[[method]][["all_tables"]])[position]
        ## This index is not there, hopefully it is an extra.
        if (is.na(test_name) || length(test_name) == 0) {
          keeper_table_map[[method]][[name]][["data"]] <- NULL
          keeper_table_map[[method]][[name]][["tablename_orientation"]] <- "none"
        } else if (test_name == keeper_table_map[[method]][[name]][["string"]]) {
          keeper_table_map[[method]][[name]][["data"]] <- datum[[method]][["all_tables"]][[position]]
          keeper_table_map[[method]][[name]][["tablename_orientation"]] <- "expected"
        } else if (test_name == keeper_table_map[[method]][[name]][["reverse_string"]]) {
          keeper_table_map[[method]][[name]][["data"]] <- datum[[method]][["all_tables"]][[position]]
          keeper_table_map[[method]][[name]][["tablename_orientation"]] <- "unexpected"
        } else {
          message("The name did not match the keeper_table_map[[name]].")
          keeper_table_map[[method]][[name]][["data"]] <- NULL
          keeper_table_map[[method]][[name]][["tablename_orientation"]] <- "none"
        }
      }
    }
  }
  class(keeper_table_map) <- "mapped_keepers"
  return(keeper_table_map)
}

#' Print a set of mapped keepers from combine_de_tables()
#'
#' @param x List full of kept information.
#' @param ... Other args to match the generic.
#' @export
print.mapped_keepers <- function(x, ...) {
  message("The set of mappings between the wanted data and available data in a pairwise comparison.")
  print(names(x))
  return(invisible(x))
}

#' When a list of 'keeper' contrasts is specified, extract it from the data.
#'
#' This is the most interesting of the extract_keeper functions.  It must check
#' that the numerators and denominators match the desired contrast and flip the
#' signs in the logFCs when appropriate.
#'
#' @param extracted Tables extracted from the all_pairwise data.
#' @param keepers In this case, one may assume either NULL or 'all'.
#' @param table_names The set of tables produced by all_pairwise().
#' @param all_coefficients The set of all experimental conditions in the
#'  experimental metadata.
#' @param apr The result from all_pairwise(), containing the
#'  limma/edger/deseq/etc data.
#' @param adjp Pull out the adjusted p-values from the data?
#' @param annot_df What annotations should be added to the table?
#' @param includes List of predicates by method.
#' @param excludes Set of genes to exclude.
#' @param padj_type Choose a specific p adjustment.
#' @param min_genes Minimum number of genes below which no plots will
#'  be generated.
#' @param fancy Include larger pdf/svg plots with the xlsx output?
#' @param loess Add a loess to plots?
#' @param lfc_cutoff Passed for volcano/MA plots and defining 'significant'
#' @param p_cutoff Passed for volcano/MA plots and defining 'significant'
#' @param format_sig Number of significant digits for stuff like
#'  pvalues.
#' @param plot_colors Define what colors should be used for
#'  'up'/'down'
#' @param z Define significantly away from the identity line in a
#'  coefficient plot.
#' @param alpha Use this alpha transparency for plots.
#' @param z_lines Include lines denoting significant z-scores?
#' @param label When not NULL, label this many genes.
#' @param label_column Try using this column for labeling genes.
#' @param wanted_genes Genes explicited desired.
#' @param scale_p Scale the pvalues to compare against dream?
#' @return The extracted, but with more stuff at the end!
extract_keepers <- function(extracted, keepers, table_names,
                            all_coefficients, apr,
                            adjp, annot_df, includes,
                            excludes, padj_type, min_genes = 10,
                            fancy = FALSE, loess = FALSE,
                            lfc_cutoff = 1.0, p_cutoff = 0.05,
                            format_sig = 4, plot_colors = NULL,
                            z = 1.5, alpha = 0.4, z_lines = FALSE,
                            label = 10, label_column = "hgnc_symbol",
                            wanted_genes = NULL, scale_p = FALSE) {
  ## First check that your set of keepers is in the data
  ## all_keepers is therefore the set of all coefficients in numerators/denominators
  message("Looking for subscript invalid names, start of extract_keepers.")
  all_keepers <- as.character(unlist(keepers))

  ## keeper_names is the set of arbitrarily chosen names for the written sheets
  keeper_names <- names(keepers)
  found_keepers <- sum(all_keepers %in% all_coefficients)
  ## Add a little logic in case we have weirdo contrasts, in which case
  ## the name of the table should remain constant and the
  ## numerator/denominator will be ignored because they should have
  ## been defined when setting up the weirdo contrast.

  ## The next few lines are looking to find 'extra' tables, these are the set of possible
  ## non-pairwise tables that one might consider e.g. more complex combinations of factors.
  ## If these complex tables are in the data, then I assume the names will exactly match the
  ## table names, otherwise it is too confusing to find them.
  ## This only considers deseq/edger/limma because those are (at least currently)
  ## the 'first class' methods I consider; my basic method, ebseq, noiseq, dream
  ## do not get as much consideration at this time.
  extra_names_edger_idx <- names(apr[["edger"]][["all_tables"]]) %in% keeper_names
  extra_names_deseq_idx <- names(apr[["deseq"]][["all_tables"]]) %in% keeper_names
  extra_names_limma_idx <- names(apr[["limma"]][["all_tables"]]) %in% keeper_names

  extra_names <- c()
  if (sum(extra_names_edger_idx) > 0) {
    extra_names <- names(apr[["edger"]][["all_tables"]])[extra_names_edger_idx]
  }
  if (sum(extra_names_deseq_idx) > 0) {
    extra_names <- c(extra_names,
                     names(apr[["deseq"]][["all_tables"]])[extra_names_deseq_idx])
  }
  if (sum(extra_names_limma_idx) > 0) {
    extra_names <- c(extra_names,
                     names(apr[["limma"]][["all_tables"]])[extra_names_limma_idx])
  }
  extra_names <- unique(extra_names)
  ## Double check that the _vs_ was not removed.
  ## On a few occasions I forgot about the _vs_
  if (is.null(extra_names)) {
    no_vs <- gsub(x = keeper_names, pattern = "_vs_", replacement = "_")
    novs_names_edger <- gsub(x = names(apr[["edger"]][["all_tables"]]),
                             pattern = "_vs_",
                             replacement = "_")
    novs_names_limma <- gsub(x = names(apr[["limma"]][["all_tables"]]),
                             pattern = "_vs_",
                             replacement = "_")
    novs_names_deseq <- gsub(x = names(apr[["deseq"]][["all_tables"]]),
                             pattern = "_vs_",
                             replacement = "_")
    extra_names_edger_idx <- novs_names_edger %in% no_vs
    extra_names_limma_idx <- novs_names_limma %in% no_vs
    extra_names_deseq_idx <- novs_names_deseq %in% no_vs
    if (sum(extra_names_edger_idx) > 0) {
      extra_names <- names(apr[["edger"]][["all_tables"]])[extra_names_edger_idx]
    }
    if (sum(extra_names_deseq_idx) > 0) {
      extra_names <- c(extra_names,
                       names(apr[["deseq"]][["all_tables"]])[extra_names_deseq_idx])
    }
    if (sum(extra_names_limma_idx) > 0) {
      extra_names <- c(extra_names,
                       names(apr[["limma"]][["all_tables"]])[extra_names_limma_idx])
    }
    extra_names <- unique(extra_names)
  }

  ## Just make sure we have something to work with.
  if (sum(found_keepers, length(extra_names)) == 0) {
    message("The keepers has no elements in the coefficients.")
    message("Here are the keepers: ", toString(all_keepers))
    message("Here are the coefficients: ", toString(all_coefficients))
    stop("Unable to find the set of contrasts to keep, fix this and try again.")
  }
  ## I do not think this unique() should be needed.
  ## but it appears that sometimes the above if() statements are causing duplicates.
  ## I changed the table_names to a list by method, primarily because ebseq is currently
  ## generating its tables in a different fashion and so cannot be assumed to have the
  ## same set of names (e.g. some may be reversed compared to others)
  for (type in names(table_names)) {
    table_names[[type]] <- c(table_names[[type]], extra_names)
  }

  datum <- list(
    "basic" = apr[["basic"]],
    "deseq" = apr[["deseq"]],
    "dream" = apr[["dream"]],
    "ebseq" = apr[["ebseq"]],
    "edger" = apr[["edger"]],
    "limma" = apr[["limma"]],
    "noiseq" = apr[["noiseq"]])
  keeper_table_map <- map_keepers(keepers, table_names, datum)
  rekeyed_map <- list()
  for (m in seq_along(names(keeper_table_map))) {
    method <- names(keeper_table_map)[m]
    element <- keeper_table_map[[method]]
    for (n in seq_along(element)) {
      contrast_name <- names(element)[n]
      contrast_data <- element[[contrast_name]]
      rekeyed_map[[contrast_name]][[method]] <- contrast_data
    }
  }

  numerators <- list()
  denominators <- list()
  mapped <- length(rekeyed_map)
  for (en in seq_along(rekeyed_map)) {
    entry <- rekeyed_map[[en]]
    entry_name <- names(rekeyed_map)[en]
    numerators[[entry_name]] <- list()
    denominators[[entry_name]] <- list()
    mesg("Starting on ", en, ": ", entry_name, ".")
    for (m in seq_along(entry)) {
      method <- names(entry)[m]
      internal <- entry[[method]]
      found_table <- internal[["string"]]
      wanted_numerator <- internal[["wanted_numerator"]]
      wanted_denominator <- internal[["wanted_denominator"]]
      numerators[[entry_name]][[method]] <- wanted_numerator
      denominators[[entry_name]][[method]] <- wanted_denominator
      if (isFALSE(internal[["string"]])) {
        warning("The table for ", entry_name, " using ", method,
                " does not appear in the pairwise data.")
      }
      if (length(internal[["idx"]]) > 1) {
        warning("There appear to be multiple tables for ", entry_name, " choosing the first.")
      }
    } ## End iterating over the methods
    combined <- combine_mapped_table(
      entry, includes, adjp = adjp, padj_type = padj_type,
      annot_df = annot_df, excludes = excludes,
      lfc_cutoff = lfc_cutoff, p_cutoff = p_cutoff, format_sig = format_sig,
      wanted_genes = wanted_genes, scale_p = scale_p)

    ## I am not sure if this should be the set of inverted tables yet.
    invert_colors <- FALSE
    invert_call <- sum(unlist(combined[["inverts"]]))
    if (invert_call > 0) {
      invert_colors <- TRUE
    }
    extracted[["data"]][[entry_name]] <- combined[["data"]]
    extracted[["table_names"]][[entry_name]] <- combined[["summary"]][["table"]]
    ## extracted[["kept"]] <- kept_tables
    extracted[["keepers"]] <- keepers
    plot_inputs <- list()
    for (i in seq_along(combined[["includes"]])) {
      plot_method <- names(combined[["includes"]])[i]
      plot_bool <- combined[["includes"]][[i]]
      if (isTRUE(plot_bool)) {
        plot_inputs[[plot_method]] <- datum[[plot_method]]
      }
    }

    ## Using the isTRUE() in case label_column is NULL, in which case
    ## NULL %in% stuff returns logical(0) and will fail an if().
    if (!isTRUE(label_column %in% colnames(combined[["data"]]))) {
      label_column <- NULL
    }

    ## Do not bother plotting if there are less than min_genes
    some_plots <- NULL
    if (nrow(combined[["data"]]) >= min_genes) {
      ## Changing this to a try() for when we have weirdo extra_contrasts.
      some_plots <- suppressWarnings(combine_extracted_plots(
        entry_name, combined, wanted_denominator, wanted_numerator, plot_inputs,
        loess = loess, lfc_cutoff = lfc_cutoff, p_cutoff = p_cutoff,
        adjp = adjp, found_table = found_table, p_type = padj_type,
        plot_colors = plot_colors, fancy = fancy,
        do_inverse = FALSE, invert_colors = invert_colors,
        z = z, alpha = alpha, z_lines = z_lines,
        label = label, label_column = label_column))
    }
    extracted[["plots"]][[entry_name]] <- some_plots
    extracted[["summaries"]] <- rbind(extracted[["summaries"]],
                                      as.data.frame(combined[["summary"]]))
  } ## Ending the for loop of elements in the keepers list.
  extracted[["numerators"]] <- numerators
  extracted[["denominators"]] <- denominators
  message("Looking for subscript invalid names, end of extract_keepers.")
  return(extracted)
}
setGeneric("extract_keepers")

#' @describeIn extract_keepers Use a character vector instead of a list.
#' @export
setMethod(
  "extract_keepers", signature = signature(extracted = "list", keepers = "character"),
  definition = function(extracted, keepers, table_names,
                        all_coefficients, apr,
                        adjp, annot_df, includes,
                        excludes, padj_type,
                        fancy = FALSE, loess = FALSE,
                        lfc_cutoff = 1.0, p_cutoff = 0.05,
                        format_sig = 4, plot_colors = plot_colors,
                        z = 1.5, alpha = 0.4, z_lines = FALSE,
                        label = 10, label_column = "hgnc_symbol",
                        scale_p = FALSE) {
    if (keepers[1] == "all") {
      new_keepers <- list()
      numerators <- denominators <- c()
      ## Note, I changed table_names to be sorted by method.  I can
      ## either iterate over every method, take the union of all, or
      ## arbitrarily choose a method...
      possible_names <- c()
      ## Limma and edger are the most likely to have extra contrasts,
      ## so check one of them first.
      if (!is.null(table_names[["limma"]])) {
        possible_names <- table_names[["limma"]]
      } else {
        possible_names <- table_names[[1]]
      }
      for (a in seq_along(possible_names)) {
        name <- possible_names[[a]]
        splitted <- strsplit(x = name, split = "_vs_")
        denominator <- splitted[[1]][2]
        numerator <- splitted[[1]][1]
        new_keepers[[name]] <- c(numerator, denominator)
      }
    } else {
      splitted <- strsplit(x = keepers, split = "_vs_")
      numerator <- splitted[[1]][1]
      denominator <- splitted[[1]][2]
      new_keepers <- list(splitted = c(numerator, denominator))
    }
    extract_keepers(extracted, new_keepers, table_names,
                    all_coefficients, apr,
                    adjp, annot_df, includes,
                    excludes, padj_type,
                    fancy = fancy, loess = loess,
                    lfc_cutoff = lfc_cutoff, p_cutoff = p_cutoff,
                    format_sig = format_sig, plot_colors = plot_colors,
                    z = z, alpha = alpha, z_lines = z_lines,
                    label = label, label_column = label_column,
                    scale_p = scale_p)
  })


#' Extract the sets of genes which are significantly more abundant than the rest.
#'
#' Given the output of something_pairwise(), pull out the genes for each contrast
#' which are the most/least abundant.  This is in contrast to extract_significant_genes().
#' That function seeks out the most changed, statistically significant genes.
#'
#' @param pairwise Output from _pairwise()().
#' @param according_to What tool(s) define 'most?'  One may use deseq, edger,
#'  limma, basic, all.
#' @param n How many genes to pull?
#' @param z Instead take the distribution of abundances and pull those past the
#'  given z score.
#' @param unique One might want the subset of unique genes in the top-n which
#'  are unique in the set of available conditions.  This will attempt to
#'  provide that.
#' @param excel Excel file to write.
#' @param ... Arguments passed into arglist.
#' @return The set of most/least abundant genes by contrast/tool.
#' @seealso \pkg{openxlsx}
#' @export
extract_abundant_genes <- function(pairwise, according_to = "deseq", n = 100,
                                   z = NULL, unique = FALSE,
                                   excel = "excel/abundant_genes.xlsx", ...) {
  arglist <- list(...)
  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]
  abundant_lists <- list()
  final_list <- list()

  data <- NULL
  if (according_to[[1]] == "all") {
    according_to <- c("limma", "deseq", "edger", "basic")
  }

  for (type in according_to) {
    datum <- pairwise[[type]]
    abundant_lists[[type]] <- get_abundant_genes(datum, type = type, n = n, z = z,
                                                 unique = unique)
    ## Something to note: the coefficient df from deseq (for example) includes coefficients
    ## across every factor in the model.  As a result, if one considers the default model of
    ## condition + batch, there will likely be coefficients which are not explicitly of interest.
  }

  if (class(excel)[1] == "character") {
    mesg("Writing a legend of columns.")
    legend <- data.frame(rbind(
      c("The first ~3-10 columns of each sheet:",
        "are annotations provided by our chosen annotation source for this experiment."),
      c("Next column", "The most/least abundant genes.")),
      stringsAsFactors = FALSE)
    colnames(legend) <- c("column name", "column definition")
    xls_result <- write_xlsx(wb, data = legend, sheet = "legend", rownames = FALSE,
                             title = "Columns used in the following tables.")
  }

  ## Now make the excel sheet for each method/coefficient
  ## Given that the set of coefficients may be more than we want, using
  ## names(abundant_lists[[according]][["high"]]) is probably not wise.
  ## Instead we should just be taking the levels of the condition factor.
  wanted_coefficients <- levels(pairwise[["deseq"]][["conditions"]])
  for (according in names(abundant_lists)) {
    for (coef in wanted_coefficients) {
      sheetname <- glue("{according}_high_{coef}")
      annotations <- rowData(pairwise[["input"]])
      high_abundances <- abundant_lists[[according]][["high"]][[coef]]
      kept_annotations <- names(high_abundances)
      kept_idx <- rownames(annotations) %in% kept_annotations
      kept_annotations <- annotations[kept_idx, ]
      high_data <- data.frame()
      if (nrow(annotations) > 0 && ncol(annotations) > 0) {
        high_df <- as.data.frame(high_abundances)
        rownames(high_df) <- names(high_abundances)
        high_data <- merge(high_df, annotations,
                           by = "row.names", all.x = TRUE)
        rownames(high_data) <- high_data[["Row.names"]]
        high_data[["Row.names"]] <- NULL
      } else {
        high_data <- as.data.frame(high_abundances)
      }
      start_row <- 1
      if (class(excel)[1] == "character") {
        title <- glue("Table SXXX: High abundance genes in {coef} according to {according}.")
        xls_result <- write_xlsx(data = high_data, wb = wb, sheet = sheetname, title = title)
        start_row <- start_row + xls_result[["end_row"]] + 2
      }

      sheetname <- glue("{according}_low_{coef}")
      low_abundances <- abundant_lists[[according]][["low"]][[coef]]
      kept_annotations <- names(low_abundances)
      kept_idx <- rownames(annotations) %in% kept_annotations
      kept_annotations <- annotations[kept_idx, ]
      low_data <- data.frame()
      if (nrow(annotations) > 0 && ncol(annotations) > 0) {
        low_df <- as.data.frame(low_abundances)
        rownames(low_df) <- names(low_abundances)
        low_data <- merge(data.frame(low_df), annotations,
                          by = "row.names", all.x = TRUE)
        rownames(low_data) <- low_data[["Row.names"]]
        low_data[["Row.names"]] <- NULL
      } else {
        low_data <- as.data.frame(low_abundances)
      }
      if (class(excel)[1] == "character") {
        title <- glue("Table SXXX: Low abundance genes in {coef} according to {according}.")
        xls_result <- write_xlsx(data = low_data, wb = wb, sheet = sheetname, title = title)
      }
    } ## End the for loop
  }

  if (class(excel)[1] == "character") {
    excel_ret <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  }
  ret <- list(
    "with_annotations" = final_list,
    "abundances" = abundant_lists)
  return(ret)
}

#' Alias for extract_significant_genes because I am dumb.
#'
#' @param ... The parameters for extract_significant_genes()
#' @return It should return a reminder for me to remember my function names or
#'   change them to something not stupid.
#' @export
extract_siggenes <- function(...) {
  extract_significant_genes(...)
}
#' Extract the sets of genes which are significantly up/down regulated
#' from the combined tables.
#'
#' Given the output from combine_de_tables(), extract the genes in
#' which we have the greatest likely interest, either because they
#' have the largest fold changes, lowest p-values, fall outside a
#' z-score, or are at the top/bottom of the ranked list.
#'
#' @param combined Output from combine_de_tables().
#' @param according_to What tool(s) decide 'significant?'  One may use
#'  the deseq, edger, limma, basic, meta, or all.
#' @param lfc Log fold change to define 'significant'.
#' @param p (Adjusted)p-value to define 'significant'.
#' @param sig_bar Add bar plots describing various cutoffs of 'significant'?
#' @param z Z-score to define 'significant'.
#' @param n Take the top/bottom-n genes.
#' @param min_mean_exprs Add a minimum expression value.
#' @param exprs_column Use this column to define expression.
#' @param top_percent Use a percentage to get the top-n genes.
#' @param p_type use an adjusted p-value?
#' @param invert_barplots Invert the significance barplots as per Najib's request?
#' @param excel Write the results to this excel file, or NULL.
#' @param fc_column Column in the DE data containing the foldchange values.
#' @param p_column Column in the DE data containing the pvalues.
#' @param column_suffix Used to help determine which columns are used
#'  to find significant genes via logfc/p-value.
#' @param gmt Write a gmt file using this result?
#' @param category When writing gmt files, set the category here.
#' @param fancy Write fancy plots with the xlsx file?
#' @param lfc_cutoffs 3 cutoffs for significant bar plots.
#' @param phenotype_name When writing gmt files, set the phenotype flag here.
#' @param set_name When writing gmt files, assign the set here.
#' @param current_id Choose the current ID type for an output gmt file.
#' @param comparison The cutoff may be '>|<' or '<=|>='.
#' @param required_id Choose the desired ID type for an output gmt file.
#' @param min_gmt_genes Define the minimum number of genes in a gene
#'  set for writing a gmt file.
#' @param ... Arguments passed into arglist.
#' @return The set of up-genes, down-genes, and numbers therein.
#' @seealso \code{\link{combine_de_tables}}
#' @export
extract_significant_genes <- function(combined, according_to = "all", lfc = 1.0,
                                      p = 0.05, sig_bar = TRUE, z = NULL, n = NULL,
                                      min_mean_exprs = NULL, exprs_column = NULL,
                                      top_percent = NULL, p_type = "adj",
                                      invert_barplots = FALSE, excel = NULL, fc_column = NULL,
                                      p_column = NULL, column_suffix = TRUE, gmt = FALSE,
                                      category = "category", fancy = FALSE, lfc_cutoffs = NULL,
                                      phenotype_name = "phenotype", set_name = "set",
                                      current_id = "ENSEMBL",
                                      comparison = "orequal", required_id = "ENTREZID",
                                      min_gmt_genes = 10, ...) {
  arglist <- list(...)
  image_files <- c()  ## For cleaning up tmp image files after saving the xlsx file.

  if (is.null(lfc_cutoffs)) {
    if (lfc == 0) {
      lfc_cutoffs <- c(0, 1, 2)
    } else {
      lfc_cutoffs <- c(0, lfc, lfc * 2)
    }
  }
  ## The following two if statements are a little silly, but will be used later
  ## To check for the existence of a column in the input data via p_column %in% colnames()
  if (is.null(fc_column)) {
    fc_column <- ""
  }
  if (is.null(p_column)) {
    p_column <- ""
  }

  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]
  num_tables <- 0
  table_names <- NULL
  all_tables <- NULL
  table_mappings <- NULL
  if (class(combined)[1] == "data.frame") {
    ## Then this is just a data frame.
    all_tables[["all"]] <- combined
    table_names <- "all"
    num_tables <- 1
    table_mappings <- table_names
  } else if (class(combined)[1] == "combined_de") {
    ## Then this is the result of combine_de_tables()
    num_tables <- length(names(combined[["data"]]))
    table_names <- combined[["table_names"]]
    all_tables <- combined[["data"]]
    table_mappings <- table_names
  } else if (!is.null(combined[["contrast_list"]])) {
    ## Then this is the result of all_pairwise()
    num_tables <- length(combined[["contrast_list"]])
    ## Extract the names of the tables which filled combined
    table_names <- names(combined[["data"]])
    ## Pull the table list
    all_tables <- combined[["data"]]
    ## Get the mappings of contrast_name -> table_name
    table_mappings <- combined[["keepers"]]
  } else {
    ## Then this is just a data frame.
    all_tables[["all"]] <- combined
    table_names <- "all"
    num_tables <- 1
    table_mappings <- table_names
  }

  if (!is.null(top_percent)) {
    n <- floor(nrow(all_tables[[1]]) * (top_percent / 100))
    mesg("Setting n to ", n)
  }
  if (!is.null(n)) {
    lfc <- NULL
    p <- NULL
  }

  logfc_suffix <- "_logfc"
  p_suffix <- "_p"
  adjp_suffix <- "_adjp"
  if (!isTRUE(column_suffix)) {
    logfc_suffix <- ""
    p_suffix <- ""
    adjp_suffix <- ""
  }

  trimmed_up <- list()
  trimmed_down <- list()
  up_titles <- list()
  down_titles <- list()
  sig_list <- list()
  title_append <- ""
  if (!is.null(lfc)) {
    title_append <- glue("{title_append} |log2fc| >= {lfc}")
  }
  if (!is.null(p)) {
    title_append <- glue("{title_append} p <= {p}")
  }
  if (!is.null(min_mean_exprs)) {
    title_append <- glue("{title_append} minimum expression >= {min_mean_exprs}")
  }
  if (!is.null(z)) {
    title_append <- glue("{title_append} |z| >= {z}")
  }
  if (!is.null(n)) {
    title_append <- glue("{title_append} top|bottom n={n}")
  }

  if (according_to[[1]] == "all") {
    according_to <- c("limma", "edger", "deseq", "ebseq", "basic")
  }

  if ("character" %in% class(excel)) {
    written <- write_sig_legend(wb)
    xls_result <- written[["xls_result"]]
  }

  ret <- list()
  sheet_count <- 0
  according_kept <- according_to
  chosen_columns <- list()
  ## Iterate over the according_to entries to see if there are valid p-value and logfc columns
  for (summary_count in seq_along(according_to)) {
    according <- according_to[summary_count]
    test_fc_param <- fc_column  ## Check if a column was provided by the user
    ## Otherwise make a column name from the method employed followed by the suffix.
    test_fc_column <- glue("{according}{logfc_suffix}")
    test_p_param <- p_column
    test_p_column <- glue("{according}{p_suffix}")
    test_adjp_column <- glue("{according}{adjp_suffix}")
    skip <- TRUE
    if (test_fc_param %in% colnames(combined[["data"]][[1]])) {
      chosen_columns[[according]][["fc"]] <- test_fc_param
      if (test_p_param %in% colnames(combined[["data"]][[1]])) {
        skip <- FALSE
        chosen_columns[[according]][["p"]] <- test_p_param
      }
      ## Adjusted p takes precedence.
      if (test_adjp_column %in% colnames(combined[["data"]][[1]])) {
        skip <- FALSE
        chosen_columns[[according]][["p"]] <- test_adjp_column
      }
    } else if (test_fc_column %in% colnames(combined[["data"]][[1]])) {
      skip <- FALSE
      chosen_columns[[according]][["fc"]] <- test_fc_column
      if (test_p_param %in% colnames(combined[["data"]][[1]])) {
        skip <- FALSE
        chosen_columns[[according]][["p"]] <- test_p_param
      }

      if (test_adjp_column %in% colnames(combined[["data"]][[1]])) {
        skip <- FALSE
        chosen_columns[[according]][["p"]] <- test_adjp_column
      }
    }
    ## We have given some chances to not skip this group, if skip never got set to FALSE
    ## then this is probably not a useful criterion for searching.
    if (isTRUE(skip)) {
      mesg("Did not find the ", test_fc_column, ", skipping ", according, ".")
      according_kept <- according_to[!according == according_to]
      next
    }
  } ## End checking set of according_tos for valid entries.
  ## Our list of chosen_columns should have two levels, the first by method
  ## The second with logfc and p, which contain the columns of interest.

  according <- NULL
  according_to <- according_kept
  for (summary_count in seq_along(according_to)) {
    according <- according_to[summary_count]
    ret[[according]] <- list()
    change_counts_up <- list()
    change_counts_down <- list()
    this_fc_column <- chosen_columns[[according]][["fc"]]
    this_p_column <- chosen_columns[[according]][["p"]]
    for (table_count in seq_along(table_names)) {
      table_name <- names(table_names)[table_count]
      plot_name <- as.character(table_names)[table_count]
      plot_name <- gsub(pattern = "-inverted", replacement = "", x = plot_name)

      this_table <- all_tables[[table_name]]
      ## Added this try() because I am not sure how I want to deal with
      ## extra contrasts, and until I decide it is easier to skip them
      trimming <- get_sig_genes(
        this_table, lfc = lfc, p = p, z = z, n = n, column = this_fc_column,
        min_mean_exprs = min_mean_exprs, exprs_column = exprs_column,
        comparison = comparison, p_column = this_p_column)
      if ("try-error" %in% class(trimming) || is.null(trimming)) {
        trimmed_up[[table_name]] <- data.frame()
        trimmed_down[[table_name]] <- data.frame()
        change_counts_up[[table_name]] <- 0
        change_counts_down[[table_name]] <- 0
        up_titles[[table_name]] <- "This table was not processed."
        down_titles[[table_name]] <- "This table was not processed."
      } else {
        trimmed_up[[table_name]] <- trimming[["up_genes"]]
        change_counts_up[[table_name]] <- nrow(trimmed_up[[table_name]])
        trimmed_down[[table_name]] <- trimming[["down_genes"]]
        change_counts_down[[table_name]] <- nrow(trimmed_down[[table_name]])
        up_title <- glue("Table SXXX: Genes deemed significantly up in \\
                       {table_name} with {title_append} according to {according}.")
        up_titles[[table_name]] <- up_title
        down_title <- glue("Table SXXX: Genes deemed significantly down in \\
                         {table_name} with {title_append} according to {according}.")
        down_titles[[table_name]] <- down_title
      }
    } ## End extracting significant genes for loop

    change_counts <- as.data.frame(cbind(as.numeric(change_counts_up),
                                         as.numeric(change_counts_down)))
    colnames(change_counts) <- c("up", "down")
    rownames(change_counts)[table_count] <- table_name

    summary_title <- glue("Counting the number of changed genes by contrast according to \\
                          {according} with {title_append}.")

    ret[[according]] <- list(
      "ups" = trimmed_up,
      "downs" = trimmed_down,
      "counts" = change_counts,
      "up_titles" = up_titles,
      "down_titles" = down_titles,
      "counts_title" = summary_title)

    ## I want to start writing out msigdb compatible gmt files and therefore
    ## want to start creating gene set collections from our data.

    do_excel <- TRUE
    if (is.null(excel)) {
      do_excel <- FALSE
    }
    if (isFALSE(excel)) {
      do_excel <- FALSE
    }

    if (isTRUE(do_excel)) {
      mesg("Printing significant genes to the file: ", excel)
      xlsx_ret <- print_ups_downs(ret[[according]], wb, excel_basename, according = according,
                                  summary_count = summary_count)
      image_files <- c(image_files, xlsx_ret[["image_files"]])
      ## This is in case writing the sheet resulted in it being shortened.
      ## wb <- xlsx_ret[["workbook"]]
    } ## End of an if whether to print the data to excel
  } ## End list of according_to's

  ## the extraneous message() statements and instead fill that information into
  ## this data frame.
  name_element <- according_to[1]
  summary_df <- data.frame(row.names = names(ret[[name_element]][["ups"]]))

  sig_bar_plots <- NULL
  if (isTRUE(do_excel) && isTRUE(sig_bar)) {
    ## This needs to be changed to get_sig_genes()
    sig_bar_plots <- significant_barplots(
      combined, lfc_cutoffs = lfc_cutoffs, invert = invert_barplots,
      p = p, z = z, p_type = p_type,
      according_to = according_to,
      ...)
    plot_row <- 1
    plot_col <- 1
    mesg("Adding significance bar plots.")

    num_tables <- length(according_to)
    plot_row <- plot_row + ((nrow(change_counts) + 1) * num_tables) + 4
    ## The +4 is for the number of tools.
    ## I know it is silly to set the row in this very explicit fashion, but I
    ## want to make clear the fact that the table has a title, a set of
    ## headings, a length corresponding to the number of contrasts,  and then
    ## the new stuff should be added. Now add in a table summarizing the numbers
    ## in the plot. The information required to make this table is in
    ## sig_bar_plots[["ups"]][["limma"]] and sig_bar_plots[["downs"]][["limma"]]
    plot_row <- plot_row + 3
    ## I messed up something here.  The plots and tables
    ## at this point should start:
    ## 5(blank spaces and titles) + 4(table headings) + 4 * the number of contrasts.
    checked <- check_xlsx_worksheet(wb, "number_changed")

    for (according in according_to) {
      tmp_df <- ret[[according]][["counts"]]
      rownames(tmp_df) <- names(ret[[according]][["ups"]])
      colnames(tmp_df) <- paste0(according, "_", colnames(tmp_df))
      summary_df <- cbind(summary_df, tmp_df)
      sig_message <- as.character(glue("Significant {according} genes."))
      xls_result <- openxlsx::writeData(
        wb = wb, sheet = "number_changed", x = sig_message,
        startRow = plot_row, startCol = plot_col)
      plot_row <- plot_row + 1
      plotname <- glue("sigbar_{according}")
      try_result <- xlsx_insert_png(
        a_plot = sig_bar_plots[[according]], wb = wb, sheet = "number_changed",
        plotname = plotname, savedir = excel_basename, width = 9, height = 6,
        start_row = plot_row, start_col = plot_col, fancy = fancy)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
      summary_row <- plot_row
      summary_col <- plot_col + 11
      de_summary <- summarize_ups_downs(sig_bar_plots[["ups"]][[according]],
                                        sig_bar_plots[["downs"]][[according]])
      xls_summary <- write_xlsx(
        data = de_summary, wb = wb, sheet = "number_changed", rownames = TRUE,
        start_row = summary_row, start_col = summary_col)
      plot_row <- plot_row + 30
    } ## End for loop writing out significance bar plots
  } ## End if we want significance bar plots
  ret[["sig_bar_plots"]] <- sig_bar_plots
  summary_df[["rownames"]] <- NULL
  ret[["summary_df"]] <- summary_df
  ret[["lfc"]] <- lfc
  ret[["p"]] <- p
  ret[["z"]] <- z
  ret[["n"]] <- n
  ret[["according"]] <- according_to
  ret[["p_type"]] <- p_type

  if (isTRUE(do_excel)) {
    excel_ret <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  }

  for (img in image_files) {
    removed <- try(suppressWarnings(file.remove(img)), silent = TRUE)
  }
  class(ret) <- c("sig_genes", "list")

  ## If there is no gmt fielname provided, but an excel filename is provided,
  ## save a .gmt with the xlsx file's basename
  if (is.null(gmt)) {
    if (!is.null(excel_basename)) {
      output_dir <- xlsx[["dirname"]]
      output_base <- paste0(xlsx[["basename"]], ".gmt")
      gmt <- file.path(output_dir, output_base)
    }
  }

  if (!isFALSE(gmt)) {
    mesg("Going to attempt to create gmt files from these results.")
    annotation_name <- annotation(combined[["input"]][["input"]])
    gsc <- try(make_gsc_from_significant(ret, according_to = according_to, orgdb = annotation_name,
                                         pair_names = c("ups", "downs"), category_name = category,
                                         phenotype_name = phenotype_name, set_name = set_name,
                                         current_id = current_id, required_id = required_id,
                                         min_gmt_genes = min_gmt_genes))
    if (! "try-error" %in% class(gsc)) {
      types <- c("both", "up", "down")
      for (t in types) {
        if (!is.null(gsc[[t]])) {
          for (g in seq_along(gsc[[t]])) {
            datum <- gsc[[t]][[g]]
            contrast_name <- names(gsc[[t]])[g]
            dname <- dirname(gmt)
            bname <- gsub(x = basename(gmt), pattern = "\\.gmt", replacement = "")
            newname <- paste0(bname, "_", t, "_", contrast_name, ".gmt")
            write_to <- file.path(dname, newname)
            written <- GSEABase::toGmt(datum, write_to)
          }
        }
      }
    }
  }
  return(ret)
}

#' Print some significantly differentially expressed genes.
#'
#' @param x List containing the parameters used, gene subset tables,
#'  plots, xlsx output file, etc.
#' @param ... Other args to match the generic.
#' @export
print.sig_genes <- function(x, ...) {
  message("A set of genes deemed significant according to ", toString(x[["according"]]), ".")
  params <- ""
  if (!is.null(x[["lfc"]])) {
    params <- glue("LFC cutoff: {x[['lfc']]}")
  }
  if (!is.null(x[["p"]])) {
    params <- glue("{params} {x[['p_type']]} P cutoff: {x[['p']]}")
  }
  if (!is.null(x[["z"]])) {
    params <- glue("{params} Z cutoff: {x[['z']]}")
  }
  if (!is.null(x[["n"]])) {
    params <- glue("{params} topN: {x[['n']]}")
  }
  message("The parameters defining significant were:")
  message(params)
  print(x[["summary_df"]])
  if (sum(rowSums(x[["summary_df"]])) > 0) {
    plot(x[["sig_bar_plots"]][["deseq"]])
  }
  return(invisible(x))
}

#' Create a summary table of the ranges of fold-change values of potential interest.
#'
#' The columns have names with explicit lfc values, but the numbers which get put in them
#' may represent any arbitrary cutoff employed by the caller.
#'
#' @param ups The set of ups!
#' @param downs and downs!
summarize_ups_downs <- function(ups, downs) {
  ## The ups and downs tables have 1 row for each contrast, 3 columns of numbers named
  ## 'a_up_inner', 'b_up_middle', 'c_up_outer'.
  ups <- ups[, -1]
  downs <- downs[, -1]
  ups[[1]] <- as.numeric(ups[[1]])
  ups[[2]] <- as.numeric(ups[[2]])
  ups[[3]] <- as.numeric(ups[[3]])
  ups[["up_sum"]] <- rowSums(ups)
  downs[[1]] <- as.numeric(downs[[1]])
  downs[[2]] <- as.numeric(downs[[2]])
  downs[[3]] <- as.numeric(downs[[3]])
  downs[["down_sum"]] <- rowSums(downs)
  summary_table <- as.data.frame(cbind(ups, downs))
  summary_table <- summary_table[, c(1, 2, 3, 5, 6, 7, 4, 8)]
  colnames(summary_table) <- c("up_from_0_to_2", "up_from_2_to_4", "up_gt_4",
                               "down_from_0_to_2", "down_from_2_to_4", "down_gt_4",
                               "sum_up", "sum_down")
  summary_table[["up_gt_2"]] <- summary_table[["up_from_2_to_4"]] +
    summary_table[["up_gt_4"]]
  summary_table[["down_gt_2"]] <- summary_table[["down_from_2_to_4"]] +
    summary_table[["down_gt_4"]]
  summary_table_idx <- rev(rownames(summary_table))
  summary_table <- summary_table[summary_table_idx, ]
  return(summary_table)
}

#' Find the sets of intersecting significant genes
#'
#' Use extract_significant_genes() to find the points of agreement between
#' limma/deseq/edger.
#'
#' @param combined Result from combine_de_tables().
#' @param lfc Define significant via fold-change.
#' @param p Or p-value.
#' @param padding_rows How much space to put between groups of data?
#' @param z Use a z-score filter?
#' @param p_type Use normal or adjusted p-values.
#' @param selectors List of methods to intersect.
#' @param order When set to the default 'inverse', go from the set with the most
#' least intersection to the most. E.g. Start with abc,bc,ac,c,ab,b,a as
#'  opposed to a,b,ab,c,ac,bc,abc.
#' @param excel An optional excel workbook to which to write.
#' @param ... Extra arguments for extract_significant_genes() and friends.
#' @return List containing the intersections between the various DE methods for
#'  both the up and down sets of genes.  It should also provide some venn
#'  diagrams showing the degree of similarity between the methods.
#' @examples
#'  \dontrun{
#'   expt <- create_expt(metadata="some_metadata.xlsx", gene_info=funkytown)
#'   big_result <- all_pairwise(expt, model_batch=FALSE)
#'   pretty <- combine_de_tables(big_result, excel="excel/combined_expt.xlsx")
#'   intersect <- intersect_significant(pretty, excel="excel/intersecting_genes.xlsx")
#'  }
#' @export
intersect_significant <- function(combined, lfc = 1.0, p = 0.05, padding_rows = 2,
                                  z = NULL, p_type = "adj", selectors = c("limma", "deseq", "edger"),
                                  order = "inverse", excel = "excel/intersect_significant.xlsx",
                                  ...) {
  ## Check the set of first->sixth and see that they exist in the combined table.
  arglist <- list(...)
  image_files <- c()
  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]
  chosen_selectors <- c()
  extract_selectors <- c()
  alternate_selectors <- c()
  ## Use the column names from the first table to make this decision.
  possible_columns <- colnames(combined[["data"]][[1]])
  for (c in seq_along(selectors)) {
    rawname <- selectors[c]
    conjugate <- glue("{rawname}_logfc")
    if (rawname %in% possible_columns) {
      chosen_selectors <- c(chosen_selectors, rawname)
      alternate_selectors <- c(alternate_selectors, rawname)
    } else if (conjugate %in% possible_columns) {
      chosen_selectors <- c(chosen_selectors, rawname)
      extract_selectors <- c(extract_selectors, rawname)
    } else if (arglist[["fc_column"]] %in% possible_columns &&
                 arglist[["p_column"]] %in% possible_columns) {
      chosen_selectors <- c(chosen_selectors, "alternate")
      alternate_selectors <- c(alternate_selectors, "alternate")
    } else {
      mesg("Skipping ", rawname)
    }
  }
  if (length(chosen_selectors) < 2) {
    stop("This requires more selectors in order to create an intersection.")
  }

  sig_genes <- sm(extract_significant_genes(combined, lfc = lfc, p = p,
                                            z = z, p_type = p_type, excel = NULL))
  if (length(alternate_selectors) > 0) {
    alt_genes <- extract_significant_genes(combined, according_to = "alternate", lfc = lfc,
                                           fc_column = arglist[["fc_column"]],
                                           p_column = arglist[["p_column"]], ma = FALSE,
                                           p = p, z = z, p_type = p_type, excel = FALSE)
    sig_genes[["alternate"]] <- alt_genes[["alternate"]]
  }

  xls_result <- NULL
  ## Set up the base data structure, a list of ups and a list of downs.
  lst <- list("ups" = list(), "downs" = list())
  set_names <- c()
  venn_row <- 1
  venn_col <- 1
  ## Arbitrarily pull the set of ups from the first 'chooser'
  table_lst <- names(sig_genes[[1]][["ups"]])
  for (t in seq_along(table_lst)) {
    table <- table_lst[t]
    lst[["ups"]][[table]] <- list()
    lst[["downs"]][[table]] <- list()
    for (dir in c("ups", "downs")) {
      for (i in seq_along(chosen_selectors)) {
        name <- chosen_selectors[i]
        lst[[dir]][[table]][[name]] <- rownames(sig_genes[[name]][[dir]][[table]])
      } ## End pulling the significants by selectors.
      sets <- Vennerable::Venn(Sets = lst[[dir]][[table]])
      intersections <- sets@IntersectionSets
      tmp_file <- tmpmd5file(pattern = "venn", fileext = ".png")
      this_plot <- png(filename = tmp_file)
      controlled <- dev.control("enable")
      plt <- Vennerable::plot(sets, doWeights = FALSE)
      rec <- grDevices::recordPlot()
      plotted <- dev.off()
      removed <- file.remove(tmp_file)
      lst[[dir]][[table]][["sets"]] <- sets
      lst[[dir]][[table]][["intersections"]] <- intersections
      lst[[dir]][[table]][["plot"]] <- rec
      print_order <- names(lst[[dir]][[table]][["intersections"]])
      if (order == "inverse") {
        print_order <- rev(print_order)
      }
      set_names <- list()
      inner_count <- 0
      for (symbolic in print_order) {
        inner_count <- inner_count + 1
        symbols <- strsplit(as.character(symbolic), "")[[1]]
        name <- c()
        for (s in seq_along(symbols)) {
          symbol <- symbols[s]
          if (symbol == 1) {
            name <- c(name, chosen_selectors[s])
          }
        }
        name <- toString(name)
        set_names[[inner_count]] <- name
      } ## End for symbolic in names(elements_per_set)
      set_names[length(set_names)] <- "none"
      set_names[1] <- "all"

      ## This has been the source of some confusion for me.
      ## I think that the
      if (order == "inverse") {
        names(set_names) <- rev(names(lst[[dir]][[table]][["intersections"]]))
      } else {
        set_names <- rev(set_names)
        names(invert_names) <- names(lst[[dir]][[table]][["intersections"]])
      }

      lst[[dir]][[table]][["set_names"]] <- set_names
      table_rows <- combined[["data"]][[table]]
      xlsx_row <- 1
      xlsx_table <- ""
      lst[[dir]][[table]][["data"]] <- list()
      for (s in seq_along(set_names)) {
        sname <- set_names[s]
        set_selection_name <- names(sname)
        clean_sname <- sname
        clean_sname <- gsub(pattern = "\\, ", replacement = "_", x = sname)

        lst[[dir]][[table]][["data"]][[clean_sname]] <- data.frame()
        ## This next operation was previously being done on the list index
        ## number, which is exactly wrong.
        if (length(lst[[dir]][[table]][["intersections"]][[set_selection_name]]) > 0) {
          idx <- lst[[dir]][[table]][["intersections"]][[set_selection_name]]
          table_subset <- table_rows[idx, ]
          lst[[dir]][[table]][["data"]][[clean_sname]] <- table_subset
          if (dir == "ups") {
            text_dir <- "up"
          } else {
            text_dir <- "down"
          }
          xlsx_table <- glue("{text_dir}_{table}")
          xlsx_title <- glue("Genes deemed {text_dir} significant via logFC: {lfc}\\
                             , p-value: {p}; by {clean_sname}.")
          if (!is.null(excel)) {
            xl_result <- write_xlsx(data = table_subset, wb = wb, sheet = xlsx_table,
                                    start_row = xlsx_row, title = xlsx_title)
            xlsx_row <- xlsx_row + nrow(table_subset) + padding_rows + 2
          } ## End checking to write the excel file.
        }
      } ## End iterating over writing the set names.
    } ## End ups vs. downs

    up_plot <- lst[["ups"]][[table]][["plot"]]
    down_plot <- lst[["downs"]][[table]][["plot"]]
    venn_title <- glue("Summary of intersections among {toString(chosen_selectors)} \\
                       for {table}.")

    summary_df <- rbind(t(Vennerable::Weights(lst[["ups"]][[table]][["sets"]])),
                        t(Vennerable::Weights(lst[["downs"]][[table]][["sets"]])))
    rownames(summary_df) <- c("up", "down")
    tmp_colnames <- rev(lst[["ups"]][[table]][["set_names"]])
    tmp_colnames[length(tmp_colnames)] <- "all"
    colnames(summary_df) <- tmp_colnames
    summary_df <- summary_df[, -1]
    lst[["summary"]] <- summary_df
    if (!is.null(excel)) {
      xl_result <- write_xlsx(wb = wb,
                              data = summary_df,
                              sheet = "summary", title = venn_title,
                              start_row = venn_row, start_col = venn_col)
      venn_row <- venn_row + 4
      try_result <- xlsx_insert_png(
        up_plot, wb = wb, sheet = "summary", width = 6, height = 6,
        start_col = venn_col, start_row = venn_row)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
      venn_col <- venn_col + 12
      try_result <- xlsx_insert_png(
        down_plot, wb = wb, sheet = "summary", width = 6, height = 6,
        start_col = venn_col, start_row = venn_row)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
      venn_row <- venn_row + 26
    }
  }  ## End iterating over the tables

  if (!is.null(excel)) {
    excel_ret <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  }
  for (img in image_files) {
    removed <- try(suppressWarnings(file.remove(img)), silent = TRUE)
  }

  class(lst) <- c("sig_intersect", "list")
  return(lst)
}

#' Print the intersection of significant genes from multiple analyses.
#'
#' @param x List containing some venn diagrams, summaries of
#'  intersections, subsets of the intersections, etc.
#' @param ... Other args to match the generic.
#' @export
print.sig_intersect <- function(x, ...) {
  message("A set of genes deemed significant by multiple methods.")
  return(invisible(x))
}


#' Reprint the output from extract_significant_genes().
#'
#' I found myself needing to reprint these excel sheets because I
#' added some new information. This shortcuts that process for me.
#'
#' @param upsdowns Output from extract_significant_genes().
#' @param wb Workbook object to use for writing, or start a new one.
#' @param excel_basename Used when including plots in the xlsx sheet.
#' @param according Use limma, deseq, or edger for defining 'significant'.
#' @param summary_count For spacing sequential tables one after another.
#' @param ma Include ma plots?
#' @param fancy Print fancy plots with the xlsx file?x
#' @return Return from write_xlsx.
#' @seealso \code{\link{combine_de_tables}}
#' @export
print_ups_downs <- function(upsdowns, wb, excel_basename, according = "limma",
                            summary_count = 1, ma = FALSE, fancy = FALSE) {
  local_image_files <- c()
  ups <- upsdowns[["ups"]]
  downs <- upsdowns[["downs"]]
  up_titles <- upsdowns[["up_titles"]]
  down_titles <- upsdowns[["down_titles"]]
  summary <- as.data.frame(upsdowns[["counts"]])
  summary_title <- upsdowns[["counts_title"]]
  ma_plots <- upsdowns[["ma_plots"]]
  table_count <- 0
  summary_count <- summary_count - 1
  num_tables <- length(names(ups))
  summary_start <- ((num_tables + 2) * summary_count) + 1
  xls_summary_result <- write_xlsx(wb = wb, data = summary, start_col = 1,
                                   start_row = summary_start,
                                   sheet = "number_changed", title = summary_title)
  xls_result <- NULL
  for (table_count in seq_along(names(ups))) {
    base_name <- names(ups)[table_count]
    up_name <- glue("up_{according}_{base_name}")
    down_name <- glue("down_{according}_{base_name}")
    up_table <- ups[[table_count]]
    up_title <- up_titles[[table_count]]
    if (nrow(up_table) > 0) {
      mesg(table_count, "/", num_tables, ": Creating significant table ", up_name)
      xls_result <- write_xlsx(data = up_table, wb = wb, sheet = up_name, title = up_title)
      ## This is in case the sheet name is past the 30 character limit.
      sheet_name <- xls_result[["sheet"]]
      if (isTRUE(ma)) {
        ma_row <- 1
        ma_col <- xls_result[["end_col"]] + 1
        if (!is.null(ma_plots[[base_name]])) {
          plot_name <- glue("ma_{according}_{base_name}")
          try_result <- xlsx_insert_png(ma_plots[[base_name]], wb = wb, sheet = sheet_name,
                                        plotname = plot_name, savedir = excel_basename,
                                        start_row = ma_row, start_col = ma_col, fancy = fancy)
          if (! "try-error" %in% class(try_result)) {
            local_image_files <- c(local_image_files, try_result[["filename"]])
          }
        }
      }
    } else {
      mesg("The up table ", base_name, " is empty.")
    }

    down_table <- downs[[table_count]]
    down_title <- down_titles[[table_count]]
    if (nrow(down_table) > 0) {
      xls_result <- write_xlsx(data = down_table, wb = wb, sheet = down_name, title = down_title)
    } else {
      mesg("The down table ", base_name, " is empty.")
    }
  } ## End for each name in ups
  retlist <- list(
    "result" = xls_result,
    "image_files" = local_image_files)
  return(retlist)
}

#' Write the legend of an excel file for combine_de_tables()
#'
#' @param wb Workbook to write
#' @param excel_basename Where to write it
#' @param plot_dim Default plot size.
#' @param apr The all_pairwise() result.
#' @param includes List of methods to include in the legend.
#' @param includes List of booleans defining which methods to examine.
#' @param padj_type P-adjustment employed.
#' @param fancy Write fancy plots with the xlsx file?
write_combined_legend <- function(wb, excel_basename, plot_dim, apr,
                                  includes, padj_type, fancy = FALSE) {
  ## I want to print a string reminding the user what kind of model was used in
  ## the analysis. Do that here.  Noting that if 'batch' is actually from a
  ## surrogate variable, then we will not have TRUE/FALSE but instead a matrix.
  do_excel <- TRUE
  if (length(excel_basename) == 0) {
    do_excel <- FALSE
  }
  if (is.null(wb)) {
    do_excel <- FALSE
  }
  reminder_model_cond <- apr[["model_cond"]]
  reminder_model_batch <- apr[["model_batch"]]
  reminder_extra <- apr[["extra_contrasts"]]
  reminder_string <- NULL
  if (class(reminder_model_batch)[1] == "matrix") {
    reminder_string <- "The contrasts were performed using surrogates from sva/ruv/etc."
  } else if (isTRUE(reminder_model_batch) && isTRUE(reminder_model_cond)) {
    reminder_string <- "The contrasts were performed with condition and batch in the model."
  } else if (isTRUE(reminder_model_cond)) {
    reminder_string <- "The contrasts were performed with only condition in the model."
  } else {
    reminder_string <- "The contrasts were performed in a strange way, beware!"
  }

  ## The next large set of data.frame() calls create the first sheet, containing a legend.
  mesg("Writing a legend of columns.")
  legend <- data.frame(rbind(
    c("The first ~3-10 columns of each sheet:",
      "are annotations provided by our chosen annotation source for this experiment."),
    c("Next 6 columns", "The logFC and p-values reported by limma, edger, and deseq2.")
  ), stringsAsFactors = FALSE)
  basic_legend <- data.frame(rbind(
    c("The next 8 columns", "Statistics generated by the basic analysis."),
    c("basic_nummed", "log2 median values of the numerator for this comparison."),
    c("basic_denmed", "log2 median values of the denominator for this comparison."),
    c("basic_numvar", "Variance observed in the numerator values."),
    c("basic_denvar", "Variance observed in the denominator values."),
    c("basic_logfc", "The log2 fold change observed by the basic analysis."),
    c("basic_t", "T-statistic from basic."),
    c("basic_p", "Resulting p-value."),
    c(glue("basic_adjp_{padj_type}"), glue("p-value adjusted with {padj_type}")),
    c("basic_adjp", "BH correction of the p-value.")
  ), stringsAsFactors = FALSE)
  deseq_legend <- data.frame(rbind(
    c("The next 7 columns", "Statistics generated by DESeq2."),
    c("deseq_logfc_rep", "The log2 fold change reported by DESeq2, again."),
    c("deseq_adjp_rep", "The adjusted-p value reported by DESeq2, again."),
    c("deseq_basemean", "The base mean of all samples according to DESeq2."),
    c("deseq_lfcse", "The standard error observed given the log2 fold change."),
    c("deseq_stat", "T-statistic reported by DESeq2 given the log2FC and observed variances."),
    c("deseq_p", "Resulting p-value."),
    c(glue("deseq_adjp_{padj_type}"), glue("p-value adjusted with {padj_type}")),
    c("deseq_q", "False-positive corrected p-value."),
    c("deseq_num", "Numerator coefficients from deseq for this contrast."),
    c("deseq_den", "Denominator coefficients from deseq for this contrast.")
  ), stringsAsFactors = FALSE)
  dream_legend <- data.frame(rbind(
    c("The next columns", "Statistics generated by variancepartition/dream."),
    c("dream_logfc", "The log2 fold change reported by dream."),
    c("deseq_adjp", "The adjusted-p value reported by dream.")
  ), stringsAsFactors = FALSE)
  ebseq_legend <- data.frame(rbind(
    c("The next 6 columns", "Statistics generated by EBSeq."),
    c("ebseq_FC", "The fold change reported by EBSeq."),
    c("ebseq_logfc", "The log2 fold change from EBSeq."),
    c("ebseq_postfc", "The post-probability fold change."),
    c("ebseq_mean", "Mean of the EBSeq values."),
    c("PPEE", "Post-probability that the numerator/denominator are equivalent."),
    c("PPDE",  "Post-probability that the numerator/denominator are different."),
    c("ebseq_adjp",  "Attempt at FDR correction of the PPEE, this is just copied from PPEE
and is in _no_ way statistically valid, but added as a plotting conveinence.")
  ), stringsAsFactors = FALSE)
  edger_legend <- data.frame(rbind(
    c("The next 7 columns", "Statistics generated by edgeR."),
    c("edger_logfc_rep", "The log2 fold change reported by edgeR, again."),
    c("edger_adjp_rep", "The adjusted-p value reported by edgeR, again."),
    c("edger_logcpm",
      "Similar DESeq2's basemean, except only including the samples in the comparison."),
    c("edger_lr", "Undocumented, I think it is the T-statistic calculated by edgeR."),
    c("edger_p", "The observed p-value from edgeR."),
    c(glue("edger_adjp_{padj_type}"), glue("p-value adjusted with {padj_type}")),
    c("edger_q", "The observed corrected p-value from edgeR.")
  ), stringsAsFactors = FALSE)
  limma_legend <- data.frame(rbind(
    c("The next 7 columns", "Statistics generated by limma."),
    c("limma_logfc_rep", "The log2 fold change reported by limma, again."),
    c("limma_adjp_rep", "The adjusted-p value reported by limma, again."),
    c("limma_ave", "Average log2 expression observed by limma across all samples."),
    c("limma_t", "T-statistic reported by limma given the log2FC and variances."),
    c("limma_p", "Derived from limma_t, the p-value asking 'is this logfc significant?'"),
    c(glue("limma_adjp_{padj_type}"), glue("p-value adjusted with {padj_type}")),
    c("limma_b", "Bayesian estimate of the log-odds significance."),
    c("limma_q", "A q-value FDR adjustment of the p-value above.")
  ), stringsAsFactors = FALSE)
  noiseq_legend <- data.frame(rbind(
    c("The next 7 columns", "Statistics generated by the basic analysis."),
    c("noiseq_numean", "log2 median values of the numerator for this comparison."),
    c("noiseq_denmean", "log2 median values of the denominator for this comparison."),
    c("noiseq_theta", "The theta value from noiseq, I forgot what this means (hey you reread that paper)."),
    c("noiseq_prob", "I assume this is the t statistic produced by noiseq, again I need to reread the paper."),
    c("noiseq_logfc", "The log2 fold change observed by noiseq analysis."),
    c("noiseq_p", "The p-value from noiseq."),
    c("noiseq_adjp", "The adjusted p-value, again I haven't spent enough time to know the details.")
  ), stringsAsFactors = FALSE)
  summary_legend <- data.frame(rbind(
    c("The next 5 columns", "Summaries of the limma/deseq/edger results."),
    c("lfc_meta", "The mean fold-change value of limma/deseq/edger."),
    c("lfc_var", "The variance between limma/deseq/edger."),
    c("lfc_varbymed", "The ratio of the variance/median (lower means better agreement.)"),
    c("p_meta", "A meta-p-value of the mean p-values."),
    c("p_var", "Variance among the 3 p-values."),
    c("The last columns: top plot left",
      "Venn of genes with logFC > 0 and p-value <= 0.05 for limma/DESeq/Edger."),
    c("The last columns: top plot right",
      "Venn of genes with logFC < 0 and p-value <= 0.05 for limma/DESeq/Edger."),
    c("The last columns: second plot",
      "Scatter plot of the voom-adjusted/normalized counts for each coefficient."),
    c("The last columns: third plot",
      "Scatter plot of the adjusted/normalized counts for each coefficient from edgeR."),
    c("The last columns: fourth plot",
      "Scatter plot of the adjusted/normalized counts for each coefficient from DESeq."),
    c("", "If this data was adjusted with sva, check for a sheet 'original_pvalues' at the end.")
  ), stringsAsFactors = FALSE)

  ## Here we including only those columns which are relevant to the analysis performed.
  if (isTRUE(includes[["limma"]])) {
    legend <- rbind(legend,
                    c("limma_logfc", "The log2 fold change reported by limma."),
                    c("limma_adjp", "The adjusted-p value reported by limma."))
  }
  if (isTRUE(includes[["deseq"]])) {
    legend <- rbind(legend,
                    c("deseq_logfc", "The log2 fold change reported by DESeq2."),
                    c("deseq_adjp", "The adjusted-p value reported by DESeq2."))
  }
  if (isTRUE(includes[["edger"]])) {
    legend <- rbind(legend,
                    c("edger_logfc", "The log2 fold change reported by edgeR."),
                    c("edger_adjp", "The adjusted-p value reported by edgeR."))
  }
  if (isTRUE(includes[["edger"]])) {
    legend <- rbind(legend,
                    c("edger_logfc", "The log2 fold change reported by edgeR."),
                    c("edger_adjp", "The adjusted-p value reported by edgeR."))
  }

  ## Simplify the code a little by bringing the table definitions up here.
  ## Figure out the set of possible contrasts by name.
  ## I am changing this logic to have it include the union of table names
  ## rather than just pick up the first set that worked.
  table_names <- list()
  if (isTRUE(includes[["basic"]])) {
    legend <- rbind(legend, basic_legend)
    table_names[["basic"]] <- apr[["basic"]][["contrasts_performed"]]
  } else {
    basic <- NULL
  }
  if (isTRUE(includes[["deseq"]])) {
    legend <- rbind(legend, deseq_legend)
    table_names[["deseq"]] <- apr[["deseq"]][["contrasts_performed"]]
  } else {
    deseq <- NULL
  }
  if (isTRUE(includes[["dream"]])) {
    legend <- rbind(legend, dream_legend)
    table_names[["dream"]] <- apr[["dream"]][["contrasts_performed"]]
  } else {
    dream <- NULL
  }
  if (isTRUE(includes[["ebseq"]])) {
    legend <- rbind(legend, ebseq_legend)
    table_names[["ebseq"]] <- names(apr[["ebseq"]][["all_tables"]])
  } else {
    ebseq <- NULL
  }
  if (isTRUE(includes[["edger"]])) {
    legend <- rbind(legend, edger_legend)
    table_names[["edger"]] <- apr[["edger"]][["contrasts_performed"]]
  } else {
    edger <- NULL
  }
  if (isTRUE(includes[["limma"]])) {
    legend <- rbind(legend, limma_legend)
    table_names[["limma"]] <- apr[["limma"]][["contrasts_performed"]]
  } else {
    limma <- NULL
  }
  if (isTRUE(includes[["noiseq"]])) {
    legend <- rbind(legend, noiseq_legend)
    table_names[["noiseq"]] <- apr[["noiseq"]][["contrasts_performed"]]
  } else {
    noiseq <- NULL
  }

  ## Make sure there were no errors and die if things went catastrophically wrong.
  if (sum(unlist(includes)) < 1) {
    stop("None of the DE tools appear to have worked.")
  }
  if (length(table_names) == 0) {
    stop("Could not find the set of table names.")
  }

  if (isTRUE(includes[["limma"]]) && isTRUE(includes[["deseq"]]) &&
        isTRUE(includes[["edger"]]) && isTRUE(includes[["basic"]])) {
    legend <- rbind(legend, summary_legend)
  }
  colnames(legend) <- c("column name", "column definition")
  xls_result <- NULL
  full_title <- paste0("Columns used in the following tables.  ", reminder_string)

  if (isTRUE(do_excel)) {
    xls_result <- write_xlsx(
      wb, data = legend, sheet = "legend", rownames = FALSE,
      title = full_title)
  }

  ## Some folks have asked for some PCA showing the before/after surrogates.
  ## Put that on the first sheet, then.
  ## This if (isTRUE()) is a little odd, perhaps it should be removed or moved up.
  image_files <- c()
  if (isTRUE(do_excel)) {
    mesg("Printing pca plots before and after surrogate|batch estimation.")
    ## Add PCA before/after
    chosen_estimate <- apr[["batch_type"]]
    xl_result <- openxlsx::writeData(
      wb = wb, sheet = "legend",
      x = "PCA plot before surrogate estimation.",
      startRow = 1, startCol = 10)
    try_result <- xlsx_insert_png(
      apr[["pre_batch"]][["plot"]], wb = wb, sheet = "legend", start_row = 2,
      width = (plot_dim * 3/2), height = plot_dim, start_col = 10,
      plotname = "pre_pca", savedir = excel_basename,
      fancy = fancy)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }
    xl_result <- openxlsx::writeData(
      wb = wb, sheet = "legend",
      x = as.character(glue("PCA after surrogate estimation with: {chosen_estimate}")),
      startRow = 36, startCol = 10)
    try_result <- xlsx_insert_png(
      apr[["post_batch"]][["plot"]], wb = wb, sheet = "legend", start_row = 37,
      width = (plot_dim * 3/2), height = plot_dim, start_col = 10,
      plotname = "pre_pca", savedir = excel_basename,
      fancy = fancy)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }
    pre_table <- write_xlsx(
      wb, data = apr[["pre_batch"]][["table"]],
      sheet = "legend", title = "Pre-Batch PCA table.",
      start_row = 66, start_col = 10)
    xls_result <- write_xlsx(
      wb, data = apr[["post_batch"]][["table"]],
      sheet = "legend", title = "Pre-Batch PCA table.",
      start_row = pre_table[["end_row"]] + 2, start_col = 10)
  }
  retlist <- list(
    "image_files" = image_files,
    "table_names" = table_names,
    "xls_result" = xls_result)
  return(retlist)
}

#' Internal function to write a summary of some combined data
#'
#' @param wb xlsx workbook to which to write.
#' @param excel_basename basename for printing plots.
#' @param apr a pairwise result
#' @param extracted table extracted from the pairwise result
#' @param compare_plots series of plots to print out.
#' @param lfc_cutoff Used for volcano/MA plots.
#' @param p_cutoff Used for volcano/MA plots.
#' @param fancy Write fancy plots with the xlsx file?
write_combined_summary <- function(wb, excel_basename, apr, extracted, compare_plots,
                                   lfc_cutoff = 1, p_cutoff = 0.05, fancy = FALSE) {
  image_files <- c()
  xl_results <- c()
  if (length(apr[["comparison"]]) == 0) {
    compare_plots <- FALSE
  }
  if (isTRUE(compare_plots)) {
    sheetname <- "pairwise_summary"
    xls_result <- write_xlsx(
      wb, data = extracted[["summaries"]], sheet = sheetname,
      title = glue("Summary of contrasts (lfc cutoff:{lfc_cutoff} p cutoff: {p_cutoff})."))
    xl_results <- c(xl_results, xls_result)
    new_row <- xls_result[["end_row"]] + 2
    xls_result <- write_xlsx(
      wb, data = apr[["comparison"]][["comp"]], sheet = sheetname, start_row = new_row,
      title = "Pairwise correlation coefficients among differential expression tools.")
    xl_results <- c(xl_results, xls_result)
    new_row <- xls_result[["end_row"]] + 2
    if (class(apr[["comparison"]][["heat"]])[1] == "recordedplot") {
      try_result <- xlsx_insert_png(
        apr[["comparison"]][["heat"]], wb = wb, sheet = sheetname, plotname = "pairwise_summary",
        savedir = excel_basename, start_row = new_row + 1, start_col = 1, fancy = fancy)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
    }
    logfc_comparisons <- try(compare_logfc_plots(extracted), silent = TRUE)
    if (class(logfc_comparisons) != "try-error") {
      logfc_names <- names(logfc_comparisons)
      new_row <- new_row + 2
      for (c in seq_along(logfc_names)) {
        lname <- logfc_names[c]
        new_row <- new_row + 32
        le <- logfc_comparisons[[c]][["le"]]
        ld <- logfc_comparisons[[c]][["ld"]]
        de <- logfc_comparisons[[c]][["de"]]
        tmpcol <- 1
        if (!is.null(le)) {
          xls_result <- openxlsx::writeData(
            wb = wb, sheet = sheetname,
            startRow = new_row - 2, startCol = tmpcol,
            x = glue("Comparing DE tools for the \\
                                          comparison of: {logfc_names[c]}"))
          xl_results <- c(xl_results, xls_result)
          xls_result <- openxlsx::writeData(
            wb = wb, sheet = sheetname,
            startRow = new_row - 1, startCol = tmpcol,
            x="Log2FC(Limma vs. EdgeR)")
          xl_results <- c(xl_results, xls_result)
          try_result <- xlsx_insert_png(
            le, wb = wb, sheet = "pairwise_summary", plotname = "compare_le",
            savedir = excel_basename, start_row = new_row, start_col = tmpcol,
            fancy = fancy)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
          tmpcol <- 8
        }
        if (!is.null(ld)) {
          xls_result <- openxlsx::writeData(
            wb = wb, sheet = sheetname,
            startRow = new_row - 1, startCol = tmpcol,
            x = "Log2FC(Limma vs. DESeq2)")
          xl_results <- c(xl_results, xls_result)
          try_result <- xlsx_insert_png(
            ld, wb = wb, sheet = sheetname, plotname = "compare_ld", savedir = excel_basename,
            start_row = new_row, start_col = tmpcol, fancy = fancy)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
          tmpcol <- 15
        }
        if (!is.null(de)) {
          xls_result <- openxlsx::writeData(
            wb = wb, sheet = sheetname,
            startRow = new_row - 1, startCol = tmpcol,
            x = "Log2FC(DESeq2 vs. EdgeR)")
          xl_results <- c(xl_results, xls_result)
          try_result <- xlsx_insert_png(
            de, wb = wb, sheet = sheetname, plotname = "compare_ld", savedir = excel_basename,
            start_row = new_row, start_col = tmpcol, fancy = fancy)
          if (! "try-error" %in% class(try_result)) {
            image_files <- c(image_files, try_result[["filename"]])
          }
        }
      } ## End iterating over the comparison of logfc plots.
    } ## End checking if printing the logfc comparison plots worked.
  } ## End if compare_plots is TRUE
  retlist  <- list(
    "image_files" = image_files,
    "xl_results" = xl_results)
  return(retlist)
}

#' Writes out the results of a single pairwise comparison.
#'
#' However, this will do a couple of things to make one's life easier:
#' 1.  Make a list of the output, one element for each comparison of the
#'     contrast matrix.
#' 2.  Write out the results() output for them in separate sheets in excel.
#' 3.  Since I have been using qvalues a lot for other stuff, add a column.
#'
#' Tested in test_24deseq.R
#' Rewritten in 2016-12 looking to simplify combine_de_tables().  That function
#' is far too big, this should become a template for that.
#'
#' @param data Output from results().
#' @param type Which DE tool to write.
#' @param coef Add coefficients
#' @param table_type Contrasts or identities
#' @param excel Filename into which to save the xlsx data.
#' @param n Limit to the top n genes.
#' @param ... Parameters passed downstream, dumped into arglist and passed,
#'  notably the number of genes (n), the coefficient column (coef)
#' @return List of data frames comprising the toptable output for each
#'  coefficient, I also added a qvalue entry to these toptable() outputs.
#' @seealso \code{\link{write_xlsx}}
#' @examples
#' \dontrun{
#'  finished_comparison <- eBayes(deseq_output)
#'  data_list <- write_deseq(finished_comparison, workbook="excel/deseq_output.xls")
#' }
#' @export
write_de_table <- function(data, type = "limma", coef = NULL, table_type = "contrasts",
                           excel = "de_table.xlsx", n = 0, ...) {
  arglist <- list(...)
  if (!is.null(data[[type]])) {
    data <- data[[type]]
  }

  if (is.null(n)) {
    n <- 0
  }
  ## Figure out the number of genes if not provided
  if (n == 0) {
    n <- nrow(data[["coefficients"]])
  }

  if (is.null(coef)) {
    coef <- data[["contrasts_performed"]]
    if (table_type %in% coef) {
      coef <- data[["contrasts_performed"]][[table_type]]
    }
  } else {
    coef <- as.character(coef)
  }
  annot <- NULL
  if (!is.null(data[["input_data"]])) {

    annot <- rowData(data[["input_data"]])
  }

  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]

  return_data <- list()
  end <- length(coef)
  for (c in seq_len(end)) {
    comparison <- coef[c]
    mesg("Writing ", c, "/", end, ": table: ", comparison, ".")
    table <- data[["all_tables"]][[c]]
    if (!is.null(annot)) {
      table <- merge(annot, table, by = "row.names")
      rownames(table) <- table[["Row.names"]]
      table[["Row.names"]] <- NULL
    }


    written <- try(write_xlsx(
      data = table, wb = wb, sheet = comparison,
      title = glue("{type} results for: {comparison}.")))
  }

  save_result <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  return(save_result)
}

#' Internal function to write a legend for significant gene tables.
#'
#' @param wb xlsx workbook object from openxlsx.
write_sig_legend <- function(wb) {
  legend <- data.frame(rbind(
    c("The first ~3-10 columns of each sheet:",
      "are annotations provided by our chosen annotation source for this experiment."),
    c("Next 6 columns", "The logFC and p-values reported by limma, edger, and deseq2."),
    c("limma_logfc", "The log2 fold change reported by limma."),
    c("deseq_logfc", "The log2 fold change reported by DESeq2."),
    c("edger_logfc", "The log2 fold change reported by edgeR."),
    c("limma_adjp", "The adjusted-p value reported by limma."),
    c("deseq_adjp", "The adjusted-p value reported by DESeq2."),
    c("edger_adjp", "The adjusted-p value reported by edgeR."),
    c("The next 5 columns", "Statistics generated by limma."),
    c("limma_ave", "Average log2 expression observed by limma across all samples."),
    c("limma_t", "T-statistic reported by limma given the log2FC and variances."),
    c("limma_p", "Derived from limma_t, the p-value asking 'is this logfc significant?'"),
    c("limma_b", "Use a Bayesian estimate to calculate log-odds significance."),
    c("limma_q", "A q-value FDR adjustment of the p-value above."),
    c("The next 5 columns", "Statistics generated by DESeq2."),
    c("deseq_basemean", "Analagous to limma's ave column, the base mean of all samples."),
    c("deseq_lfcse", "The standard error observed given the log2 fold change."),
    c("deseq_stat", "T-statistic reported by DESeq2 given the log2FC and observed variances."),
    c("deseq_p", "Resulting p-value."),
    c("deseq_q", "False-positive corrected p-value."),
    c("deseq_num", "Numerator coefficients for this contrast."),
    c("deseq_den", "Denominator coefficients for this contrast."),
    c("The next 4 columns", "Statistics generated by edgeR."),
    c("edger_logcpm",
      "Similar DESeq2's basemean, only including the samples in the comparison."),
    c("edger_lr", "Undocumented, I think it is the T-statistic calculated by edgeR."),
    c("edger_p", "The observed p-value from edgeR."),
    c("edger_q", "The observed corrected p-value from edgeR."),
    c("The next 7 columns", "Statistics generated by ebseq."),
    c("ebseq_fc", "Fold-change reported by ebseq."),
    c("ebseq_logfc", "Oddly, ebseq also reports a log2fc."),
    c("ebseq_postfc", "Post-analysis fold change from ebseq."),
    c("ebseq_mean", "Mean values in the ebseq analysis."),
    c("ebseq_ppee", "Prior probability that the numerators and denominators are the same."),
    c("ebseq_ppde", "... that they are different (eg. 1-ppee)."),
    c("ebseq_adjp", "Currently just a copy of ppee until I figure them out."),
    c("The next 8 columns", "Statistics generated by the basic analysis."),
    c("basic_nummed", "log2 median values of the numerator (like edgeR's basemean)."),
    c("basic_denmed", "log2 median values of the denominator for this comparison."),
    c("basic_numvar", "Variance observed in the numerator values."),
    c("basic_denvar", "Variance observed in the denominator values."),
    c("basic_logfc", "The log2 fold change observed by the basic analysis."),
    c("basic_t", "T-statistic from basic."),
    c("basic_p", "Resulting p-value."),
    c("basic_adjp", "BH correction of the p-value."),
    c("The next 5 columns", "Summaries of the limma/deseq/edger results."),
    c("lfc_meta", "The mean fold-change value of limma/deseq/edger."),
    c("lfc_var", "The variance between limma/deseq/edger."),
    c("lfc_varbymed", "The ratio of the variance/median (lower means better agreement.)"),
    c("p_meta", "A meta-p-value of the mean p-values."),
    c("p_var", "Variance among the 3 p-values."),
    c("The last columns: top plot left",
      "Venn diagram of the genes with logFC > 0 and p-value <= 0.05 for limma/DESeq/Edger."),
    c("The last columns: top plot right",
      "Venn diagram of the genes with logFC < 0 and p-value <= 0.05 for limma/DESeq/Edger."),
    c("The last columns: second plot",
      "Scatter plot of the voom-adjusted/normalized counts for each coefficient."),
    c("The last columns: third plot",
      "Scatter plot of the adjusted/normalized counts for each coefficient from edgeR."),
    c("The last columns: fourth plot",
      "Scatter plot of the adjusted/normalized counts for each coefficient from DESeq."),
    c("", "If this data was adjusted with sva, look 'original_pvalues' at the end.")
  ),
  stringsAsFactors = FALSE)

  colnames(legend) <- c("column name", "column definition")
  xls_result <- write_xlsx(wb, data = legend, sheet = "legend", rownames = FALSE,
                           title = "Columns used in the following tables.")
  retlist <- list("wb" = wb, xls_result = xls_result)
  return(retlist)
}

#' Put the metadata at the end of combined_de_tables()
#'
#' For the moment this is a stupidly short function.  I am betting we will
#' elaborate on this over time.
#'
#' @param wb workbook object.
#' @param apr Pairwise result.
write_sample_design <- function(wb, apr) {
  meta_df <- colData(apr[["input"]])
  xls_meta_result <- write_xlsx(wb = wb, data = meta_df,
                                sheet = "metadata", title = "Experiment metadata.")
  return(xls_meta_result)
}

write_plots_de_xlsx <- function(includes, extracted, sheetname, current_row,
                                current_column, tab, xls_result, wb, plot_dim,
                                excel_basename, image_files,
                                plot_rows = 31, plot_columns = 10) {
  ## Now add the coefficients, ma, and volcanoes below the venns.
  ## Text on row 18, plots from 19-49 (30 rows)
  favorites <- c()
  if (isTRUE(includes[["deseq"]])) {
    favorites <- c(favorites, "deseq")
  }
  if (isTRUE(includes[["edger"]])) {
    favorites <- c(favorites, "edger")
  }
  if (isTRUE(includes[["limma"]])) {
    favorites <- c(favorites, "limma")
  }
  if (isTRUE(includes[["dream"]])) {
    favorites <- c(favorites, "dream")
  }
  if (isTRUE(includes[["ebseq"]])) {
    favorites <- c(favorites, "ebseq")
  }
  if (isTRUE(includes[["noiseq"]])) {
    favorites <- c(favorites, "noiseq")
  }
  if (isTRUE(includes[["basic"]])) {
    favorites <- c(favorites, "basic")
  }

  for (t in seq_along(favorites)) {
    num_plotted <- 0
    type <- favorites[t]
    sc <- paste0(type, "_scatter_plots")
    ma <- paste0(type, "_ma_plots")
    vo <- paste0(type, "_vol_plots")
    pp <- paste0(type, "_p_plots")
    ap <- paste0(type, "_adjp_plots")
    short <- substr(type, 0, 2)
    cap <- R.utils::capitalize(type)
    plt <- extracted[["plots"]][[sheetname]][[sc]]
    ma_plt <- extracted[["plots"]][[sheetname]][[ma]]
    vol_plt <- extracted[["plots"]][[sheetname]][[vo]]
    p_plt <- extracted[["plots"]][[sheetname]][[pp]]
    adjp_plt <- extracted[["plots"]][[sheetname]][[ap]]
    current_row <- current_row + 2
    current_column <- xls_result[["end_col"]] + 2

    ## Note that these are lists now.
    if (class(plt)[1] != "try-error" && length(plt) > 0) {
      ## The summary is now included with the plot.
      printme <- as.character(glue("{cap} expression coefficients for {tab}"))
      xl_result <- openxlsx::writeData(
        wb = wb, sheet = sheetname, x = printme,
        startRow = current_row, startCol = current_column)
      plotname <- paste0(short, "scatter")
      try_result <- xlsx_insert_png(
        plt[["scatter"]], wb = wb, sheet = sheetname,
        width = plot_dim, height = plot_dim, start_col = current_column,
        plotname = plotname, savedir = excel_basename, start_row = current_row + 1)
      if (!"try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
        num_plotted <- num_plotted + 1
        current_column <- current_column + plot_columns
      }
    }

    if (class(ma_plt)[1] != "try-error" && length(ma_plt) > 0) {
      xl_result <- openxlsx::writeData(
        wb = wb, sheet = sheetname, x = paste0(type, " MA plot"),
        startRow = current_row, startCol = current_column)
      plotname <- paste0(short, "ma")
      ## I think something here is failing and leaving behind a tempfile.
      try_ma_result <- xlsx_insert_png(
        ma_plt, wb = wb, sheet = sheetname, width = plot_dim,
        height = plot_dim, start_col = current_column, plotname = plotname,
        savedir = excel_basename, start_row = current_row + 1)
      if (!"try-error" %in% class(try_ma_result)) {
        image_files <- c(image_files, try_ma_result[["filename"]])
        num_plotted <- num_plotted + 1
        current_column <- current_column + plot_columns
      }
    }

    if (class(vol_plt)[1] != "try-error" && length(vol_plt) > 0) {
      xl_result <- openxlsx::writeData(
        wb = wb, sheet = sheetname, x = paste0(type, " volcano plot"),
        startRow = current_row, startCol = current_column)
      plotname <- paste0(short, "vol")
      try_vol_result <- xlsx_insert_png(
        vol_plt, wb = wb, sheet = sheetname, width = plot_dim,
        height = plot_dim, start_col = current_column, pltname = plotname,
        savedir = excel_basename, start_row = current_row + 1)
      if (!"try-error" %in% class(try_vol_result)) {
        image_files <- c(image_files, try_vol_result[["filename"]])
        num_plotted <- num_plotted + 1
        current_column <- current_column + plot_columns
      }
    }

    if (class(p_plt)[1] != "try-error" && length(p_plt) > 0) {
      xl_result <- openxlsx::writeData(
        wb = wb, sheet = sheetname, x = paste0(type, " p-value plot"),
        startRow = current_row, startCol = current_column)
      plotname <- paste0(short, "p")
      try_p_result <- xlsx_insert_png(
        p_plt, wb = wb, sheet = sheetname, width = plot_dim,
        height = plot_dim, start_col = current_column, pltname = plotname,
        savedir = excel_basename, start_row = current_row + 1)
      if (!"try-error" %in% class(try_p_result)) {
        image_files <- c(image_files, try_p_result[["filename"]])
        num_plotted <- num_plotted + 1
        current_column <- current_column + plot_columns
      }
    }

    if (class(adjp_plt)[1] != "try-error" && length(adjp_plt) > 0) {
      xl_result <- openxlsx::writeData(
        wb = wb, sheet = sheetname, x = paste0(type, " adjp-value plot"),
        startRow = current_row, startCol = current_column)
      plotname <- paste0(short, "adjp")
      try_adjp_result <- xlsx_insert_png(
        adjp_plt, wb = wb, sheet = sheetname, width = plot_dim,
        height = plot_dim, start_col = current_column, pltname = plotname,
        savedir = excel_basename, start_row = current_row + 1)
      if (!"try-error" %in% class(try_adjp_result)) {
        image_files <- c(image_files, try_adjp_result[["filename"]])
        num_plotted <- num_plotted + 1
        current_column <- current_column + plot_columns
      }
    }

    if (num_plotted > 0) {
      current_row <- current_row + plot_rows
    }

  } ## End adding limma, deseq, and edger plots.
  ret <- list(
    "image_files" = image_files,
    "wb" = wb,
    "current_row" = current_row,
    "current_column" = current_column)
  return(ret)
} ## End checking whether to add plots

write_venns_de_xlsx <- function(written_table, tab, wb, sheetname,
                                current_row, current_column, excel_basename,
                                plot_dim, image_files, venn_rows = 16,
                                venn_columns = 4, lfc_cutoff = 1.0, p_cutoff = 0.05,
                                includes = NULL, plot_columns = 10) {
  ## Make some venn diagrams comparing deseq/limma/edger!
  venns <- list()
  if (is.null(includes)) {
    includes <- list("deseq" = TRUE, "edger" = TRUE, "limma" = TRUE)
  }
  starting_column <- current_column
  venn_nop_lfc0 <- try(de_venn(written_table, lfc = 0, adjp = FALSE, p = 1.0))
  venn_nop <- try(de_venn(written_table, lfc = lfc_cutoff, adjp = FALSE, p = 1.0))
  venn_list <- try(de_venn(written_table, lfc = 0, adjp = p_cutoff))
  venn_sig_list <- try(de_venn(written_table, lfc = lfc_cutoff, adjp = p_cutoff))
  venns[[tab]] <- list(venn_nop_lfc0, venn_nop, venn_list, venn_sig_list)
  names(venns[[tab]]) <- c("nop_lfc0", "nop_lfc1", "p_lfc0", "p_lfc1")
  ## If they worked, add them to the excel sheets after the data,
  ## but make them smaller than other graphs.
  if (class(venn_list)[1] != "try-error") {
    ## First row of plots all going up
    xl_result <- openxlsx::writeData(
      wb = wb, sheet = sheetname, x = "Venn of all genes, lfc > 0.",
      startRow = current_row, startCol = current_column)
    up_plot <- venn_nop_lfc0[["up_venn"]]
    try_result <- xlsx_insert_png(
      up_plot, wb = wb, sheet = sheetname, width = (plot_dim / 2), height = (plot_dim / 2),
      start_col = current_column, plotname = "lfc0upvennnop", savedir = excel_basename,
      start_row = current_row + 1, doWeights = FALSE)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }
    current_column <- current_column + venn_columns
    venn_string <- paste0("Venn of all genes, lfc > ", lfc_cutoff, ".")
    xl_result <- openxlsx::writeData(
      wb = wb, sheet = sheetname, x = venn_string,
      startRow = 1, startCol = current_column)
    up_plot <- venn_nop[["up_venn"]]
    try_result <- xlsx_insert_png(
      up_plot, wb = wb, sheet = sheetname, width = (plot_dim / 2), height = (plot_dim / 2),
      start_col = current_column, plotname = "upvennnop", savedir = excel_basename,
      start_row = current_row + 1, doWeights = FALSE)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }
    current_column <- current_column + venn_columns
    xl_result <- openxlsx::writeData(
      wb = wb, sheet = sheetname, x = "Venn of p-value up genes, lfc > 0.",
      startRow = 1, startCol = current_column)
    up_plot <- venn_list[["up_venn"]]
    try_result <- xlsx_insert_png(
      up_plot, wb = wb, sheet = sheetname, width = (plot_dim / 2), height = (plot_dim / 2),
      start_col = current_column, plotname = "upvenn", savedir = excel_basename,
      start_row = current_row + 1, doWeights = FALSE)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }
    current_column <- current_column + venn_columns
    venn_string <- paste0("Venn of p-value up genes, lfc > ", lfc_cutoff, ".")
    xl_result <- openxlsx::writeData(
      wb = wb, sheet = sheetname, x = venn_string,
      startRow = current_row, startCol = current_column)
    sig_up_plot <- venn_sig_list[["up_venn"]]
    try_result <- xlsx_insert_png(
      sig_up_plot, wb = wb, sheet = sheetname, width = (plot_dim / 2),
      height =  (plot_dim / 2), start_col = current_column, plotname = "upvenn",
      savedir = excel_basename, start_row = current_row + 1, doWeights = FALSE)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }

    sig_methods <- c()
    if (isTRUE(includes[["limma"]])) {
      sig_methods <- c("limma", sig_methods)
    }
    if (isTRUE(includes[["edger"]])) {
      sig_methods <- c("edger", sig_methods)
    }
    if (isTRUE(includes[["deseq"]])) {
      sig_methods <- c("deseq", sig_methods)
    }
    siggene_lst <- try(plot_num_siggenes(written_table, methods =  sig_methods))
    current_column <- current_column + venn_columns
    if (class(siggene_lst)[1] != "try-error") {
      xl_result <- openxlsx::writeData(
        wb = wb, sheet = sheetname, x = "Significant genes by fc going up.",
        startRow = current_row, startCol = current_column)
      try_result <- xlsx_insert_png(
        siggene_lst[["up"]], wb = wb, sheet =  sheetname, width = plot_dim,
        height = (plot_dim / 2), start_col =  current_column, plotname =  "siggenesup",
        savedir = excel_basename, start_row = current_row + 1, doWeights = FALSE)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
      current_column <- current_column + plot_columns
      xl_result <- openxlsx::writeData(
        wb = wb, sheet = sheetname, x = "Significant genes by p going up.",
        startRow = 1, startCol = current_column)
      try_result <- xlsx_insert_png(
        siggene_lst[["pup"]], wb = wb, sheet = sheetname, width = plot_dim,
        height = (plot_dim / 2), start_col = current_column, plotname = "siggenespup",
        savedir = excel_basename, start_row = current_row + 1, doWeights = FALSE)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
    }
    ## Plot down venns etc, so reset the column but not the rows!
    current_column <- starting_column
    current_row <- current_row + venn_rows
    venn_string <- paste0("Venn of all genes, lfc < ", lfc_cutoff, ".")
    xl_result <- openxlsx::writeData(
      wb = wb, sheet = sheetname, x = venn_string,
      startRow = current_row, startCol = current_column)
    down_plot <- venn_nop_lfc0[["down_venn"]]
    try_result <- xlsx_insert_png(
      down_plot, wb = wb, sheet = sheetname, width = plot_dim / 2, height = plot_dim / 2,
      start_col = current_column, plotname = "lfc0downvennnop", savedir = excel_basename,
      start_row = current_row + 1, doWeights = FALSE)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }
    current_column <- current_column + venn_columns
    venn_string <- paste0("Venn of all genes, lfc < -", lfc_cutoff, ".")
    xl_result <- openxlsx::writeData(
      wb = wb, sheet = sheetname, x = venn_string,
      startRow = current_row, startCol = current_column)
    down_plot <- venn_nop[["down_venn"]]
    try_result <- xlsx_insert_png(
      down_plot, wb = wb, sheet = sheetname, width = plot_dim / 2, height = plot_dim / 2,
      start_col = current_column, plotname = "downvennnop", savedir = excel_basename,
      start_row = current_row + 1, doWeights = FALSE)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }
    current_column <- current_column + venn_columns
    xl_result <- openxlsx::writeData(
      wb = wb, sheet = sheetname, x = "Venn of p-value down genes, lfc < 0.",
      startRow = current_row, startCol = current_column)
    down_plot <- venn_list[["down_venn"]]
    try_result <- xlsx_insert_png(
      down_plot, wb = wb, sheet = sheetname, width = plot_dim / 2, height = plot_dim / 2,
      start_col = current_column, plotname = "downvenn", savedir = excel_basename,
      start_row = current_row + 1, doWeights = FALSE)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }
    current_column <- current_column + venn_columns
    xl_result <- openxlsx::writeData(
      wb = wb, sheet = sheetname, x = "Venn of p-value down genes, lfc < -1.",
      startRow = current_row, startCol = current_column)
    down_plot <- venn_sig_list[["down_venn"]]
    try_result <- xlsx_insert_png(
      down_plot, wb = wb, sheet = sheetname, width = plot_dim / 2, height = plot_dim / 2,
      start_col = current_column, plotname = "downvenn", savedir = excel_basename,
      start_row = current_row + 1, doWeights = FALSE)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }
    current_column <- current_column + venn_columns
    if (class(siggene_lst)[1] != "try-error") {
      xl_result <- openxlsx::writeData(
        wb = wb, sheet = sheetname, x = "Significant genes by fc going down.",
        startRow = current_row, startCol = current_column)
      try_result <- xlsx_insert_png(
        siggene_lst[["down"]], wb = wb, sheet = sheetname,
        width = plot_dim, height = (plot_dim / 2),
        start_col = current_column, plotname = "siggenesup",
        savedir = excel_basename, start_row = current_row + 1, doWeights = FALSE)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
      current_column <- current_column + plot_columns
      xl_result <- openxlsx::writeData(
        wb = wb, sheet = sheetname, x = "Significant genes by p going down.",
        startRow = current_row, startCol = current_column)
      try_result <- xlsx_insert_png(
        siggene_lst[["pdown"]], wb = wb, sheet = sheetname,
        width = plot_dim, height = (plot_dim / 2),
        start_col = current_column, plotname = "siggenespdown",
        savedir = excel_basename, start_row = current_row + 1, doWeights = FALSE)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
    }
    current_row <- current_row + venn_rows
  }
  ret <- list(
    "image_files" = image_files,
    "wb" = wb,
    "current_row" = current_row,
    "current_column" = current_column)
  return(ret)
}

summarize_combined <- function(comb, up_fc, down_fc, p_cutoff, adjp = TRUE) {
  limma_p_column <- "limma_p"
  deseq_p_column <- "deseq_p"
  edger_p_column <- "edger_p"
  basic_p_column <- "basic_p"
  ebseq_p_column <- "ebseq_p"
  if (isTRUE(adjp)) {
    limma_p_column <- "limma_adjp"
    deseq_p_column <- "deseq_adjp"
    edger_p_column <- "edger_adjp"
    basic_p_column <- "basic_adjp"
    ebseq_p_column <- "ebseq_adjp"
  }
  ret <- list(
    "total" = nrow(comb),
    "limma_up" = sum(comb[["limma_logfc"]] >= up_fc),
    "limma_sigup" = sum(
      comb[["limma_logfc"]] >= up_fc & as.numeric(comb[[limma_p_column]]) <= p_cutoff),
    "deseq_up" = sum(comb[["deseq_logfc"]] >= up_fc),
    "deseq_sigup" = sum(
      comb[["deseq_logfc"]] >= up_fc & as.numeric(comb[[deseq_p_column]]) <= p_cutoff),
    "edger_up" = sum(comb[["edger_logfc"]] >= up_fc),
    "edger_sigup" = sum(
      comb[["edger_logfc"]] >= up_fc & as.numeric(comb[[edger_p_column]]) <= p_cutoff),
    "basic_up" = sum(comb[["basic_logfc"]] >= up_fc),
    "basic_sigup" = sum(
      comb[["basic_logfc"]] >= up_fc & as.numeric(comb[["basic_p"]]) <= p_cutoff),
    "limma_down" = sum(comb[["limma_logfc"]] <= down_fc),
    "limma_sigdown" = sum(
      comb[["limma_logfc"]] <= down_fc & as.numeric(comb[[limma_p_column]]) <= p_cutoff),
    "deseq_down" = sum(comb[["deseq_logfc"]] <= down_fc),
    "deseq_sigdown" = sum(
      comb[["deseq_logfc"]] <= down_fc & as.numeric(comb[[deseq_p_column]]) <= p_cutoff),
    "edger_down" = sum(comb[["edger_logfc"]] <= down_fc),
    "edger_sigdown" = sum(
      comb[["edger_logfc"]] <= down_fc & as.numeric(comb[[edger_p_column]]) <= p_cutoff),
    "basic_down" = sum(comb[["basic_logfc"]] <= down_fc),
    "basic_sigdown" = sum(
      comb[["basic_logfc"]] <= down_fc & as.numeric(comb[["basic_p"]]) <= p_cutoff),
    "meta_up" = sum(comb[["fc_meta"]] >= up_fc),
    "meta_sigup" = sum(
      comb[["lfc_meta"]] >= up_fc & as.numeric(comb[["p_meta"]]) <= p_cutoff),
    "meta_down" = sum(comb[["lfc_meta"]] <= down_fc),
    "meta_sigdown" = sum(
      comb[["lfc_meta"]] <= down_fc & as.numeric(comb[["p_meta"]]) <= p_cutoff))
  return(ret)
}

write_upset_groups <- function(retlist, which = "all", excel = "excel/test.xlsx") {
  lU_letters <- c(letters, LETTERS)
  group_name <- paste0(which, "_groups")
  upsetr_groups <- retlist[[group_name]]
  list_name <- paste0(which, "_list")
  upsetr_list <- retlist[[list_name]]
  sig_name <- paste0(which, "_sig")
  sig_upsetr <- retlist[[sig_name]]
  plot_name <- paste0(which, "_plot")
  upsetr_plot <- retlist[[plot_name]]

  start_names <- names(upsetr_groups)
  categories <- unique(unlist(strsplit(x = start_names, split = ":")))
  names(categories) <- lU_letters[1:length(categories)]
  sheet_names <- start_names
  sheet_name_df <- data.frame(row.names = sheet_names)
  for (n in seq_along(categories)) {
    start <- as.character(categories[n])
    end <- names(categories)[n]
    sheet_names <- gsub(x = sheet_names, pattern = start, replacement = end)
  }
  sheet_name_df[["short_name"]] <- sheet_names
  sheet_name_df[["long_name"]] <- as.character(start_names)

  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]
  do_excel <- TRUE
  if (is.null(wb)) {
    do_excel <- FALSE
  }
  xlsx_result <- write_xlsx(wb, data = sheet_name_df, sheet = "legend", rownames = FALSE,
                            title = "Names of the sheets and their corresponding contrasts.")
  try_result <- xlsx_insert_png(start_row = 3, start_col = 12,
                                width = 10, height = 10,
                                a_plot = upsetr_plot, wb = wb, sheet = "legend")
  all_groups <- rownames(sheet_name_df)
  datum_subsets <- list()
  for (n in seq_along(all_groups)) {
    group_name <- all_groups[n]
    sheet_name <- sheet_names[n]
    group_pieces <- strsplit(x = group_name, split = ":")[[1]]
    first_group <- group_pieces[1]

    ## I think something is wrong here
    ids <- overlap_geneids(upsetr_groups, group_name)
    datum <- sig_upsetr[[first_group]]
    datum_subset <- datum[ids, ]
    datum_subsets[[group_name]] <- datum_subset
    xlsx_name <- gsub(x = sheet_name, pattern = ":", replacement = "_")
    xlsx_result <- write_xlsx(wb, data = datum_subset, sheet = xlsx_name,
                              title = paste0("Genes in group: ", group_name, "."))
  }
  save_result <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  return(datum_subsets)
}

## EOF
