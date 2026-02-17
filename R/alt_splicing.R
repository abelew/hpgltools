## A group of functions to examine the results from alternative splicing
## methods. The hope is to remove some corner cases from tools like
## suppa/miso/etc and somewhat standardize the resulting outputs.

gather_suppa_files <- function(file_prefix = NULL, numerator = "d15", denominator = "d10") {
  type_suffixes <- c("A3", "A5", "AF", "AL", "RI", "MX", "SE")
  if (is.null(file_prefix)) {
    file_prefix <- file.path("outputs", "90suppa_mm38_100_time")
  }
  comparison_prefix <- glue("{numerator}_{denominator}");
  dpsi_prefix <- file.path(file_prefix, "diff")
  tpm_prefix <- file.path(file_prefix, "tpm")
  event_prefix <- file.path(file_prefix, "events")
  psi_prefix <- file.path(file_prefix, "psi")

  type_dpsi_files <- file.path(dpsi_prefix, paste0(comparison_prefix, "_",
                                                   type_suffixes, ".dpsi"))
  tpm_files <- file.path(tpm_prefix, paste0(c(numerator, denominator), ".tpm"))
  type_events_files <- file.path(event_prefix, paste0("local_as_events_",
                                                     type_suffixes, "_strict.ioe"))
  type_psi_files <- file.path(dpsi_prefix, paste0(comparison_prefix, "_",
                                                  type_suffixes, ".psivec"))
  type_tpm_files <- file.path(tpm_prefix, paste0(c(numerator, denominator), ".tpm"))

  tx_dpsi_file <- file.path(dpsi_prefix, glue("{comparison_prefix}_ioi_empirical.dpsi"))
  tx_events_file <- file.path(event_prefix, "transcript_events.ioi")
  tx_psi_file <- file.path(dpsi_prefix, glue("{comparison_prefix}_ioi_classical.psivec"))
  tx_tpm_file <- file.path(dpsi_prefix, glue("{comparison_prefix}_ioi_empirical_avglogtpm.tab"))
  psi_numerator_matrices <- Sys.glob(file.path(psi_prefix, paste0(numerator, "*")))
  psi_denominator_matrices <- Sys.glob(file.path(psi_prefix, paste0(denominator, "*")))
  psi_matrix_files <- c(psi_numerator_matrices, psi_denominator_matrices)

  retlist <- list(
    "tx_dpsi" = tx_dpsi_file,
    "tx_events" = tx_events_file,
    "tx_psi" = tx_psi_file,
    "tx_tpm" = tx_tpm_file,
    "type_dpsi" = type_dpsi_files,
    "type_events" = type_events_files,
    "type_psi" = type_psi_files,
    "type_tpm" = type_tpm_files,
    "psi_matrices" = psi_matrix_files)
  return(retlist)
}

gather_rmats_files <- function(file_prefix = NULL, numerator = "t4h", denominator = "no") {
  type_suffixes <- c("A3", "A5", "AF", "AL", "RI", "MX", "SE")
  if (is.null(file_prefix)) {
    file_prefix <- file.path("outputs", "90rmats_hg38_111_infect_state")
  }
  retlist <- list(
    ## Taken from the GTF file
    "gtf_a3ss" = file.path(file_prefix, "fromGTF.A3SS.txt"),
    "gtf_a5ss" = file.path(file_prefix, "fromGTF.A5SS.txt"),
    "gtf_mex" = file.path(file_prefix, "fromGTF.MXE.txt"),
    "se_mex" = file.path(file_prefix, "fromGTF.SE.txt"),
    "ri_mex" = file.path(file_prefix, "fromGTF.RI.txt"),
    "gtf_a3_novel" = file.path(file_prefix, "fromGTF.novelJunction.A3SS.txt"),
    "gtf_a5_novel" = file.path(file_prefix, "fromGTF.novelJunction.A5SS.txt"),
    "gtf_mex_novel" = file.path(file_prefix, "fromGTF.novelJunction.MXE.txt"),
    "gtf_ri_novel" = file.path(file_prefix, "fromGTF.novelJunction.RI.txt"),
    "gtf_se_novel" = file.path(file_prefix, "fromGTF.novelJunction.SE.txt"),
    "gtf_a3_novelss" = file.path(file_prefix, "fromGTF.novelSpliceSite.A3SS.txt"),
    "gtf_a5_novelss" = file.path(file_prefix, "fromGTF.novelSpliceSite.A5SS.txt"),
    "gtf_mxe_novelss" = file.path(file_prefix, "fromGTF.novelSpliceSite.MXE.txt"),
    "gtf_ri_novelss" = file.path(file_prefix, "fromGTF.novelSpliceSite.RI.txt"),
    "gtf_se_novelss" = file.path(file_prefix, "fromGTF.novelSpliceSite.SE.txt"),
    ## Junction and exon reads
    "a5ss_both" = file.path(file_prefix, "A5SS.MATS.JCEC.txt"),
    "a3ss_both" = file.path(file_prefix, "A3SS.MATS.JCEC.txt"),
    "mxe_both" = file.path(file_prefix, "MXE.MATS.JCEC.txt"),
    "ri_both" = file.path(file_prefix, "RI.MATS.JCEC.txt"),
    "se_both" = file.path(file_prefix, "SE.MATS.JCEC.txt"),
    "a3_input_both" = file.path(file_prefix, "JCEC.raw.input.A3SS.txt"),
    "a5_input_both" = file.path(file_prefix, "JCEC.raw.input.A5SS.txt"),
    "mxe_input_both" = file.path(file_prefix, "JCEC.raw.input.MXE.txt"),
    "ri_input_both" = file.path(file_prefix, "JCEC.raw.input.RI.txt"),
    "se_input_both" = file.path(file_prefix, "JCEC.raw.input.SE.txt"),
    ## Final outputs, junction reads
    "a3ss_junct" = file.path(file_prefix, "A3SS.MATS.JC.txt"),
    "a5ss_junct" = file.path(file_prefix, "A5SS.MATS.JC.txt"),
    "mxe_junct" = file.path(file_prefix, "MXE.MATS.JC.txt"),
    "ri_junct" = file.path(file_prefix, "RI.MATS.JC.txt"),
    "se_junct" = file.path(file_prefix, "SE.MATS.JC.txt"),
    "a3_input_junct" = file.path(file_prefix, "JC.raw.input.A3SS.txt"),
    "a5_input_junct" = file.path(file_prefix, "JC.raw.input.A5SS.txt"),
    "mxe_input_junct" = file.path(file_prefix, "JC.raw.input.MXE.txt"),
    "ri_input_junct" = file.path(file_prefix, "JC.raw.input.RI.txt"),
    "se_input_junct" = file.path(file_prefix, "JC.raw.input.SE.txt"))
  return(retlist)
}


#' Given some psi and tpm data, make a pretty plot!
#'
#' This should take either a dataframe or filename for the psi data from suppa,
#' along with the same for the average log tpm data (acquired from suppa
#' diffSplice with --save_tpm_events)
#'
#' @param file_prefix Directory containing the various requisite input files.
#' @param file_list Vector of filenames.
#' @param type Either transcript or 'type' referring to the type of DPSI analysis performed.
#' @param annot Dataframe of annotations.
#' @param annot_column Column in annot with the IDs to cross reference against.
#' @param sig_threshold Use this significance threshold.
#' @param label_type Choose a type of event to label.
#' @param alpha How see-through should the points be in the plot?
#' @param numerator Name of the desired comparison's numerator.
#' @param denominator Name of the desired comparison's denominator.
#' @return List containing the plot and some of the requisite data.
#' @seealso [plot_rmats()]
#' @examples
#'  \dontrun{
#'  suppa_plot <- plot_suppa(dpsi_file, tmp_file)
#' }
#' @export
plot_suppa <- function(file_prefix, file_list = NULL, type = "type", annot = NULL,
                       annot_column = NULL, sig_threshold = 0.05, label_type = NULL, alpha = 0.3,
                       numerator = "infected", denominator = "uninfected") {
  ## A couple variable declarations to keep R CMD check happy when I use NSE syntax with dplyr
  transcripts_1 <- NULL
  event <- NULL

  if (!is.null(annot)) {
    if (! annot_column %in% colnames(annot)) {
      warning("The column ", annot_column,
              " is not in the colnames of the annotations, skipping.")
      annot <- NULL
    }
  }
  if (is.null(file_list)) {
    file_list <- gather_suppa_files(file_prefix = file_prefix,
                                    numerator = numerator,
                                    denominator = denominator)
  }
  dpsi <- NULL
  if (type == "type") {
    dpsi <- file_list[["type_dpsi"]]
  } else {
    dpsi <- file_list[["tx_dpsi"]]
  }
  events <- NULL
  if (type == "type") {
    events <- file_list[["type_events"]]
  } else {
    events <- file_list[["tx_events"]]
  }
  psi <- NULL
  if (type == "type") {
    psi <- file_list[["type_psi"]]
  } else {
    psi <- file_list[["tx_psi"]]
  }
  tpm <- NULL
  if (type == "type") {
    tpm <- file_list[["type_tpm"]]
  } else {
    tpm <- file_list[["tx_tpm"]]
  }

  dpsi_data <- data.frame()
  mesg("Merging dpsi data from ", length(dpsi), " files.")
  for (p in dpsi) {
    tmp_data <- read.table(p, sep = "\t", skip = 1)
    dpsi_data <- rbind(dpsi_data, tmp_data)
  }
  dpsi_data <- as.data.frame(dpsi_data)
  rownames(dpsi_data) <- dpsi_data[[1]]
  dpsi_data[[1]] <- NULL
  colnames(dpsi_data) <- c("dpsi", "pvalue")

  events_data <- data.frame()
  mesg("Merging events from ", length(events), " event files.")
  for (e in events) {
    tmp_data <- read.table(e, sep = "\t", skip = 1)
    events_data <- rbind(events_data, tmp_data)
  }
  events_data <- as.data.frame(events_data)
  colnames(events_data) <- c("number", "gene", "event", "transcripts_1", "transcripts_2")

  psi_data <- data.frame()
  mesg("Merging percentage-splice-included data from ", length(psi), " files.")
  fake_colnames <- c()
  for (p in psi) {
    tmp_data <- read.table(p, sep = "\t")
    if (is.null(fake_colnames[1])) {
      fake_colnames <- colnames(tmp_data)
    } else {
      colnames(tmp_data) <- fake_colnames
    }
    psi_data <- rbind(psi_data, tmp_data)
  }
  psi_data <- as.data.frame(psi_data)

  tpm_data <- data.frame()
  if (length(tpm) == 1) {
    tpm_data <- as.data.frame(read.table(tpm, sep = "\t"))
  } else {
    ## This provides the tpm values by transcript, I need to merge that
    ## against the set of events.
    count <- 0
    start_tpm <- data.frame()
    for (t in tpm) {
      count <- count + 1
      tmp_data <- read.table(t, sep = "\t")
      if (count == 1) {
        start_tpm <- tmp_data
      } else {
        start_tpm <- merge(start_tpm, tmp_data, by = "row.names")
        rownames(start_tpm) <- start_tpm[["Row.names"]]
        start_tpm[["Row.names"]] <- NULL
      }
    }
    mean_tpm_df <- data.frame(row.names = rownames(start_tpm))
    mean_tpm_df[["mean_tpm"]] <- log2(rowMeans(start_tpm + 1))
    events_by_tx <- events_data %>%
      tidyr::separate_rows(transcripts_1, sep = ",") %>%
      dplyr::select(transcripts_1, event)
    mean_tpm_df <- merge(mean_tpm_df, events_by_tx, by.x = "row.names", by.y = "transcripts_1")
    colnames(mean_tpm_df) <- c("transcript", "mean_tpm", "event")
    tpm_data <- mean_tpm_df[, c("event", "mean_tpm")]
    keep_events <- !duplicated(tpm_data[["event"]])
    tpm_data <- tpm_data[keep_events, ]
  }
  colnames(tpm_data) <- c("event", "avglogtpm")

  numerator_regex <- paste0("^", numerator)
  numerator_samples <- grepl(pattern = numerator_regex, x = colnames(psi_data))
  denominator_regex <- paste0("^", denominator)
  denominator_samples <- grepl(pattern = denominator_regex, x = colnames(psi_data))
  numerator_names <- c("event", paste0("numerator", seq_len(sum(numerator_samples))))
  colnames(psi_data)[numerator_samples] <- numerator_names
  denominator_names <- c("event", paste0("denominator", seq_len(sum(denominator_samples))))
  colnames(psi_data)[denominator_samples] <- denominator_names

  plotting_data <- merge(dpsi_data, tpm_data, by.x = "row.names", by.y = "event")
  plotting_data[["event"]] <- plotting_data[["Row.names"]]
  rownames(plotting_data) <- plotting_data[["Row.names"]]
  plotting_data[["Row.names"]] <- NULL

  ## Now we have the basis for everything we are likely to plot.
  ## Add some categorizations of the data.
  plotting_data[["dpsi"]] <- as.numeric(plotting_data[["dpsi"]])
  plotting_data[["log10pval"]] <- -1.0 * log10(as.numeric(plotting_data[["pvalue"]]))
  plotting_data[["adjp"]] <- p.adjust(plotting_data[["pvalue"]], method = "fdr")
  plotting_data[["psig"]] <- FALSE
  plotting_data[["adjpsig"]] <- FALSE
  plotting_data[["gene_name"]] <- gsub(x = rownames(plotting_data),
                                       pattern = "^(.*);.*$",
                                       replacement = "\\1")
  if (!is.null(annot)) {
    plotting_data <- merge(plotting_data, annot, by.x = "gene_name",
                           by.y = "row.names", all.x = TRUE)
    rownames(plotting_data) <- plotting_data[["event"]]
    plotting_data[["Row.names"]] <- NULL
  }
  if (is.null(annot_column)) {
    plotting_data[["label"]] <- plotting_data[["gene_name"]]
  } else {
    plotting_data[["label"]] <- plotting_data[[annot_column]]
  }

  plotting_data[["category"]] <- gsub(x = plotting_data[["event"]],
                                      pattern = "^(.*);(.{2}):.*$",
                                      replacement = "\\2")
  plotting_data[["coordinates"]] <- gsub(x = plotting_data[["event"]],
                                         pattern = "^(.*);(.{2}):(.*)$",
                                         replacement = "\\3")
  plotting_data[["plot_cat"]] <- plotting_data[["category"]]
  plotting_data[["plot_cat"]] <- ifelse(
    test = plotting_data[["plot_cat"]] == "SE",
    yes = "Skipping exon",
    no = ifelse(
      test = plotting_data[["plot_cat"]] == "MX",
      yes = "Mutually exclusive exons",
      no = ifelse(
        test = plotting_data[["plot_cat"]] == "A5",
        yes = "Alternate 5 prime",
        no = ifelse(
          test = plotting_data[["plot_cat"]] == "A3",
          yes = "Alternate 3 prime",
          no = ifelse(
            test = plotting_data[["plot_cat"]] == "RI",
            yes = "Retained intron",
            no = ifelse(plotting_data[["plot_cat"]] == "AF",
                        yes = "Alternate first exon",
                        no = ifelse(plotting_data[["plot_cat"]] == "AL",
                                    yes = "Alternate last exon",
                                    no = "Transcripts")))))))
  plotting_data[["category"]] <- plotting_data[["plot_cat"]]
  insig_idx <- plotting_data[["pvalue"]] > 0.05
  plotting_data[insig_idx, "plot_cat"] <- "Insignificant"
  ## If, somehow something is observed as unknown, make it insignificant.
  unknown_idx <- plotting_data[["plot_cat"]] == "Unknown"
  plotting_data[unknown_idx, "plot_cat"] <- "Insignificant"

  level_names <- c("Skipping exon", "Mutually exclusive exons", "Alternate 5 prime",
                   "Alternate 3 prime", "Retained intron", "Alternate first exon",
                   "Alternate last exon", "Insignificant", "Transcripts")
  plotting_data[["category"]] <- factor(plotting_data[["category"]], levels = level_names,
                                        labels = level_names)
  plotting_data[["event"]] <- rownames(plotting_data)

  psig_idx <- plotting_data[["pvalue"]] <= sig_threshold
  plotting_data[psig_idx, "psig"] <- TRUE
  adjpsig_idx <- plotting_data[["adjp"]] <= sig_threshold
  plotting_data[adjpsig_idx, "adjpsig"] <- TRUE

  ## A quick volcano plot, which should be made prettier soon.
  sig_splicing_volplot <- ggplot(plotting_data,
                                 aes(x = .data[["dpsi"]], y = .data[["log10pval"]],
                                     color = .data[["psig"]])) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = c(-0.5, 0.5), size = 1) +
    ggplot2::geom_hline(yintercept = 1.3, size = 1)

  ## Now a somewhat more involved ma plot, first dropping the super-low tpm stuff
  plotting_subset_idx <- plotting_data[["avglogtpm"]] >= 0
  plotting_subset <- plotting_data[plotting_subset_idx, ]
  label_subset_idx <- plotting_subset[["psig"]] == TRUE
  label_subset <- plotting_subset[label_subset_idx, ]
  if (!is.null(label_type)) {
    type_subset_idx <- label_subset[["plot_cat"]] == label_type
    label_subset <- label_subset[type_subset_idx, ]
  }

  color_values <- c("Skipping exon" = "#92250E",
                    "Mutually exclusive exons" = "#717600",
                    "Alternate 5 prime" = "#5B095A",
                    "Alternate 3 prime" = "#1D6E72",
                    "Retained intron" = "#31326D",
                    "Alternate first exon" = "black",
                    "Alternate last exon" = "#003300",
                    "Insignificant" = "#666666")

  ## Now make a quick and dirty ma plot.
  sig_splicing_maplot <- ggplot(plotting_subset,
                                aes(x = .data[["avglogtpm"]], y = .data[["dpsi"]],
                                    color = .data[["plot_cat"]], fill = .data[["plot_cat"]])) +
    ggplot2::geom_point(alpha = alpha) +
    ggplot2::scale_shape_manual(values = 21) +
    ggplot2::scale_fill_manual(name = "Category",
                               guide = "legend", ## keep my preferred order.
                               breaks = levels(plotting_data[["category"]]),
                               values = color_values) +
    ggplot2::scale_color_manual(name = "Category",
                                values = color_values, ## keep my preferred order.
                                breaks = levels(plotting_data[["category"]]),
                                guide = ggplot2::guide_legend(override.aes = list(size = 5))) +
    ggrepel::geom_text_repel(data = label_subset, alpha = 1.0,
                             show.legend = FALSE,
                             arrow = ggplot2::arrow(length = ggplot2::unit(0.01, "npc")),
                             aes(x = .data[["avglogtpm"]], y = .data[["dpsi"]],
                                 label = .data[["label"]])) +
    ggplot2::xlab("Average log(transcripts per million).") +
    ggplot2::ylab("Delta PSI calculated by Suppa.") +
    ggplot2::theme_bw(base_size = base_size)

  ## If available, add the event data to the table
  colnames(events_data) <- c("number", "gene", "event_id", "alternative_transcripts", "total_transcripts")
  if (!is.null(events_data)) {
    events_data <- events_data[, c("event_id", "alternative_transcripts", "total_transcripts")]
    plotting_data <- merge(plotting_data,
                           events_data,
                           by.x = "row.names",
                           by.y = "event_id")
    rownames(plotting_data) <- plotting_data[["Row.names"]]
    plotting_data[["Row.names"]] <- NULL
  }
  ## If available, add all the psi data to the table
  if (!is.null(psi_data)) {
    plotting_data <- merge(plotting_data,
                           psi_data,
                           by = "row.names")
    rownames(plotting_data) <- plotting_data[["Row.names"]]
    plotting_data[["Row.names"]] <- NULL
  }

  sig_idx <- plotting_data[["psig"]] == TRUE
  sig_data <- plotting_data[sig_idx, ]


  retlist <- list(
    "volcano" = sig_splicing_volplot,
    "ma" = sig_splicing_maplot,
    "data" = plotting_data,
    "sig" = sig_data,
    "files" = file_list)
  return(retlist)
}

#' Take a set of results from suppa and attempt to write it to a pretty xlsx file.
#'
#' Suppa provides a tremendous amount of output, this attempts to standardize
#' those results and print them to an excel sheet.
#'
#' @param table Result table from suppa.
#' @param annotations Set of annotation data to include with the suppa result.
#' @param by_table Use this column to merge the annotations and data tables from
#'  the perspective of the data table.
#' @param by_annot Use this column to merge the annotations and data tables
#'  from the perspective of the annotations.
#' @param columns Choose a subset of columns to include, or leave the defaults.
#' @param excel Provide an excel file to write.
#' @return Data frame of the merged data.
#' @seealso [write_xlsx()]
#' @examples
#'  \dontrun{
#'  prettier_table <- write_suppa_table(suppa_result_file,
#'                                      annotations = gene_info,
#'                                      excel = "excel/pretty_suppa_table.xlsx")
#' }
#' @export
write_suppa_table <- function(table, annotations = NULL, by_table = "gene_name",
                              by_annot = "ensembl_gene_id",
                              columns = "default", excel = "excel/suppa_table.xlsx") {

  default_columns <- c("event", "dpsi", "pvalue", "adjp", "avglogtpm", "category",
                       "coordinates", "alternative_transcripts", "total_transcripts",
                       "denominator1", "denominator2", "denominator3", "numerator1",
                       "numerator2", "numerator3", "numerator4", "numerator5", "numerator6")
  full_table <- data.frame()
  if (is.null(annotations)) {
    full_table <- data.table::as.data.table(table)
  } else {
    annot <- data.table::as.data.table(annotations)
    annot[[by_annot]] <- make.names(annot[[by_annot]], unique = TRUE)
    tab <- data.table::as.data.table(table)
    full_table <- merge(annot, tab, by.x = by_annot, by.y = by_table, all.y = TRUE)
  }

  xls_data <- as.data.frame(full_table)
  ## Now coerce numeric columns
  xlsx_result <- write_xlsx(data = xls_data, excel = excel)
  return(xls_data)
}

#' Given some psi and tpm data from rMATS, make a pretty plot!
#'
#' This should take either a dataframe or filename for the psi data from rMATS.
#' This was mostly copy/pasted from plot_suppa().
#'
#' @param se Table of skipped exon data from rmats.
#' @param a5ss Table of alternate 5p exons.
#' @param a3ss Table of alternate 3p exons.
#' @param mxe Table of alternate exons.
#' @param ri Table of retained introns.
#' @param sig_threshold Use this significance threshold.
#' @param dpsi_threshold Use a delta threshold.
#' @param label_type Choose a type of event to label.
#' @param alpha How see-through should the points be in the plot?
#' @return List containing the plot and some of the requisite data.
#' @seealso [plot_supps()]
#' @examples
#'  \dontrun{
#'  rmats_plot <- plot_rmats(se_table, a5_table, a3_table)
#' }
#' @export
plot_rmats <- function(se = NULL, a5ss = NULL, a3ss = NULL, mxe = NULL, ri = NULL,
                       sig_threshold = 0.05, dpsi_threshold = 0.7,
                       label_type = NULL, alpha = 0.7) {
  ijc_numerator <- ijc_denominator <- num_ijc_mean <- den_ijc_mean <- NULL ## R CMD check
  numerator_inclusion <- denominator_inclusion <- num_mean <- den_mean <- NULL ## Ibid
  if (is.null(se) && is.null(a5ss) && is.null(a3ss) &&
        is.null(mxe) && is.null(ri)) {
    stop("No data was provided.")
  }

  se_data <- data.frame()
  if (class(se) == "character") {
    se_data <- read.table(se, header = TRUE)
    mesg(se, " has ", nrow(se_data), " rows.")
  } else if (class(se) == "data.frame" | class(se) == "NULL") {
    se_data <- se_data
  } else {
    stop("I only understand filenames and data frames, your psi are neither.")
  }

  a5ss_data <- data.frame()
  if (class(a5ss) == "character") {
    a5ss_data <- read.table(a5ss, header = TRUE)
    mesg(a5ss, " has ", nrow(a5ss_data), " rows.")
  } else if (class(a5ss)[1] == "data.frame" || class(a5ss)[1] == "NULL") {
    a5ss_data <- a5ss_data
  } else {
    stop("I only understand filenames and data frames, your psi are neither.")
  }

  a3ss_data <- data.frame()
  if (class(a3ss) == "character") {
    a3ss_data <- read.table(a3ss, header = TRUE)
    mesg(a3ss, " has ", nrow(a3ss_data), " rows.")
  } else if (class(a3ss)[1] == "data.frame" || class(a3ss)[1] == "NULL") {
    a3ss_data <- a3ss_data
  } else {
    stop("I only understand filenames and data frames, your psi are neither.")
  }

  ## This explodes when using read.table.
  mxe_data <- data.frame()
  if (class(mxe) == "character") {
    mxe_data <- read.table(mxe, header = TRUE)
    mesg(mxe, " has ", nrow(mxe_data), " rows.")
  } else if (class(mxe)[1] == "data.frame" || class(mxe)[1] == "NULL") {
    mxe_data <- mxe_data
  } else {
    stop("I only understand filenames and data frames, your psi are neither.")
  }

  ri_data <- data.frame()
  if (class(ri) == "character") {
    ri_data <- read.table(ri, header = TRUE)
    mesg(ri, " has ", nrow(ri_data), " rows.")
  } else if (class(ri)[1] == "data.frame" || class(ri)[1] == "NULL") {
    ri_data <- ri_data
  } else {
    stop("I only understand filenames and data frames, your psi are neither.")
  }

  all_data <- data.frame()
  shared_columns <- c("ID", "GeneID", "geneSymbol", "chr", "strand", "IJC_SAMPLE_1",
                      "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncFormLen",
                      "SkipFormLen", "PValue", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference")
  new_colnames <- c("id", "gene_id", "gene_symbol", "chr", "strand", "ijc_numerator",
                    "sjc_numerator", "ijc_denominator", "sjc_denominator", "inclusion_length",
                    "skip_length", "p_value", "fdr", "numerator_inclusion", "denominator_inclusion",
                    "inclusion_difference", "event")
  if (nrow(se_data > 0)) {
    se_subset <- se_data[, shared_columns]
    se_subset[["event"]] <- "SE"
    colnames(se_subset) <- new_colnames
    all_data <- rbind(all_data, se_subset)
  }
  if (nrow(a5ss_data) > 0) {
    a5ss_subset <- a5ss_data[, shared_columns]
    a5ss_subset[["event"]] <- "A5"
    colnames(a5ss_subset) <- new_colnames
    all_data <- rbind(all_data, a5ss_subset)
  }
  if (nrow(a3ss_data) > 0) {
    a3ss_subset <- a3ss_data[, shared_columns]
    a3ss_subset[["event"]] <- "A3"
    colnames(a3ss_subset) <- new_colnames
    all_data <- rbind(all_data, a3ss_subset)
  }
  if (nrow(mxe_data) > 0) {
    mxe_subset <- mxe_data[, shared_columns]
    mxe_subset[["event"]] <- "MX"
    colnames(mxe_subset) <- new_colnames
    all_data <- rbind(all_data, mxe_subset)
  }
  if (nrow(ri_data) > 0) {
    ri_subset <- ri_data[, shared_columns]
    ri_subset[["event"]] <- "RI"
    colnames(ri_subset) <- new_colnames
    all_data <- rbind(all_data, ri_subset)
  }

  all_data <- data.table::as.data.table(all_data)
  all_data[["id"]] <- seq_len(nrow(all_data))
  ## Adding stub variables to make dplyr/tidyr NSE operations not throw warnings
  ## when running R CMD check
  l1a <- l1b <- l2a <- l2b <- id <- l1mean <- l2mean <- NULL
  mesg("Getting numerator/denominator mean values, this is slow.")
  plotting_data <- all_data %>%
    group_by(id) %>%
    dplyr::mutate(
      num_ijc_mean = suppressWarnings(mean(as.numeric(strsplit(ijc_numerator, ",")[[1]]), na.rm = TRUE)),
      den_ijc_mean = suppressWarnings(mean(as.numeric(strsplit(ijc_denominator, ",")[[1]]), na.rm = TRUE)),
      all_ijc_mean = mean(c(num_ijc_mean, den_ijc_mean), na.rm = TRUE),
      num_mean = suppressWarnings(mean(as.numeric(strsplit(numerator_inclusion, ",")[[1]]), na.rm = TRUE)),
      den_mean = suppressWarnings(mean(as.numeric(strsplit(denominator_inclusion, ",")[[1]]), na.rm = TRUE)),
      all_mean = mean(c(num_mean, den_mean), na.rm = TRUE))
  mesg("Adding comparison columns and log p-values.")
  plotting_data[["subtraction"]] <- plotting_data[["num_mean"]] - plotting_data[["den_mean"]]
  plotting_data[["log10pval"]] <- -1.0 * log10(plotting_data[["p_value"]])
  plotting_data[["log10adjpval"]] <- -1.0 * log10(plotting_data[["fdr"]])
  plotting_data[["log2ijc"]] <- log2(plotting_data[["all_ijc_mean"]] + 1)
  plotting_data[["psig"]] <- FALSE
  plotting_data[["adjpsig"]] <- FALSE
  plotting_data[["plot_cat"]] <- plotting_data[["event"]]
  plotting_data[["plot_cat"]] <- ifelse(
    test = plotting_data[["plot_cat"]] == "SE",
    yes = "Skipping exon",
    no = ifelse(
      test = plotting_data[["plot_cat"]] == "MX",
      yes = "Mutually exclusive exons",
      no = ifelse(
        test = plotting_data[["plot_cat"]] == "A5",
        yes = "Alternate 5 prime",
        no = ifelse(
          test = plotting_data[["plot_cat"]] == "A3",
          yes = "Alternate 3 prime",
          no = ifelse(
            test = plotting_data[["plot_cat"]] == "RI",
            yes = "Retained intron",
            no = ifelse(plotting_data[["plot_cat"]] == "AF",
                        yes = "Alternate first exon",
                        no = ifelse(plotting_data[["plot_cat"]] == "AL",
                                    yes = "Alternate last exon",
                                    no = "Unknown")))))))
  plotting_data[["category"]] <- plotting_data[["plot_cat"]]
  insig_idx <- plotting_data[["fdr"]] > sig_threshold
  plotting_data[insig_idx, "plot_cat"] <- "Insignificant"
  ## If, somehow something is observed as unknown, make it insignificant.
  unknown_idx <- plotting_data[["plot_cat"]] == "Unknown"
  plotting_data[unknown_idx, "plot_cat"] <- "Insignificant"

  level_names <- c("Skipping exon", "Mutually exclusive exons", "Alternate 5 prime",
                   "Alternate 3 prime", "Retained intron", "Alternate first exon",
                   "Alternate last exon", "Insignificant")
  plotting_data[["category"]] <- factor(plotting_data[["category"]], levels = level_names,
                                        labels = level_names)
  plotting_data[["event"]] <- rownames(plotting_data)

  psig_idx <- plotting_data[["fdr"]] <= sig_threshold
  plotting_data[psig_idx, "psig"] <- TRUE
  adjpsig_idx <- plotting_data[["fdr"]] <= sig_threshold
  plotting_data[adjpsig_idx, "adjpsig"] <- TRUE

  ## A quick volcano plot, which should be made prettier soon.
  sig_splicing_volplot <- ggplot(plotting_data,
                                 aes(x = .data[["dpsi"]], y = .data[["log10pval"]],
                                     color = .data[["psig"]])) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = c(-0.5, 0.5), linewidth = 1) +
    ggplot2::geom_hline(yintercept = 1.3, linewidth = 1)

  ## Now a somewhat more involved ma plot, first dropping the super-low tpm stuff
  label_subset_idx <- plotting_data[["psig"]] == TRUE
  label_subset <- plotting_data[label_subset_idx, ]
  label_subset_idx <- abs(label_subset[["subtraction"]]) > dpsi_threshold ## |
  ##    abs(label_subset[["all_mean"]] - 0.5) > (dpsi_threshold / 2)
  label_subset <- label_subset[label_subset_idx, ]
  ##label_subset_idx <- label_subset[["plot_cat"]] != "Insignificant"
  ##label_subset <- label_subset[label_subset_idx, ]
  if (!is.null(label_type)) {
    type_subset_idx <- label_subset[["plot_cat"]] == label_type
    label_subset <- label_subset[type_subset_idx, ]
  }

  color_values <- c("Skipping exon" = "#92250E",
                    "Mutually exclusive exons" = "#717600",
                    "Alternate 5 prime" = "#5B095A",
                    "Alternate 3 prime" = "#1D6E72",
                    "Retained intron" = "#31326D",
                    "Alternate first exon" = "black",
                    "Alternate last exon" = "#003300",
                    "Insignificant" = "#666666")

  ## Now make a quick and dirty ma plot.
  sig_splicing_maplot <- ggplot(
    plotting_data,
    aes(x = .data[["log2ijc"]], y = .data[["subtraction"]],
        color = .data[["plot_cat"]], fill = .data[["plot_cat"]])) +
    ggplot2::geom_point(alpha = alpha) +
    ggplot2::scale_shape_manual(values = 21) +
    ggplot2::scale_fill_manual(name = "Category",
                               guide = "legend", ## keep my preferred order.
                               breaks = levels(plotting_data[["category"]]),
                               values = color_values) +
    ggplot2::scale_color_manual(name = "Category",
                                values = color_values, ## keep my preferred order.
                                breaks = levels(plotting_data[["category"]]),
                                guide = ggplot2::guide_legend(override.aes = list(size = 5))) +
    ggrepel::geom_text_repel(
      data = label_subset,
      show.legend = FALSE,
      arrow = ggplot2::arrow(length = ggplot2::unit(0.01, "npc")),
      aes(x = .data[["all_mean"]], y = .data[["subtraction"]], label = .data[["gene_id"]])) +
    ggplot2::xlab("Average Inclusion.") +
    ggplot2::ylab("Delta PSI calculated by rMATS.") +
    ggplot2::theme_bw(base_size = base_size)

  retlist <- list(
    "volcano" = sig_splicing_volplot,
    "ma" = sig_splicing_maplot,
    "data" = plotting_data)
  return(retlist)
}

## EOF
