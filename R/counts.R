#' Sum the reads/gene for multiple sequencing runs of a single condition/batch.
#'
#' On occasion we have multiple technical replicates of a sequencing run.  This
#' can use a column in the experimental design to identify those replicates and
#' sum the counts into a single column in the count tables.
#'
#' Untested as of 2016-12-01, but used in a couple of projects where sequencing
#' runs got repeated.
#'
#' @param expt Experiment class containing the requisite metadata and count tables.
#' @param column Column of the design matrix used to specify which samples are replicates.
#' @return Expt with the concatenated counts, new design matrix, batches, conditions, etc.
#' @seealso [Biobase] [exprs()] [fData()] [pData()] [create_expt()]
#' @examples
#' \dontrun{
#'  compressed <- concatenate_runs(expt)
#' }
#' @export
concatenate_runs <- function(expt, column = "replicate") {
  design <- pData(expt)
  message("The original expressionset has ", nrow(design), " samples.")
  replicates <- levels(as.factor(design[[column]]))
  final_expt <- expt
  final_data <- NULL
  final_design <- NULL
  column_names <- list()
  colors <- list()
  conditions <- list()
  batches <- list()
  samplenames <- list()
  for (rep in replicates) {
    ## expression <- paste0(column, "=='", rep, "'")
    expression <- glue("{column} == '{rep}'")
    tmp_expt <- subset_expt(expt, expression)
    tmp_data <- rowSums(exprs(tmp_expt))
    tmp_design <- pData(tmp_expt)[1, ]
    final_data <- cbind(final_data, tmp_data)
    final_design <- rbind(final_design, tmp_design)
    column_names[[rep]] <- as.character(tmp_design[, "sampleid"])
    colors[[rep]] <- as.character(tmp_expt[["colors"]][1])
    batches[[rep]] <- as.character(pData(tmp_expt)[["batch"]][1])
    conditions[[rep]] <- as.character(pData(tmp_expt)[["condition"]][1])
    samplenames[[rep]] <- paste(conditions[[rep]], batches[[rep]], sep = "-")
    colnames(final_data) <- column_names
  }
  metadata <- new("AnnotatedDataFrame", final_design)
  sampleNames(metadata) <- colnames(final_data)
  feature_data <- new("AnnotatedDataFrame", fData(expt))
  featureNames(feature_data) <- rownames(final_data)
  experiment <- new("ExpressionSet", exprs = final_data,
                    phenoData = metadata, featureData = feature_data)
  final_expt[["expressionset"]] <- experiment
  final_expt[["samples"]] <- final_design
  final_expt[["colors"]] <- as.character(colors)
  final_expt[["samplenames"]] <- as.character(samplenames)
  message("The final expressionset has ", nrow(pData(final_expt)), " samples.")
  return(final_expt)
}
setGeneric("concatenate_runs")

#' Sum the reads/gene for multiple sequencing runs of a single condition/batch.
#'
#' On occasion we have multiple technical replicates of a sequencing run.  This
#' can use a column in the experimental design to identify those replicates and
#' sum the counts into a single column in the count tables.
#'
#' Untested as of 2016-12-01, but used in a couple of projects where sequencing
#' runs got repeated.
#'
#' @param expt Experiment class containing the requisite metadata and count tables.
#' @param column Column of the design matrix used to specify which samples are replicates.
#' @return Expt with the concatenated counts, new design matrix, batches, conditions, etc.
#' @seealso [Biobase] [exprs()] [fData()] [pData()] [create_expt()]
#' @examples
#' \dontrun{
#'  compressed <- concatenate_runs(expt)
#' }
#' @export
setMethod(
  "concatenate_runs", signature(expt = "SummarizedExperiment", column = "character"),
  definition = function(expt, column = "replicate") {
    design <- colData(expt)
    message("The original SE has ", nrow(design), " samples.")
    replicates <- levels(as.factor(design[[column]]))
    final_expt <- expt
    final_data <- NULL
    final_design <- NULL
    column_names <- list()
    colors <- list()
    conditions <- list()
    batches <- list()
    samplenames <- list()
    for (rep in replicates) {
      idx <- design[[column]] == rep
      sub <- expt[, idx]
      tmp_data <- rowSums(assay(sub))
      tmp_design <- colData(sub)[1, ]
      final_data <- cbind(final_data, tmp_data)
      final_design <- rbind(final_design, tmp_design)
      column_names[[rep]] <- as.character(tmp_design[, "sampleid"])
      colors[[rep]] <- as.character(get_colors(sub)[1])
      batches[[rep]] <- as.character(batches(sub)[1])
      conditions[[rep]] <- as.character(conditions(sub)[1])
      samplenames[[rep]] <- paste(conditions[[rep]], batches[[rep]], sep = "-")
      colnames(final_data) <- column_names
    }

    metadata <- new("AnnotatedDataFrame", as.data.frame(final_design))
    sampleNames(metadata) <- colnames(final_data)
    feature_data <- new("AnnotatedDataFrame", as.data.frame(rowData(expt)))
    featureNames(feature_data) <- rownames(final_data)
    experiment <- SummarizedExperiment(assays = final_data,
                                       colData = as.data.frame(final_design),
                                       rowData = as.data.frame(rowData(expt)))
    colors(experiment) <- colors
    message("The final SE has ", nrow(colData(experiment)), " samples.")
    return(experiment)
  })

#' Read a bunch of count tables and create a usable data frame from them.
#'
#' It is worth noting that this function has some logic intended for the elsayed
#' lab's data storage structure. It shouldn't interfere with other usages, but
#' it attempts to take into account different ways the data might be stored.
#'
#' Used primarily in create_expt()
#' This is responsible for reading count tables given a list of filenames.  It
#' tries to take into account upper/lowercase filenames and uses data.table to
#' speed things along.
#'
#' @param ids List of experimental ids.
#' @param files List of files to read.
#' @param header Whether or not the count tables include a header row.
#' @param include_summary_rows Whether HTSeq summary rows should be included.
#' @param all.x When merging (as opposed to join), choose the x data column.
#' @param all.y When merging (as opposed to join), choose the y data column.
#' @param merge_type Choose one, merge or join.
#' @param suffix Optional suffix to add to the filenames when reading them.
#' @param countdir Optional count directory to read from.
#' @param tx_gene_map Dataframe which provides a mapping between
#'  transcript IDs and gene IDs.
#' @param file_type Short circuit the file format autodetection.
#' @param ignore_tx_version Pass along TRUE to tximport's parameter
#'  ignoreTxIds to alleviate the headaches associated with salmon's
#'  stupid transcript ID .x suffix.
#' @param ... More options for happy time!
#' @return Data frame of count tables.
#' @seealso [data.table] [create_expt()] [tximport]
#' @examples
#' \dontrun{
#'  count_tables <- read_counts(as.character(sample_ids), as.character(count_filenames))
#' }
#' @importFrom data.table as.data.table
#' @importFrom S4Vectors mcols
#' @export
read_counts <- function(ids, files, header = FALSE, include_summary_rows = FALSE,
                        all.x = TRUE, all.y = FALSE, merge_type = "merge",
                        suffix = NULL, countdir = NULL, tx_gene_map = NULL,
                        file_type = NULL, ignore_tx_version = TRUE, ...) {
  ## load first sample
  arglist <- list(...)
  retlist <- list()
  ## Add an optional directory if I don't feel like specifying in the sample sheet.
  if (!is.null(countdir)) {
    files <- file.path(countdir, files)
  }
  skippers <- (files == "" | files == "undef" | is.null(files) | is.na(files))
  files <- files[!skippers]
  ids <- ids[!skippers]
  skippers <- (ids == "" | ids == "undef" | is.null(ids) | is.na(ids))
  if (sum(skippers) > 0) {
    message("Checking for NULL/undefined/NA filenames resulted in skipping ",
            sum(skippers), " files.")
    files <- files[!skippers]
    ids <- ids [!skippers]
  }

  retlist[["kept_ids"]] <- ids
  retlist[["kept_files"]] <- files
  ## lower_filenames <- files
  dirs <- dirname(files)

  count_table <- NULL
  for (f in seq_along(files)) {
    ## Get rid of any lurking spaces
    files[f] <- gsub(pattern = " ", replacement = "", x = files[f])
    ## Check if these are relative or absolute paths and standardize them to absolute.
    if (!grepl(pattern = "^\\/", x = files[f])) {
      files[f] <- file.path(getwd(), files[f])
    }
  }

  ## When I start using sailfish/salmon/etc, I will need to add more conditions to this
  ## and probably change it to a switch to more prettily take them into account.

  ## Likely important options:
  ## txOut: When true, do not back-convert to gene-level.
  ## Otherwise, tximport requires a tx2gene data frame with 2 columns, TXNAME and GENEID.
  ## These columns should be definition be available in my fData annotation.
  ## Therefore, I will set the flags tx2gene and txOut accordingly.
  mesg("Reading count tables.")
  txout <- TRUE
  if (!is.null(tx_gene_map)) {
    mesg("Using the transcript to gene mapping.")
    message("In some cases, (notably salmon) the format of the IDs used by this can be tricky.
It is likely to require the transcript ID followed by a '.' and the ensembl column:
'transcript_version', which is explicitly different than the gene version column.
If this is not correctly performed, very few genes will be observed")
    txout <- FALSE
  }

  if (is.null(file_type)) {
    if (grepl(pattern = "\\.tsv|\\.h5", x = files[1])) {
      file_type <- "kallisto"
    } else if (grepl(pattern = "\\.genes\\.results", x = files[1])) {
      file_type <- "rsem"
    } else if (grepl(pattern = "\\.sf", x = files[1])) {
      file_type <- "salmon"
    } else if (grepl(pattern = "fcounts\\.csv", x = files[1])) {
      file_type <- "featureCounts"
    } else if (grepl(pattern = "\\.bw\\.*.*", x = files[1])) {
      file_type = "bigwig"
    } else {
      file_type <- "table"
    }
  }

  if (file_type == "kallisto") {
    mesg("Reading kallisto data with tximport.")
    ## This hits if we are using the kallisto outputs.
    names(files) <- ids
    if (!all(file.exists(files))) {
      warning(files)
    }
    import <- NULL
    import_scaled <- NULL
    if (is.null(tx_gene_map)) {
      warning("Check that these count column names are correct, it may be the case that tximport failed and we need to set the colnames as per salmon.")
      import <- sm(tximport::tximport(files = files, type = "kallisto", txOut = txout))
      import_scaled <- sm(tximport::tximport(
        files = files, type = "kallisto",
        txOut = txout, countsFromAbundance = "lengthScaledTPM"))
    } else {
      warning("Check that these count column names are correct, it may be the case that tximport failed and we need to set the colnames as per salmon.")
      import <- sm(tximport::tximport(
        files = files, type = "kallisto", tx2gene = tx_gene_map, txOut = txout))
      import_scaled <- sm(tximport::tximport(
        files = files, type = "kallisto", tx2gene = tx_gene_map,
        txOut = txout, countsFromAbundance = "lengthScaledTPM"))
    }
    retlist[["count_table"]] <- as.data.table(
      import[["counts"]], keep.rownames = "rownames")
    retlist[["count_table"]] <- setkey(retlist[["count_table"]], rownames)
    retlist[["tximport"]] <- import
    retlist[["tximport_scaled"]] <- import_scaled
    retlist[["source"]] <- "tximport"
  } else if (file_type == "rsem") {
    mesg("Reading rsem data with tximport.")
    names(files) <- ids
    import <- NULL
    import_scaled <- NULL
    if (is.null(tx_gene_map)) {
      warning("Check that these count column names are correct, it may be the case that tximport failed and we need to set the colnames as per salmon.")
      import <- tximport::tximport(files = files, type = "rsem", txOut = txout)
      import_scaled <- tximport::tximport(files = files, type = "rsem",
                                          txOut = txout, countsFromAbundance = "lengthScaledTPM")
    } else {
      import <- tximport::tximport(files = files, type = "rsem",
                                   tx2gene = tx_gene_map, txOut = txout)
      import_scaled <- tximport::tximport(files = files, type = "rsem", tx2gene = tx_gene_map,
                                          txOut = txout, countsFromAbundance = "lengthScaledTPM")
    }
    retlist[["count_table"]] <- as.data.table(import[["counts"]],
                                                          keep.rownames = "rownames")
    retlist[["tximport"]] <- import
    retlist[["tximport_scaled"]] <- import_scaled
    retlist[["source"]] <- "tximport"
  } else if (file_type == "salmon") {
    mesg("Reading salmon data with tximport.")
    ## This hits if we are using the salmon outputs.
    names(files) <- ids
    missing_idx <- !file.exists(files)
    if (sum(missing_idx) > 0) {
      missing_files <- files[missing_idx]
      warning("Not all the files exist: ", toString(missing_files))
    }
    import <- NULL
    import_scaled <- NULL
    if (is.null(tx_gene_map)) {
      import <- sm(tximport::tximport(files = files, type = "salmon",
                                      txOut = txout, ignoreTxVersion = ignore_tx_version))
      import_scaled <- sm(tximport::tximport(
        files = files, type = "salmon",
        txOut = txout, countsFromAbundance = "lengthScaledTPM",
        ignoreTxVersion = ignore_tx_version))
    } else {
      ## Add a little test to see how well the tx_gene_map versions
      ## match those in the results from salmon.
      if (isTRUE(ignore_tx_version)) {
        message("Rewriting the transcript<->gene map to remove tx versions.")
        tx_gene_map[["novers"]] = gsub(x = tx_gene_map[["transcript"]], pattern = "^(.*)\\..*$", replacement = "\\1")
        tx_gene_map <- tx_gene_map[, c("novers", "ensembl_gene_id")]
        colnames(tx_gene_map) <- c("transcript", "gene")
      }
      import <- tximport::tximport(
        files = as.character(files), type = "salmon",
        tx2gene = tx_gene_map, txOut = txout,
        ignoreTxVersion = ignore_tx_version)
      ## It appears that something in tximport recently changed which causes it to no longer
      ## put the actual sampleIDs as the colnames of its return.
      ## And, since I use the count table columns as the primary arbiter of the sample IDs
      ## when merging the sample annotations (metadata) and the counts, this does not end well.
      colnames(import[["abundance"]]) <- ids
      colnames(import[["counts"]]) <- ids
      colnames(import[["length"]]) <- ids
      import_scaled <- sm(tximport::tximport(
        files = files, type = "salmon", tx2gene = tx_gene_map,
        txOut = txout, countsFromAbundance = "lengthScaledTPM",
        ignoreTxVersion = ignore_tx_version))
      colnames(import_scaled[["abundance"]]) <- ids
      colnames(import_scaled[["counts"]]) <- ids
      colnames(import_scaled[["length"]]) <- ids
    }
    retlist[["count_table"]] <- as.data.table(
      import[["counts"]], keep.rownames = "rownames")
    retlist[["count_table"]] <- setkey(retlist[["count_table"]], rownames)
    retlist[["tximport"]] <- import
    retlist[["tximport_scaled"]] <- import_scaled
    retlist[["source"]] <- "tximport"
  } else if (file_type == "featureCounts") {
    mesg("Reading featureCounts with read.table().")
    ## Use this codepath when we are working with htseq
    count_table <- read.table(files[1], header = header, stringsAsFactors = FALSE, skip = 2)
    count_table <- count_table[, c(1, 7)]
    colnames(count_table) <- c("rownames", ids[1])
    ## We are going to immediately check for NA, so I think we can suppress warnings.
    count_table[, 2] <- suppressWarnings(as.numeric(count_table[, 2]))
    na_idx <- is.na(count_table[[2]])
    ## This is a bit more circuituous than I would like.
    ## I want to make sure that na_rownames does not evaluate to something
    ## annoying like 'logical(0)' which it will if there are no nas in the count
    ## table and you just ask for is.na(count_table).
    na_rownames <- ""
    if (sum(na_idx) > 0) {
      na_rownames <- count_table[na_idx, "rownames"]
    }
    keepers_idx <- count_table[["rownames"]] != na_rownames
    count_table <- count_table[keepers_idx, ]
    count_table <- as.data.table(count_table)
    count_table <- setkey(count_table, rownames)
    first_rownames <- sort(count_table[["rownames"]])
    if (class(count_table)[1] == "try-error") {
      stop("There was an error reading: ", files[1])
    }
    mesg(files[1], " contains ", length(rownames(count_table)), " rows.")
    ## iterate over and append remaining samples
    for (num in seq(from = 2, to = length(files))) {
      table <- files[num]
      if (file.exists(tolower(table))) {
        table <- tolower(table)
      }
      tmp_count <- try(read.table(table, header = header, skip = 2))
      tmp_count <- tmp_count[, c(1, 7)]
      ## Drop the rows with NAs before coercing to numeric.
      keepers_idx <- tmp_count[[1]] != na_rownames
      tmp_count <- tmp_count[keepers_idx, ]
      current_rownames <- sort(tmp_count[[1]])
      mismatched_rownames <- sum(first_rownames != current_rownames)
      if (mismatched_rownames > 0) {
        warning("The file: ", table, " has mismatched rownames.")
      }
      tmp_count[, 2] <- as.numeric(tmp_count[, 2])
      if (class(tmp_count)[1] == "try-error") {
        stop("There was an error reading: ", table)
      }
      colnames(tmp_count) <- c("rownames", ids[num])
      tmp_count <- as.data.table(tmp_count)
      pre_merge <- nrow(tmp_count)
      if (merge_type == "merge") {
        count_table <- merge(count_table, tmp_count, by = "rownames",
                             all.x = all.x, all.y = all.y)
      } else {
        count_table <- plyr::join(count_table, tmp_count)
      }
      ## rownames(count_table) <- count_table[, "Row.names"]
      ## count_table <- count_table[, -1, drop = FALSE]
      ## post_merge <- length(rownames(count_table))
      post_merge <- nrow(count_table)
      mesg(table, " contains ", pre_merge, " rows and merges to ", post_merge, " rows.")
    } ## End for loop
    retlist[["count_table"]] <- count_table
    retlist[["source"]] <- "featureCounts"
    retlist[["count_table"]] <- setkey(retlist[["count_table"]], rownames)
  } else if (file_type == "table") {
    ## Use this codepath when we are working with htseq
    count_table <- read.table(files[1], header = header, stringsAsFactors = FALSE)
    colnames(count_table) <- c("rownames", ids[1])
    ## We are going to immediately check for NA, so I think we can suppress warnings.
    count_table[, 2] <- suppressWarnings(as.numeric(count_table[, 2]))
    na_idx <- is.na(count_table[[2]])
    ## This is a bit more circuituous than I would like.
    ## I want to make sure that na_rownames does not evaluate to something
    ## annoying like 'logical(0)' which it will if there are no nas in the count
    ## table and you just ask for is.na(count_table).
    na_rownames <- ""
    if (sum(na_idx) > 0) {
      na_rownames <- count_table[na_idx, "rownames"]
    }
    keepers_idx <- count_table[["rownames"]] != na_rownames
    count_table <- count_table[keepers_idx, ]
    count_table <- as.data.table(count_table)
    count_table <- setkey(count_table, rownames)
    first_rownames <- sort(count_table[["rownames"]])
    if (class(count_table)[1] == "try-error") {
      stop("There was an error reading: ", files[1])
    }
    mesg(files[1], " contains ", length(rownames(count_table)), " rows.")
    ## iterate over and append remaining samples
    for (num in seq(from = 2, to = length(files))) {
      table <- files[num]
      if (file.exists(tolower(table))) {
        table <- tolower(table)
      }
      tmp_count <- try(read.table(table, header = header))
      ## Drop the rows with NAs before coercing to numeric.
      keepers_idx <- tmp_count[[1]] != na_rownames
      tmp_count <- tmp_count[keepers_idx, ]
      current_rownames <- sort(tmp_count[[1]])
      mismatched_rownames <- sum(first_rownames != current_rownames)
      if (mismatched_rownames > 0) {
        warning("The file: ", table, " has mismatched rownames.")
      }
      tmp_count[, 2] <- as.numeric(tmp_count[, 2])
      if (class(tmp_count)[1] == "try-error") {
        stop("There was an error reading: ", table)
      }
      colnames(tmp_count) <- c("rownames", ids[num])
      tmp_count <- as.data.table(tmp_count)
      pre_merge <- nrow(tmp_count)
      if (merge_type == "merge") {
        count_table <- merge(count_table, tmp_count, by = "rownames",
                             all.x = all.x, all.y = all.y)
      } else {
        count_table <- plyr::join(count_table, tmp_count)
      }
      ## rownames(count_table) <- count_table[, "Row.names"]
      ## count_table <- count_table[, -1, drop = FALSE]
      ## post_merge <- length(rownames(count_table))
      post_merge <- nrow(count_table)
      mesg(table, " contains ", pre_merge, " rows and merges to ", post_merge, " rows.")
    } ## End iterating over the files to read.
    ## remove summary fields added by HTSeq
    if (!isTRUE(include_summary_rows)) {
      ## Depending on what happens when the data is read in, these rows may get prefixed with 'X'
      ## In theory, only 1 of these two cases should ever be true.
      htseq_meta_rows <- c("__no_feature", "__ambiguous", "__too_low_aQual",
                           "__not_aligned", "__alignment_not_unique",
                           "X__no_feature", "X__ambiguous", "X__too_low_aQual",
                           "X__not_aligned", "X__alignment_not_unique")
      kept_rows <- !count_table[["rownames"]] %in% htseq_meta_rows
      count_table <- count_table[kept_rows, ]
    } ## End the difference between tximport and reading tables.
    retlist[["count_table"]] <- setkey(count_table, rownames)
    retlist[["source"]] <- "htseq"
  } else if (file_type == "bigwig") {
    count_df <- as.data.frame(rtracklayer::import(files[1]))
    count_df[["rownames"]] <- paste0(count_df[["seqnames"]], "_", count_df[["start"]], "_",
                                     count_df[["end"]])
    count_table <- count_df[, c("rownames", "score")]
    colnames(count_table) <- c("rownames", ids[1])
    tmp_count <- as.data.table(count_table)
    pre_merge <- nrow(count_table)
    ## We are going to immediately check for NA, so I think we can suppress warnings.
    count_table[, 2] <- suppressWarnings(as.numeric(count_table[, 2]))
    na_idx <- is.na(count_table[[2]])
    ## This is a bit more circuituous than I would like.
    ## I want to make sure that na_rownames does not evaluate to something
    ## annoying like 'logical(0)' which it will if there are no nas in the count
    ## table and you just ask for is.na(count_table).
    na_rownames <- ""
    if (sum(na_idx) > 0) {
      na_rownames <- count_table[na_idx, "rownames"]
    }
    keepers_idx <- count_table[["rownames"]] != na_rownames
    count_table <- count_table[keepers_idx, ]
    count_table <- as.data.table(count_table)
    count_table <- setkey(count_table, rownames)
    first_rownames <- sort(count_table[["rownames"]])
    if (class(count_table)[1] == "try-error") {
      stop("There was an error reading: ", files[1])
    }
    mesg(files[1], " contains ", length(rownames(count_table)), " rows.")
    ## iterate over and append remaining samples
    for (num in seq(from = 2, to = length(files))) {
      table <- files[num]
      if (file.exists(tolower(table))) {
        table <- tolower(table)
      }
      tmp_count <- as.data.frame(rtracklayer::import(table))
      tmp_count[["rownames"]] <- paste0(tmp_count[["seqnames"]], "_", tmp_count[["start"]], "_",
                                        tmp_count[["end"]])
      tmp_count <- tmp_count[, c("rownames", "score")]
      colnames(tmp_count) <- c("rownames", ids[num])
      ## Drop the rows with NAs before coercing to numeric.
      keepers_idx <- tmp_count[[1]] != na_rownames
      tmp_count <- tmp_count[keepers_idx, ]
      if (merge_type == "merge") {
        count_table <- merge(count_table, tmp_count, by = "rownames",
                             all.x = all.x, all.y = all.y)
      } else {
        count_table <- plyr::join(count_table, tmp_count)
      }
      ## rownames(count_table) <- count_table[, "Row.names"]
      ## count_table <- count_table[, -1, drop = FALSE]
      ## post_merge <- length(rownames(count_table))
      post_merge <- nrow(count_table)
      mesg(table, " contains ", pre_merge, " rows and merges to ", post_merge, " rows.")
    } ## End iterating over the files to read.
    na_idx <- is.na(count_table)
    count_table[na_idx] <- 0
    retlist[["count_table"]] <- setkey(count_table, rownames)
    retlist[["source"]] <- "bigwig"
  } else {
    stop("I do not understand this count data type.")
  }
  mesg("Finished reading count data.")
  return(retlist)
}
