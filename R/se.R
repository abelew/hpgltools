#' @export
batches <- function(se) {
  batches <- colData(se)[["batch"]]
  names(batches) <- sampleNames(se)
  return(batches)
}
setGeneric("batches")

#' @export
`batches<-` <- function(se, values) {
  colData(se)[["batch"]] <- values
  return(se)
}
setGeneric("batches<-")

#' @export
conditions <- function(se) {
  return(colData(se)[["condition"]])
}

#' @export
`conditions<-` <- function(se, values) {
  colData(se)[["condition"]] <- values
  return(se)
}

#' Create a SummarizedExperiment given some metadata
#'
#' This function was taken from create_expt() and repurposed to create SummarizedExperiments.
#'
#' @param metadata Filename or table of metadata about the samples of interest.
#' @param gene_info Annotations for the genes in the count data.
#' @param count_dataframe Optional table of counts.
#' @param sanitize_rownames Clean up unruly gene IDs?
#' @param sample_colors Specify the colors for the samples?
#' @param title Provide a title for the experiment.
#' @param notes Provide arbitrary notes.
#' @param countdir (deprecated) Directory containing count tables.
#' @param include_type Used to specify types of genes/annotations to use.
#' @param include_gff Keep a copy of the gff with the data?
#' @param file_column Metadata column containing the counts for each sample.
#' @param id_column Non-default column containing the sample IDs.
#' @param savefile Filename to which to save a rda file of the data structure.
#' @param low_files I don't remember this, I bet it is deprecated.
#' @param annotation orgDB associated with this, primarily used with gsva-like tools.
#' @param palette Color palette when auto-choosing colors for the samples.
#' @param round Round the data if/when it is not integer?
#' @param tx_gene_map When using tximport, use this to convert from transcripts to genes.
#' @param ... Extra options.
#' @importFrom SummarizedExperiment SummarizedExperiment metadata<- assays
#' @seealso [summarizedExperiment]
#' @export
create_se <- function(metadata = NULL, gene_info = NULL, count_dataframe = NULL,
                      sanitize_rownames = FALSE, sample_colors = NULL, title = NULL,
                      notes = NULL, include_type = "all", count_source = "htseq",
                      countdir = NULL, include_gff = NULL, file_column = "file", file_type = NULL,
                      savefile = NULL, low_files = FALSE, annotation = "org.Hs.eg.db",
                      palette = "Dark2", round = FALSE, tx_gene_map = NULL,
                      ...) {
  arglist <- list(...)  ## pass stuff like sep=, header=, etc here
  if (is.null(metadata)) {
    stop("This requires some metadata at minimum.")
  }
  ## I am learning about simplifying vs. preserving subsetting
  ## This is a case of simplifying and I believe one which is good because I
  ## just want the string out from my list. Lets assume that palette is in fact
  ## an element in arglist, I really don't care that the name of the resturn is
  ## 'palette'; I already knew that by asking for it.
  if (is.null(title)) {
    title <- "This is a summarized experiment."
  }
  if (is.null(notes)) {
    notes <- glue("Created on {date()}.
")
  }
  ## An expressionset needs to have a Biobase::annotation() in order for
  ## GSEABase to work with it. Reading the documentation, these are primarily
  ## used for naming the type of microarray chip used.
  ## I do not know if any work will need to be done for a SE

  gff_type <- "all"
  if (!is.null(arglist[["include_type"]])) {
    gff_type <- arglist[["include_type"]]
  }

  if (is.null(id_column)) {
    id_column <- "sampleid"
  } else {
    id_column <- tolower(id_column)
    id_column <- gsub(pattern = "[^_[:^punct:]]", replacement = "", x = id_column, perl = TRUE)
  }
  ## Read in the metadata from the provided data frame, csv, or xlsx.
  message("Reading the sample metadata.")
  sample_definitions <- extract_metadata(metadata, id_column = id_column,
                                         ...)
  ## sample_definitions <- extract_metadata(metadata)
  ## Add an explicit removal of the column named 'file' file column if the option file_column is NULL.
  ## This is a just in case measure to avoid conflicts.
  if (is.null(file_column)) {
    if (!is.null(metadata[["file"]])) {
      metadata[["previous_file_column"]] <- metadata[["file"]]
      message("file_column is NULL, moving this column to 'previous_file_column'.")
      metadata[["file"]] <- NULL
    }
  }
  message("The sample definitions comprises: ", nrow(sample_definitions),
          " rows(samples) and ", ncol(sample_definitions),
          " columns(metadata fields).")
  num_samples <- nrow(sample_definitions)
  ## Create a matrix of counts with columns as samples and rows as genes
  ## This may come from either a data frame/matrix, a list of files from the metadata
  ## or it can attempt to figure out the location of the files from the sample names.
  filenames <- NULL
  all_count_tables <- NULL
  ## This set of if() statements is too complex and requires some reworking.
  if (!is.null(count_dataframe)) {
    ## Lets set the order of the count data to that of the sample definitions.
    test_col_rownames <- all.equal(sort(colnames(count_dataframe)),
                                   sort(rownames(sample_definitions)))
    if (isTRUE(test_col_rownames)) {
      count_dataframe <- count_dataframe[, rownames(sample_definitions)]
    } else {
      message("The count table column names are: ",
              toString(sort(colnames(count_dataframe))))
      message("The  meta   data  row  names are: ",
              toString(sort(rownames(sample_definitions))))
      stop("The count table column names are not the same as the sample definition row names.")
    }
    all_count_tables <- data.table::as.data.table(count_dataframe, keep.rownames = "rownames")
    ## If neither of these cases is true, start looking for the files in the
    ## processed_data/ directory
  } else if (is.null(sample_definitions[[file_column]])) {
    message("This will fail because: ", file_column, " does not exist.")
    message("Here are the possible columns: ", toString(colnames(sample_definitions)))
    stop("This requires a column containing the input data.")
  }

  ## At this point sample_definitions$file should be filled in no matter what;
  ## so read the files.
  tximport_data <- NULL
  ## The count_data list should include a set of IDs and tables which are coherent
  ## Therefore, we will want to check in with it later.
  ## Notably, it has slots for: 'kept_ids' which should match 1:1 with the slot 'kept_files',
  ## 'source' which should remind us if the data came from htseq/tximport/etc.
  ## and count_table which should have one column for every kept_id/kept_file.
  count_data <- NULL
  if (is.null(all_count_tables)) {
    ## If all_count_tables does not exist, then we want to read the various files
    ## in the sample definitions to get them.
    filenames <- as.character(sample_definitions[[file_column]])
    sample_ids <- rownames(sample_definitions)
    count_data <- read_counts(sample_ids, filenames, countdir = countdir,
                              file_type = file_type,
                              ...)
    if (count_data[["source"]] == "tximport") {
      tximport_data <- list("raw" = count_data[["tximport"]],
                            "scaled" = count_data[["tximport_scaled"]])
    }
    all_count_tables <- count_data[["count_table"]]
  } else {
    ## if all_count_tables _did_ exist, then we already had the count tables and so
    ## count_data should have them and all ids as 'kept'.
    count_data <- list(
        "source" = "dataframe",
        "raw" = all_count_tables,
        "kept_ids" = rownames(sample_definitions))
    ## Remember that R does not like rownames to start with a number, and if they do
    ## I already changed the count table rownames to begin with 's'.
    count_data[["kept_ids"]] <- gsub(pattern = "^([[:digit:]])",
                                     replacement = "s\\1",
                                     x = count_data[["kept_ids"]])
  }
  ## Here we will prune the metadata for any files/ids which were dropped
  ## when reading in the count tables.
  kept_definitions_idx <- rownames(sample_definitions) %in% count_data[["kept_ids"]]
  if (sum(kept_definitions_idx) < length(kept_definitions_idx)) {
    warning("Some samples were removed when cross referencing the samples against the count data.")
  }
  sample_definitions <- sample_definitions[kept_definitions_idx, ]
  ## While we are removing stuff...
  ## I have had a couple data sets with incomplete counts, get rid of those rows
  ## before moving on.
  all_count_tables <- all_count_tables[complete.cases(all_count_tables), ]

  numeric_columns <- colnames(all_count_tables) != "rownames"
  for (col in colnames(all_count_tables)[numeric_columns]) {
    ## Ensure there are no stupid entries like target_id est_counts
    all_count_tables[[col]] <- as.numeric(all_count_tables[[col]])
  }
  ## Features like exon:alicethegene-1 are annoying and entirely too common in TriTrypDB data
  if (isTRUE(sanitize_rownames)) {
    all_count_tables[["rownames"]] <- gsub(pattern = "^exon:", replacement = "",
                                           x = all_count_tables[["rownames"]])
    all_count_tables[["rownames"]] <- make.names(gsub(pattern = ":\\d+", replacement = "",
                                                      x = all_count_tables[["rownames"]]),
                                                 unique = TRUE)
  }

  ## There is an important caveat here!!
  ## data.table::as.data.table(stuff, keep.rownames='column') will change the
  ## rownames to remove punctuation including ':'!  Which means that if I have a
  ## rowname that looks like 'LmjF.01.0010:mRNA', it will get changed to
  ## 'LmjF.01.0010.mRNA' Which will of course kill any downstream analyses which
  ## depend on consistent rownames between the count table and any tximport
  ## data, since the tximport data  will still have the ':'...
  ## I am not certain what the best solution is, I am thinking perhaps to recast
  ## the tximport data as a set of data tables so that whatever silly stuff it
  ## does it will at least do consistently.  That of course will have unintended
  ## consequences for other tools which use tximport data (DESeq2), but if my
  ## rownames are consistent, then other problems will be easier to handle via
  ## recasting the data to a matrix or df or whatever the downstream tool
  ## requires.
  ## In contrast, I can take a simpler but less transparent route and change the
  ## rownames of all the tximport imported data to match that returned by
  ## as.data.table().
  if (!is.null(tximport_data[["raw"]])) {
    rownames(tximport_data[["raw"]][["abundance"]]) <- gsub(
        pattern = ":", replacement = "\\.",
        x = rownames(tximport_data[["raw"]][["abundance"]]))
    rownames(tximport_data[["raw"]][["counts"]]) <- gsub(
        pattern = ":", replacement = "\\.",
        x = rownames(tximport_data[["raw"]][["counts"]]))
    rownames(tximport_data[["raw"]][["length"]]) <- gsub(
        pattern = ":", replacement = "\\.",
        x = rownames(tximport_data[["raw"]][["length"]]))
  }
  if (!is.null(tximport_data[["scaled"]])) {
    rownames(tximport_data[["scaled"]][["abundance"]]) <- gsub(
        pattern = ":", replacement = "\\.",
        x = rownames(tximport_data[["scaled"]][["abundance"]]))
    rownames(tximport_data[["scaled"]][["counts"]]) <- gsub(
        pattern = ":", replacement = "\\.",
        x = rownames(tximport_data[["scaled"]][["counts"]]))
    rownames(tximport_data[["scaled"]][["length"]]) <- gsub(
        pattern = ":", replacement = "\\.",
        x = rownames(tximport_data[["scaled"]][["length"]]))
  }

  ## Try a couple different ways of getting gene-level annotations into the se.
  annotation <- NULL
  if (is.null(gene_info)) {
    ## Including, if all else fails, just grabbing the gene names from the count tables.
    if (is.null(include_gff)) {
      gene_info <- data.table::as.data.table(all_count_tables[["rownames"]],
                                             keep.rownames = "rownames")
      names(gene_info) <- "rownames"
    } else {
      ## Or reading a gff file.
      message("create_expt(): Reading annotation gff, this is slow.")
      annotation <- load_gff_annotations(gff = include_gff, type = gff_type)
      gene_info <- data.table::as.data.table(annotation, keep.rownames = "rownames")
    }
  } else if (class(gene_info)[[1]] == "list" && !is.null(gene_info[["genes"]])) {
    ## In this case, it is using the output of reading a OrgDB instance
    gene_info <- data.table::as.data.table(gene_info[["genes"]], keep.rownames = "rownames")
  } else if (class(gene_info)[[1]] == "data.table" || class(gene_info)[[1]] == "tbl_df") {
    ## Try to make the data table usage consistent by rownames.
    ## Sometimes we take these from data which did "keep.rownames='some_column'"
    ## Sometimes we take these from data which set rownames(dt)
    ## And sometimes the rownames were never set.
    ## Therefore I will use rownames(dt) as the master, dt$rownames as secondary, and
    ## as a fallback take the first column in the data.
    if (is.null(rownames(gene_info)) && is.null(gene_info[["rownames"]])) {
      gene_info[["rownames"]] <- make.names(rownames[[1]], unique = TRUE)
      message("Both rownames() and $rownames were null.")
    }
  } else {
    gene_info <- data.table::as.data.table(gene_info, keep.rownames = "rownames")
  }

  geneid_check <- grepl(x = all_count_tables[["rownames"]], pattern = "^gene:")
  if (sum(geneid_check) > 0) {
    all_count_tables[["rownames"]] <- gsub(x = all_count_tables[["rownames"]],
                                           pattern = "^gene:", replacement = "")
  }

  ## It turns out that loading the annotation information from orgdb/etc may not set the
  ## row names. Perhaps I should do that there, but I will add a check here, too.
  found_sum <- sum(gene_info[["rownames"]] %in% all_count_tables[["rownames"]])
  if (found_sum == 0) {
    if (!is.null(gene_info[["geneid"]])) {
      gene_info[["rownames"]] <- gene_info[["geneid"]]
      found_sum <- sum(gene_info[["rownames"]] %in% all_count_tables[["rownames"]])
    }
  }
  if (found_sum == 0) {
    warning("Even after changing the rownames in gene info, they do not match the count table.")
    message("Even after changing the rownames in gene info, they do not match the count table.")
    message("Here are the first few rownames from the count tables:")
    message(toString(head(all_count_tables[["rownames"]])))
    message("Here are the first few rownames from the gene information table:")
    message(toString(head(gene_info[["rownames"]])))
  } else {
    message("Matched ", found_sum, " annotations and counts.")
  }

  ## Take a moment to remove columns which are blank
  columns_to_remove <- NULL
  for (col in seq_along(colnames(gene_info))) {
    sum_na <- sum(is.na(gene_info[[col]]))
    sum_null <- sum(is.null(gene_info[[col]]))
    sum_empty <- sum_na + sum_null
    if (sum_empty ==  nrow(gene_info)) {
      ## This column is empty.
      columns_to_remove <- append(columns_to_remove, col)
    }
    ## While we are looping through the columns,
    ## Make certain that no columns in gene_info are lists or factors.
    ## FIXME: 202104: I am no longer sure why I changed factors.
    if (class(gene_info[[col]]) == "factor" ||
        class(gene_info[[col]]) == "AsIs" ||
        class(gene_info[[col]]) == "list") {
      gene_info[[col]] <- as.character(gene_info[[col]])
    }
  }
  if (length(columns_to_remove) > 0) {
    gene_info <- gene_info[-columns_to_remove]
  }

  ## There should no longer be blank columns in the annotation data.
  ## Maybe I will copy/move this to my annotation collection toys?
  ## This temporary id number will be used to ensure that the order of features in everything
  ## will remain consistent, as we will call order() using it later.
  all_count_tables[["temporary_id_number"]] <- seq_len(nrow(all_count_tables))
  message("Bringing together the count matrix and gene information.")
  ## The method here is to create a data.table of the counts and annotation data,
  ## merge them, then split them apart.

  ## Made a small change to check for new tximport rownames in the gene information.
  ## This should automagically check and fix rownames when they would otherwise
  ## not match after using tximport.
  if (!is.null(tx_gene_map)) {
    matched_rows <- sum(rownames(gene_info) %in% tx_gene_map[[2]])
    if (matched_rows < 1) {
      message("The mapped IDs are not the rownames of your gene information, changing them now.")
      if (names(tx_gene_map)[2] %in% colnames(gene_info)) {
        new_name <- names(tx_gene_map)[2]
        rownames(gene_info) <- make.names(tx_gene_map[[new_name]], unique = TRUE)
      } else {
        warning("Cannot find an appropriate column in gene_info, refusing to use the tx_map.")
      }
    }
  }

  counts_and_annotations <- merge(all_count_tables, gene_info, by = "rownames", all.x = TRUE)
  ## In some cases, the above merge will result in columns being set to NA
  ## We should set all the NA fields to something I think.
  na_entries <- is.na(counts_and_annotations)
  if (sum(na_entries) > 0) {
    message("Some annotations were lost in merging, setting them to 'undefined'.")
  }
  counts_and_annotations[na_entries] <- "undefined"
  ## Set an incrementing id number to make absolutely paranoidly certain the
  ## order stays constant.
  counts_and_annotations <- counts_and_annotations[
      order(counts_and_annotations[["temporary_id_number"]]), ]
  ## Pull out the annotation data and convert to data frame.
  kept_columns <- colnames(counts_and_annotations) %in% colnames(gene_info)
  final_annotations <- counts_and_annotations[, kept_columns, with = FALSE]
  final_annotations <- as.data.frame(final_annotations, stringsAsFactors = FALSE)
  rownames(final_annotations) <- final_annotations[["rownames"]]
  final_kept <- colnames(final_annotations) != "rownames"
  final_annotations <- final_annotations[, final_kept]

  final_counts <- counts_and_annotations
  kept_columns <- colnames(counts_and_annotations) %in% colnames(all_count_tables) &
    colnames(counts_and_annotations) != "temporary_id_number"
  final_counts <- final_counts[, kept_columns, with = FALSE]
  final_counts <- as.data.frame(final_counts)
  rownames(final_counts) <- final_counts[["rownames"]]
  final_kept <- colnames(final_counts) != "rownames"
  final_counts <- final_counts[, final_kept]
  final_counts <- as.matrix(final_counts)

  ## I found a non-bug but utterly obnoxious behaivor in R
  ## Imagine a dataframe with 2 entries: TcCLB.511511.3 and TcCLB.511511.30
  ## Then imagine that TcCLB.511511.3 gets removed because it is low abundance.
  ## Then imagine what happens if I go to query 511511.3...
  ## Here are some copy/pasted lines illustrating it:
  ## > find_fiveeleven["TcCLB.511511.3", ]
  ##                 logFC AveExpr    t   P.Value adj.P.Val     B    qvalue
  ## TcCLB.511511.30  5.93   6.315 69.6 1.222e-25 7.911e-25 48.86 9.153e-27
  ## Here is the line in the dataframe documentation explaining this nonsense:
  ## https://stat.ethz.ch/R-manual/R-devel/library/base/html/Extract.data.frame.html
  ## Both [ and [[ extraction methods partially match row names. By default neither partially
  ## match column names, but [[ will if exact = FALSE (and with a warning if exact = NA). If you
  ## want to exact matching on row names use match, as in the examples
  ## How about you go eff yourself?  If you then look carefully at the match help, you will see
  ## that this is a feature, not a bug and that if you want truly exact matches, then the string
  ## must not end with a numeric value... oooo....kkk....
  ## Therefore, the following line replaces a terminal numeric rowname with the number and .
  ## > test_df = data.frame(a = c(1,1,1), b = c(2,2,2))
  ## > rownames(test_df) = c("TcCLB.511511.3","TcCLB.511511.30","bob")
  ## > test_df
  ##                 a b
  ## TcCLB.511511.3  1 2
  ## TcCLB.511511.30 1 2
  ## bob             1 2
  ## > rownames(test_df) <- gsub(pattern = "(\\d)$", replacement = "\\1\\.", x = rownames(test_df))
  ## > test_df
  ##                  a b
  ## TcCLB.511511.3.  1 2
  ## TcCLB.511511.30. 1 2
  ## bob              1 2
  ## This is so stupid I think I am going to go and cry.
  ##rownames(final_annotations) <- gsub(pattern = "(\\d)$", replacement = "\\1\\.",
  ##                                    x = rownames(final_annotations), perl = TRUE)
  ##rownames(final_counts) <- gsub(pattern = "(\\d)$", replacement = "\\1\\.",
  ##                                    x = rownames(final_counts), perl = TRUE)

  ##final_counts <- final_counts[, -1, drop = FALSE]

  ## If the user requests input of non-int counts, fix that here.
  if (isTRUE(round)) {
    final_counts <- round(final_counts)
    less_than <- final_counts < 0
    final_counts[less_than] <- 0
  }

  ## I moved the color choices to this area pretty late in the process to make sure that
  ## there was time to remove unused samples.
  ## Make sure we have a viable set of colors for plots
  chosen_colors <- generate_expt_colors(sample_definitions, sample_colors = sample_colors,
                                        chosen_palette = palette)

  requireNamespace("SummarizedExperiment")
  ## SummarizedExperiments vs. ExpressionSets:
  ## assays() vs. exprs()
  ## rowData()/rowRanges() vs. fData()
  ## colData() vs. pData()
  ## Samples metadata access via $ accessor
  ## Experimental metadata (e.g. publication, lab, sra, whatever) via metadata()
  ## I may need to do some reorganizing to avoid confusion between my single experimental metadata
  ## and the metadata() provided by SummarizedExperiment
  ## Note that metadata() is just a list, so anything may be dumped here.
  se <- SummarizedExperiment(assays = final_counts,
                             rowData = final_annotations,
                             colData = sample_definitions)
  metadata(se)[["notes"]] <- notes
  metadata(se)[["title"]] <- title
  metadata(se)[["annotation"]] <- annotation
  metadata(se)[["gff_file"]] <- include_gff
  ## the 'state' slot in the expt is used to keep track of how the data is modified over time.
  starting_state <- list(
      "filter" = "raw",
      "normalization" = "raw",
      "conversion" = "raw",
      "batch" = "raw",
      "transform" = "raw")
  metadata(se)[["state"]] <- starting_state
  se_conditions <- sample_definitions[["condition"]]
  names(se_conditions) <- rownames(sample_definitions)
  se_batches <- sample_definitions[["batch"]]
  names(se_batches) <- rownames(sample_definitions)
  se_libsizes <- colSums(final_counts)
  names(se_libsizes) <- rownames(sample_definitions)
  metadata(se)[["libsize"]] <- se_libsizes

  if (sum(se_libsizes == 0) > 0) {
    zero_idx <- se_libsizes == 0
    zero_samples <- names(se_libsizes)[zero_idx]
    warning("The following samples have no counts! ", zero_samples)
  }

  ## Save the chosen colors
  names(chosen_colors) <- rownames(sample_definitions)
  metadata(se)[["colors"]] <- chosen_colors
  metadata(se)[["tximport"]] <- tximport_data

  ## Save an rdata file of the se.
  if (is.null(savefile)) {
    if ("character" %in% class(metadata)) {
      savefile <- paste0(gsub(x = basename(metadata), pattern = "^(.*)\\..*",
                          replacement = "\\1"), ".rda")
    } else {
      message("Saving the summarized experiment to 'se.rda'.")
      savefile <- "se.rda"
    }
  }
  save_result <- try(save(se, file = savefile), silent = TRUE)
  if (class(save_result) == "try-error") {
    warning("Saving the summarized experiment object failed, perhaps you do not have permissions?")
  }
  message("The final summarized experiment has ", nrow(exprs(se)),
          " rows and ", ncol(colData(se)), " columns.")
  return(se)
}

#' Analagous function to make_pombe_expt()
#'
#' @param annotation Include annotations?
#' @export
make_pombe_se <- function(annotation = TRUE) {
  fission <- new.env()
  tt <- sm(requireNamespace("fission"))
  tt <- sm(try(attachNamespace("fission"), silent = TRUE))
  tt <- data(fission, envir = fission)
  ## some minor shenanigans to get around the oddities of loading from data()
  fission <- fission[["fission"]]
  meta <- as.data.frame(fission@colData)
  meta[["condition"]] <- glue::glue("{meta[['strain']]}.{meta[['minute']]}")
  meta[["batch"]] <- meta[["replicate"]]
  meta[["sample.id"]] <- rownames(meta)
  meta <- meta[, c("sample.id", "id", "strain", "minute",
                   "replicate", "condition", "batch")]

  fission_data <- fission@assays$data[["counts"]]

  annotations <- NULL
  if (isTRUE(annotation)) {
    ## Neat, it works, and even figures out that the default mart is incorrect by itself.
    pombe_annotations <- try(load_biomart_annotations(
      host = "fungi.ensembl.org", trymart = "fungi_mart",
      trydataset = "spombe_eg_gene",
      gene_requests = c("pombase_transcript", "ensembl_gene_id", "ensembl_transcript_id",
        "hgnc_symbol", "description", "gene_biotype"),
      species = "spombe", overwrite = TRUE))
    if ("try-error" %in% class(pombe_annotations)) {
      warning("There was an error downloading the pombe annotations, this will still return.")
    } else {
      pombe_mart <- pombe_annotations[["mart"]]
      annotations <- pombe_annotations[["annotation"]]
      ## As per create_pombe_expt:
      ## I think ensembl changed the IDs to match and the following line is no longer needed.
      ## rownames(annotations) <- make.names(gsub(pattern = "\\.\\d+$",
      ##                                         replacement = "",
      ##                                         x = rownames(annotations)), unique = TRUE)
    }
  }
  pombe_se <- sm(create_se(metadata = meta,
                           count_dataframe = fission_data,
                           gene_info = annotations))
  detach("package:fission")
  return(pombe_se)
}

get_se_colors <- function(se, keep_underscore = TRUE) {
  all_colors <- colors(se)
  condition_fact <- as.character(colData(se)[["condition"]])
  if (isTRUE(keep_underscore)) {
    condition_fact <- gsub(pattern="[^_[:^punct:]]", replacement = "",
                           x = condition_fact, perl = TRUE)
  } else {
    condition_fact <- gsub(x = condition_fact, pattern = "[[:punct:]]", replacement = "")
  }
  names(all_colors) <- condition_fact
  single_idx <- !duplicated(all_colors)
  all_colors <- all_colors[single_idx]
  return(all_colors)
}

set_se_batches <- function(se, fact, ids = NULL, ...) {
  arglist <- list(...)
  original_batches <- colData(se)[["batch"]]
  original_length <- length(original_batches)
  if (length(fact) == 1) {
    ## Assume it is a column in the design
    if (fact %in% colnames(colData(se))) {
      fact <- colData(se)[[fact]]
    } else {
      stop("The provided factor is not in the design matrix.")
    }
  }

  if (length(fact) != original_length) {
    stop("The new factor of batches is not the same length as the original.")
  }
  colData(se)[["batch"]] <- fact
  message("The number of samples by batch are: ")
  print(table(colData(se)[["batch"]]))
  return(se)
}

set_se_colors <- function(se, colors = TRUE,
                            chosen_palette = "Dark2", change_by = "condition") {
  condition_factor <- as.factor(colData(se)[["condition"]])

  ## Since I already have logic for named characters, just convert a list to one...
  if ("list" %in% class(colors)) {
    new_colors <- as.character(colors)
    names(new_colors) <- names(colors)
    colors <- new_colors
  }

  num_conditions <- length(levels(condition_factor))
  design <- colData(se)
  num_samples <- nrow(design)
  sample_ids <- design[["sampleid"]]
  ## chosen_colors <- se[["conditions"]]
  chosen_colors <- condition_factor
  chosen_names <- names(chosen_colors)
  sample_colors <- NULL
  if (is.null(colors) | isTRUE(colors)) {
    sample_colors <- sm(
      grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  } else if (class(colors) == "factor") {
    if (change_by == "condition") {
      mesg("The new colors are a factor, changing according to condition.")
      ## In this case, we have every color accounted for in the set of conditions.
      colors_allocated <- names(colors) %in% levels(colData(se)[["condition"]])
      if (sum(colors_allocated) < length(colors)) {
        missing_colors <- colors[!colors_allocated]
        stop("Colors for the following categories are not being used: ",
             names(missing_colors), ".")
      }
      possible_conditions <- levels(colData(se)[["condition"]])
      conditions_allocated <- possible_conditions %in% names(colors)
      if (sum(conditions_allocated) < length(possible_conditions)) {
        missing_conditions <- possible_conditions[!conditions_allocated]
        missing_samples <- c()
        for (cond in missing_conditions) {
          missing_by_condition <- colData(se)[["condition"]] == cond
          missing_samples_by_cond <- rownames(colData(se))[missing_by_condition]
          missing_samples <- c(missing_samples, missing_samples_by_cond)
        }
        warning("Some conditions do not have a color: ", missing_conditions, ".")
        warning("These samples are: ", missing_samples, ".")
      }
      mapping <- colors
      chosen_colors <- mapping[as.character(chosen_colors)]
      names(chosen_colors) <- chosen_names
    } else if (change_by == "sample") {
      mesg("The new colors are a factor, changing according to sampleID.")
      ## This is changing them by sample id.
      ## In this instance, we are changing specific colors to the provided colors.
      chosen_colors <- se[["colors"]]
      for (snum in seq_along(names(colors))) {
        sampleid <- names(colors)[snum]
        sample_color <- colors[[snum]]
        chosen_colors[[sampleid]] <- sample_color
      }
    }
    chosen_idx <- complete.cases(chosen_colors)
    chosen_colors <- chosen_colors[chosen_idx]
  } else if (class(colors) == "character") {
    if (is.null(names(colors))) {
      names(colors) <- levels(as.factor(se[["conditions"]]))
    }
    if (change_by == "condition") {
      mesg("The new colors are a character, changing according to condition.")
      ## In this case, we have every color accounted for in the set of conditions.
      mapping <- colors
      pd_factor <- as.factor(colData(se)[["condition"]])
      possible_conditions <- levels(pd_factor)
      colors_allocated <- names(colors) %in% possible_conditions
      if (sum(colors_allocated) < length(colors)) {
        missing_colors <- colors[!colors_allocated]
        warning("Colors for the following categories are not being used: ",
                names(missing_colors), ".")
      }
      conditions_allocated <- possible_conditions %in% names(colors)
      if (sum(conditions_allocated) < length(possible_conditions)) {
        missing_conditions <- possible_conditions[!conditions_allocated]
        missing_samples <- c()
        for (cond in missing_conditions) {
          missing_by_condition <- colData(se)[["condition"]] == cond
          missing_samples_by_cond <- rownames(colData(se))[missing_by_condition]
          missing_samples <- c(missing_samples, missing_samples_by_cond)
        }
        warning("Some conditions do not have a color: ", missing_conditions, ".")
        warning("These samples are: ", missing_samples, ".")
      }
      chosen_colors <- mapping[as.character(chosen_colors)]
      names(chosen_colors) <- chosen_names
    } else if (change_by == "sample") {
      mesg("The new colors are a character, changing according to sampleID.")
      ## This is changing them by sample id.
      ## In this instance, we are changing specific colors to the provided colors.
      chosen_colors <- se[["colors"]]
      for (snum in seq_along(names(colors))) {
        sampleid <- names(colors)[snum]
        sample_color <- colors[[snum]]
        chosen_colors[[sampleid]] <- sample_color
      }
    }
    chosen_idx <- complete.cases(chosen_colors)
    chosen_colors <- chosen_colors[chosen_idx]
  } else if (class(colors) == "list") {
    if (change_by == "condition") {
      mesg("The new colors are a list, changing according to condition.")
      ## In this case, we have every color accounted for in the set of conditions.
      mapping <- as.character(colors)
      names(mapping) <- names(colors)
      chosen_colors <- mapping[chosen_colors]
    } else if (change_by == "sample") {
      mesg("The new colors are a list, changing according to sampleID.")
      ## This is changing them by sample id.
      ## In this instance, we are changing specific colors to the provided colors.
      chosen_colors <- se[["colors"]]
      for (snum in seq_along(names(colors))) {
        sampleid <- names(colors)[snum]
        sample_color <- colors[[snum]]
        chosen_colors[[sampleid]] <- sample_color
        ## Set the condition for the changed samples to something unique.
        original_condition <- colData(se)[sampleid, "condition"]
        changed_condition <- glue("{original_condition}{snum}")
        ## se[["design"]][sampleid, "condition"] <- changed_condition
        tmp_pdata <- colData(se)
        old_levels <- levels(tmp_pdata[["condition"]])
        new_levels <- c(old_levels, changed_condition)
        levels(tmp_pdata[["condition"]]) <- new_levels
        tmp_pdata[sampleid, "condition"] <- changed_condition
        colData(se[["expressionset"]]) <- tmp_pdata
      }
    }
    chosen_idx <- complete.cases(chosen_colors)
    chosen_colors <- chosen_colors[chosen_idx]
  } else if (is.null(colors)) {
    mesg("Setting colors according to a color ramp.")
    colors <- sm(
      grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    ## Check that all conditions are named in the color list:
    mapping <- setNames(colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  } else {
    warning("Number of colors provided does not match the number of conditions nor samples.")
    warning("Unsure of what to do, so choosing colors with RColorBrewer.")
    sample_colors <- suppressWarnings(
      grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(num_conditions, chosen_palette))(num_conditions))
    mapping <- setNames(sample_colors, unique(chosen_colors))
    chosen_colors <- mapping[chosen_colors]
  }

  ## Catchall in case I forgot to set the names before now.
  names(chosen_colors) <- chosen_names
  metadata(se)[["colors"]] <- chosen_colors
  return(se)
}

set_se_conditions <- function(se, fact = NULL, ids = NULL, prefix = NULL,
                              null_cell = "null", colors = TRUE,
                              ...) {
  arglist <- list(...)
  if (!is.null(arglist[["factor"]])) {
    warning("I probably should change this argument to factor, but it is 'fact'.")
    fact <- arglist[["factor"]]
  }
  original_conditions <- colData(se)[["condition"]]
  original_length <- length(original_conditions)
  original_num_conditions <- length(levels(as.factor(original_conditions)))
  new_se <- se  ## Explicitly copying se to new_se
  ## because when I run this as a function call() it seems to be not properly setting
  ## the conditions and I do not know why.
  fact_vector <- NULL
  fact_name <- "condition"
  if (!is.null(ids)) {
    ## Change specific id(s) to given condition(s).
    mesg("Setting condition for ids ", toString(ids), " to ", fact, ".")
    old_pdata <- colData(se)
    old_cond <- as.character(old_pdata[["condition"]])
    names(old_cond) <- rownames(old_pdata)
    new_cond <- old_cond
    new_cond[ids] <- fact
    new_pdata <- old_pdata
    new_pdata[["condition"]] <- as.factor(new_cond)
    colData(new_se) <- new_pdata
  } else if (length(fact) == 1) {
    fact_name <- fact
    ## Assume it is a column in the design
    if (fact %in% colnames(colData(se))) {
      new_fact <- colData(se)[[fact]]
      null_ids <- is.na(new_fact) | is.null(new_fact)
      ## Only do this if there are some null entries.
      if (sum(null_ids) > 0) {
        new_fact[null_ids] <- null_cell
      }
      if (!is.null(prefix)) {
        new_fact <- paste0(prefix, new_fact)
      }
      fact_vector <- new_fact
      colData(new_se)[["condition"]] <- new_fact
      ## new_se[["design"]][["condition"]] <- new_fact
    } else {
      stop("The provided factor is not in the design matrix.")
    }
  } else if (length(fact) != original_length) {
    stop("The new factor of conditions is not the same length as the original.")
  } else {
    colData(new_se)[["condition"]] <- fact
    ## new_se[["design"]][["condition"]] <- fact
    fact_vector <- fact
  }

  message("The numbers of samples by condition are: ")
  print(table(colData(new_se)[["condition"]]))
  condition_states <- levels(as.factor(colData(new_se)[["condition"]]))
  if (class(colors)[1] == "list") {
    ## A list of colors may either be a color_choices list or
    ## a hash of states->color which could/should be a named vector.
    color_state_names <- names(colors)
    found_colors <- sum(color_state_names %in% condition_states)
    found_names <- sum(fact_name %in% color_state_names)
    ## In this first instance, the choices should be in this element.
    if (found_names > 0) {
      mesg("The colors appear to be a list delineated by state name.")
      colors <- colors[[fact]]
    } else if (found_colors > 0) {
      mesg("The colors appear to be a single list delineated by condition.")
    } else {
      message("A list of colors was provided, but element ", fact,
              " is not in it; using defaults")
      colors <- NULL
    }
  }
  new_se <- set_se_colors(new_se, colors = colors)
  return(new_se)
}

subset_se <- function(se, subset = NULL, ids = NULL,
                      nonzero = NULL, coverage = NULL) {
  starting_se <- se
  starting_metadata <- colData(se)
  starting_samples <- sampleNames(se)

  if (!is.null(ids)) {
    idx <- starting_samples %in% ids
    se <- se[, idx]
  }

  note_appended <- NULL
  subset_design <- NULL
  if (is.null(coverage) && is.null(nonzero)) {
    if (is.null(subset)) {
      subset_design <- starting_metadata
    } else {
      mesg("Using a subset expression, before subsetting there are: ", nrow(starting_metadata),
           " samples.")
      r_expression <- glue("subset(starting_metadata, {subset})")
      subset_design <- eval(parse(text = r_expression))
      subset_idx <- rownames(starting_metadata) %in% rownames(subset_design)
      note_appended <- glue("Subsetted with {subset} on {date()}.
")
      mesg("Following subsetting, there are: ", nrow(subset_design), " samples.")
      se <- se[, subset_idx]
    }
    if (nrow(subset_design) == 0) {
      stop("When the subset was taken, the resulting design has 0 members.")
    }
    subset_design <- as.data.frame(subset_design, stringsAsFactors = FALSE)
  } else if (is.null(nonzero)) {
    ## If coverage is defined, then use it to subset based on the minimal desired coverage
    ## Perhaps in a minute I will make this work for strings like '1z' to get the lowest
    ## standard deviation or somesuch...
    mesg("Subsetting given a minimal number of counts/sample.")
    coverages <- colSums(assay(se))

    if (is.null(colData(se)[["sample_coverage"]])) {
      colData(se)[["sample_coverage"]] <- coverages
    }
    subset_idx <- coverages >= as.numeric(coverage) ## In case I quote it on accident.
    lost_samples <- sampleNames(se)[!subset_idx]
    se <- se[, subset_idx]
    message("The samples removed (and read coverage) when filtering samples with less than ",
            coverage, " reads are: ")
    print(lost_samples)
    print(colData(se)[["sample_coverage"]])
  } else if (is.null(coverage)) {
    ## Remove samples with less than this number of non-zero genes.
    nonzero_idx <- assay(se) != 0
    num_nonzero <- colSums(nonzero_idx)
    if (is.null(pData(se)[["num_nonzero"]])) {
      colData(se)[["num_nonzero"]] <- num_nonzero
    }
    remove_idx <- num_nonzero < nonzero
    if (sum(remove_idx) == 0) {
      message("No samples have fewer than ", nonzero, " observed genes.")
      return(se)
    }
    samples_dropped <- num_nonzero[remove_idx]
    message("The samples (and read coverage) removed when filtering ",
            nonzero, " non-zero genes are: ")
    print(colSums(exprs(se))[remove_idx])
    print(num_nonzero[remove_idx])
    se <- se[, !remove_idx]
  } else {
    stop("Unable to determine what is being subset.")
  }

  ## This is to get around stupidity with respect to needing all factors to be
  ## in a DESeqDataSet
  notes <- se[["notes"]]
  if (!is.null(note_appended)) {
    notes <- glue("{notes}{note_appended}")
  }

  return(se)
}

## EOF
