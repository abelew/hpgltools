#' Create a SummarizedExperiment given some metadata
#'
#' This function was taken from create_expt() and repurposed to create
#' SummarizedExperiments.
#'
#' @param metadata Filename or table of metadata about the samples of interest.
#' @param gene_info Annotations for the genes in the count data.
#' @param count_dataframe Optional table of counts.
#' @param sanitize_rownames Clean up unruly gene IDs?
#' @param sample_colors Specify the colors for the samples?
#' @param title Provide a title for the experiment.
#' @param notes Provide arbitrary notes.
#' @param include_type Used to specify types of genes/annotations to use.
#' @param count_source Explicitly specify the method used to create the text table.
#' @param countdir (deprecated) Directory containing count tables.
#' @param include_gff Keep a copy of the gff with the data?
#' @param file_column Metadata column containing the counts for each sample.
#' @param file_type Force a specific source of count data instead
#'  of autodetecting it.
#' @param id_column Non-default column containing the sample IDs.
#' @param handle_na How to handle NA values in the data.
#' @param researcher When creating gene sets, set the researcher.
#' @param study_name When creating gene sets, set the study name.
#' @param feature_type When matching annotations, use this feature type.
#' @param ignore_tx_version tximport can strictly match transcript/gene
#'  versions, or not.
#' @param savefile Filename to which to save a rda file of the data structure.
#' @param low_files I don't remember this, I bet it is deprecated.
#' @param annotation orgDB associated with this, primarily used with gsva-like tools.
#' @param palette Color palette when auto-choosing colors for the samples.
#' @param round Round the data if/when it is not integer?
#' @param tx_gene_map When using tximport, use this to convert from
#'  transcripts to genes.
#' @param ... Extra options.
#' @importFrom SummarizedExperiment SummarizedExperiment metadata<- assays
#' @importFrom S4Vectors metadata
#' @seealso [summarizedExperiment]
#' @export
create_se <- function(metadata = NULL, gene_info = NULL, count_dataframe = NULL,
                      sanitize_rownames = FALSE, sample_colors = NULL, title = NULL,
                      notes = NULL, include_type = "all", count_source = "htseq",
                      countdir = NULL, include_fasta = NULL,  include_gff = NULL,
                      file_column = "file", file_type = NULL, id_column = NULL,
                      handle_na = "drop", researcher = "elsayed", study_name = NULL,
                      feature_type = "gene", ignore_tx_version = TRUE, savefile = NULL,
                      low_files = FALSE, annotation = NULL, palette = "Dark2",
                      round = FALSE, tx_gene_map = NULL, species = NULL,
                      condition_column = NULL, batch_column = NULL,
                      ...) {
  arglist <- list(...)  ## pass stuff like sep=, header=, etc here
  if (!is.null(species)) {
    include_gff <- file.path(Sys.getenv("HOME"), "libraries", "genome", "gff", paste0(species, ".gff"))
    gene_info <- include_gff
    include_fasta <- file.path(Sys.getenv("HOME"), "libraries", "genome", "fasta", paste0(species, ".fasta"))
  }
  if (is.null(metadata)) {
    stop("This requires some metadata at minimum.")
  }

  if ("gene_tx_map" %in% names(arglist)) {
    message("Hey, it is tx_gene_map, not gene_tx_map!")
    tx_gene_map <- arglist[["gene_tx_map"]]
  }
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
                                         condition_column = condition_column,
                                         batch_column = batch_column,
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
                              file_type = file_type, tx_gene_map = tx_gene_map,
                              ignore_tx_version = ignore_tx_version,
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

  test_df <- as.data.frame(all_count_tables)
  rownames(test_df) <- test_df[["rownames"]]
  test_df[["rownames"]] <- NULL
  count_nas <- is.na(test_df)
  if (sum(count_nas) > 0) {
    warning("There are some NAs in this data, the 'handle_nas' parameter may be required.")
  }
  check_counts <- colSums(test_df, na.rm = TRUE)
  zero_count_samples <- check_counts == 0
  if (sum(zero_count_samples) > 0) {
    warning("The following samples have no counts: ", toString(names(check_counts)[zero_count_samples]))
    message("If the handle_na parameter is 'drop', this will result in an empty dataset.")
  }

  if (handle_na == "drop") {
    all_count_tables <- all_count_tables[complete.cases(all_count_tables), ]
  } else {
    na_idx <- is.na(all_count_tables)
    all_count_tables[na_idx] <- 0
  }

  numeric_columns <- colnames(all_count_tables) != "rownames"
  for (col in colnames(all_count_tables)[numeric_columns]) {
    ## Ensure there are no stupid entries like target_id est_counts
    all_count_tables[[col]] <- as.numeric(all_count_tables[[col]])
  }

  ## While we are removing stuff...
  ## I have had a couple data sets with incomplete counts, get rid of those rows
  ## before moving on.
  all_count_tables <- all_count_tables[complete.cases(all_count_tables), ]

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
  ## FIXME: Move this entirely to a gr and use dispatch.
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
  } else if (class(gene_info)[[1]] == "character") { ## Assume the input is a gff
    gene_info <- gff2gr(gene_info)
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

  counts_and_annotations <- merge_counts_annotations(gene_info, all_count_tables, tx_gene_map)
  final_annotations <- counts_and_annotations[["final_annotations"]]
  final_counts <- counts_and_annotations[["final_counts"]]
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
  if (is.null(study_name)) {
    study_name <- basename(getwd())
  }
  metadata(se)[["study"]] <- study_name
  metadata(se)[["researcher"]] <- researcher
  metadata(se)[["title"]] <- title
  genome_data <- NULL
  if (!is.null(include_fasta)) {
    genome_data <- Biostrings::readDNAStringSet(include_fasta)
    metadata(se)[["genome"]] <- genome_data
  }
  grange_data <- NULL
  if (!is.null(include_gff)) {
    grange_data <- gff2gr(include_gff)
    if (!is.null(genome_data)) {
      seqinfo(grange_data) <- seqinfo(genome_data)
    }
  }
  metadata(se)[["grange"]] <- grange_data
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
#' @param host ensembl host to query
#' @export
make_pombe_se <- function(annotation = TRUE, host = "nov2020-fungi.ensembl.org") {
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
      host = host, trymart = "fungi_mart",
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

#' Use an expression to subset a summarized experiment.
#'
#' I like just passing an expression string to get subsets.
#'
#' @param se Input se.
#' @param subset expression to use to subset on the metadata.
#' @param ids Optional vector of sample IDs.
#' @param nonzero A number of nonzero genes to use instead.
#' @param coverage A minimum coverage to use instead.
#' @param print_excluded Print the sampleIDs excluded by this subset.
#' @export
subset_se <- function(se, subset = NULL, ids = NULL,
                      nonzero = NULL, coverage = NULL,
                      print_excluded = TRUE) {
  starting_se <- se
  starting_metadata <- colData(se)
  starting_samples <- sampleNames(se)
  starting_colors <- get_colors(se)
  end_colors <- starting_colors
  current_libsize <- libsize(se)
  subset_libsize <- current_libsize
  if (!is.null(ids)) {
    idx <- starting_samples %in% ids
    se <- se[, idx]
    subset_libsize <- subset_libsize[idx]
    end_colors <- starting_colors[idx]
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
      subset_libsize <- subset_libsize[subset_idx]
      end_colors <- starting_colors[subset_idx]
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
    subset_libsize <- subset_libsize[subset_idx]
    end_colors <- starting_colors[subset_idx]
    if (isTRUE(print_excluded)) {
    message("The samples removed (and read coverage) when filtering samples with less than ",
            coverage, " reads are: ")
    message(toString(lost_samples))
    print(colData(se)[["sample_coverage"]])
    }
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
    if (isTRUE(print_excluded)) {
      message("Samples removed: ", toString(samples_dropped))
    }
    se <- se[, !remove_idx]
    end_colors <- starting_colors[!remove_idx]
    subset_libsize <- subset_libsize[!remove_idx]
  } else {
    stop("Unable to determine what is being subset.")
  }
  ## This is to get around stupidity with respect to needing all factors to be
  ## in a DESeqDataSet
  notes <- se[["notes"]]
  if (!is.null(note_appended)) {
    notes <- glue("{notes}{note_appended}")
  }
  colors(se) <- end_colors
  ## If we have condition/batch factors, droplevels them.
  if (!is.null(colData(se)[["condition"]])) {
    colData(se)[["condition"]] <- droplevels(as.factor(colData(se)[["condition"]]))
  }
  if (!is.null(colData(se)[["batch"]])) {
    colData(se)[["batch"]] <- droplevels(as.factor(colData(se)[["batch"]]))
  }

  S4Vectors::metadata(se)[["libsize"]] <- subset_libsize
  return(se)
}

#' Merge gene annotations and count tables when the gene_info is a data table/dataframe
merge_counts_annotations <- function(gene_info, all_count_tables, tx_gene_map = NULL) {
  ## It turns out that loading the annotation information from orgdb/etc may not set the
  ## row names. Perhaps I should do that there, but I will add a check here, too.
  all_count_tables[["temporary_id_number"]] <- seq_len(nrow(all_count_tables))
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

  ## This should automagically check and fix rownames when they would otherwise
  ## not match after using tximport.
  if (!is.null(tx_gene_map)) {
    matched_rows <- sum(gene_info[["rownames"]] %in% tx_gene_map[[2]])
    if (matched_rows < 1) {
      message("The mapped IDs are not the rownames of your gene information, changing them now.")
      if (names(tx_gene_map)[2] %in% colnames(gene_info)) {
        new_name <- names(tx_gene_map)[2]
        gene_info[["rownames"]] <- make.names(tx_gene_map[[new_name]], unique = TRUE)
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
  final_counts[["rownames"]] <- NULL
  final_counts <- as.matrix(final_counts)
  retlist <- list(
    "final_annotations" = final_annotations,
    "final_counts" = final_counts)
  return(retlist)
}
setGeneric("merge_counts_annotations")

#' @export
setMethod(
  "merge_counts_annotations", signature = signature(gene_info = "GRanges"),
  definition = function(gene_info, all_count_tables, tx_gene_map = NULL) {
    ## First find the best match between the GR mcols and count table rownames
    assay_rownames <- all_count_tables[["rownames"]]
    all_count_tables[["temporary_id_number"]] <- seq_len(nrow(all_count_tables))
    max_matches <- 0
    gene_info_df <- as.data.frame(mcols(gene_info))
    potentials <- colnames(gene_info_df)
    chosen_column <- NULL
    for (r in potentials) {
      hits <- sum(gene_info_df[[r]] %in% assay_rownames)
      if (hits > max_matches) {
        chosen_column <- r
        max_matches <- hits
      }
    }
    if (is.null(chosen_column)) {
      stop("I did not find a matching column.")
    } else {
      gene_info_df[["rownames"]] <- gene_info_df[[chosen_column]]
    }
    ## Get rid of duplicates and NAs
    na_rownames <- is.na(gene_info_df[["rownames"]])
    gene_info_df <- gene_info_df[!na_rownames, ]
    dup_rownames <- duplicated(gene_info_df[["rownames"]])
    gene_info_df <- gene_info_df[!dup_rownames, ]

    ## Get rid of columns that were not filled in for the chosen rowtype
    remaining_rows <- nrow(gene_info_df)
    for (r in potentials) {
      num_nas <- sum(is.na(gene_info_df[[r]]))
      if (num_nas == remaining_rows) {
        mesg("Dropping ", r, " it is comprised of only NA.")
        gene_info_df[[r]] <- NULL
      }
    }

    ## This should automagically check and fix rownames when they would otherwise
    ## not match after using tximport.
    if (!is.null(tx_gene_map)) {
      matched_rows <- sum(gene_info_df[[chosen_column]] %in% tx_gene_map[[2]])
      if (matched_rows < 1) {
        message("The mapped IDs are not the rownames of your gene information, changing them now.")
        if (names(tx_gene_map)[2] %in% colnames(gene_info_df)) {
          new_name <- names(tx_gene_map)[2]
          gene_info_df[["rownames"]] <- make.names(tx_gene_map[[new_name]], unique = TRUE)
        } else {
          warning("Cannot find an appropriate column in gene_info, refusing to use the tx_map.")
        }
      }
    }
    counts_and_annotations <- merge(all_count_tables, gene_info_df, by = "rownames", all.x = TRUE)
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
    counts_and_annotations[["temporary_id_number"]] <- NULL
    kept_columns <- colnames(counts_and_annotations) %in% colnames(gene_info_df)
    final_annotations <- counts_and_annotations[, kept_columns, with = FALSE]
    final_annotations <- as.data.frame(final_annotations, stringsAsFactors = FALSE)
    rownames(final_annotations) <- final_annotations[["rownames"]]
    final_annotations[["rownames"]] <- NULL
    final_counts <- counts_and_annotations
    kept_columns <- colnames(counts_and_annotations) %in% colnames(all_count_tables)
    final_counts <- final_counts[, kept_columns, with = FALSE]
    final_counts <- as.data.frame(final_counts, stringsAsFactors = FALSE)
    rownames(final_counts) <- final_counts[["rownames"]]
    final_counts[["rownames"]] <- NULL
    final_counts <- as.matrix(final_counts)
    retlist <- list(
      "final_annotations" = final_annotations,
      "final_counts" = final_counts)
    return(retlist)
  })

#' Make pretty xlsx files of count data.
#'
#' Some folks love excel for looking at this data.  ok.
#'
#' Tested in test_03graph_metrics.R
#' This performs the following:  Writes the raw data, graphs the raw data,
#' normalizes the data, writes it, graphs it, and does a median-by-condition and
#' prints that.  I replaced the openxlsx function which writes images into xlsx
#' files with one which does not require an opening of a pre-existing plotter.
#' Instead it (optionally)opens a pdf device, prints the plot to it, opens a png
#' device, prints to that, and inserts the resulting png file.  Thus it
#' sacrifices some flexibility for a hopefully more consistent behaivor.  In
#' addition, one may use the pdfs as a set of images importable into illustrator
#' or whatever.
#'
#' @param se An expressionset to print.
#' @param excel Filename to write.
#' @param norm Normalization to perform.
#' @param violin Include violin plots?
#' @param sample_heat Include sample heatmaps?
#' @param convert Conversion to perform.
#' @param transform Transformation used.
#' @param batch Batch correction applied.
#' @param filter Filtering method used.
#' @param med_or_mean When printing mean by condition, one may want median.
#' @param color_na Color cells which were NA before imputation this color.
#' @param merge_order Used to decide whether to put the counts or annotations first when
#'  printing count tables.
#' @param ... Parameters passed down to methods called here (graph_metrics, etc).
#' @return A big honking excel file and a list including the dataframes and images created.
#' @seealso [openxlsx] [Biobase] [normalize()] [graph_metrics()]
#' @example inst/examples/se.R
#' @export
write_se <- function(se, excel = "excel/pretty_counts.xlsx", norm = "quant",
                       violin = TRUE, sample_heat = NULL, convert = "cpm", transform = "log2",
                       batch = "svaseq", filter = TRUE, med_or_mean = "mean",
                       color_na = "#DD0000", merge_order = "counts_first", ...) {
  arglist <- list(...)
  xlsx <- init_xlsx(excel)
  wb <- xlsx[["wb"]]
  excel_basename <- xlsx[["basename"]]
  new_row <- 1
  new_col <- 1

  ## Set up a vector of images to clean up when finished.
  image_files <- c()

  ## I think the plot dimensions should increase as the number of samples increase
  plot_dim <- 6
  num_samples <- ncol(exprs(se))
  if (num_samples > 12) {
    plot_dim <- ceiling(num_samples / 4)
  }
  plot_cols <- floor(plot_dim * 1.5)
  plot_rows <- ceiling(plot_dim * 5.0)

  ## Write an introduction to this foolishness.
  message("Writing the first sheet, containing a legend and some summary data.")
  sheet <- "legend"
  norm_state <- glue("{transform}({convert}({norm}({batch}({filter}(counts)))))")
  legend <- data.frame(
    "sheet" = c("1.", "2.", "3.", "4.", "5.", "6."),
    "sheet_definition" = c(
      "This sheet, including the experimental design.",
      "The raw counts and annotation data on worksheet 'raw_data'.",
      "Some graphs describing the distribution of raw data in worksheet 'raw_plots'.",
      glue("The counts normalized with: {norm_state}"),
      "Some graphs describing the distribution of the normalized data on 'norm_plots'.",
      "The median normalized counts by condition factor on 'median_data'."),
    stringsAsFactors = FALSE)
  colnames(legend) <- c("Worksheets", "Contents")
  xls_result <- write_xlsx(data = legend, wb = wb, sheet = sheet, rownames = FALSE,
                           title = "Columns used in the following tables.")
  rows_down <- nrow(legend)
  new_row <- new_row + rows_down + 3
  annot <- as.data.frame(pData(se), strinsAsFactors = FALSE)
  xls_result <- write_xlsx(data = annot, wb = wb, start_row = new_row, rownames = FALSE,
                           sheet = sheet, start_col = 1, title = "Experimental Design.")

  ## Get the library sizes and other raw plots before moving on...
  do_qq <- FALSE
  do_sample_heat <- FALSE
  if (ncol(exprs(se)) < 18) {
    do_qq <- TRUE
    do_sample_heat <- TRUE
  }

  metrics <- graph_metrics(se, qq = do_qq, gene_heat = do_sample_heat,
                           ...)
  new_row <- new_row + nrow(pData(se)) + 3
  libsizes <- as.data.frame(metrics[["libsizes"]])[, c("id", "sum", "condition")]
  xls_result <- write_xlsx(data = libsizes, wb = wb, start_row = new_row,
                           rownames = FALSE, sheet = sheet, start_col = 1,
                           title = "Library sizes.")

  new_row <- new_row + nrow(libsizes) + 3
  libsize_summary <- as.data.frame(metrics[["libsize_summary"]])
  xls_result <- write_xlsx(data = libsize_summary, wb = wb, start_row = new_row,
                           rownames = FALSE, sheet = sheet, start_col = 1,
                           title = "Library size summary.")

  ## Write the raw read data and gene annotations
  mesg("Writing the raw reads.")
  sheet <- "raw_reads"
  new_row <- 1
  new_col <- 1
  reads <- exprs(se)
  info <- fData(se)

  if (!is.null(info[["Row.names"]])) {
    ridx <- colnames(info) == "Row.names"
    message("Hey, you merged the annotation data and did not reset the column names!")
    colnames(info)[ridx] <- "old_row_names"
  }
  read_info <- merge(info, reads, by = "row.names")
  xls_result <- write_xlsx(data = read_info, wb = wb, sheet = sheet, rownames = FALSE,
                           start_row = new_row, start_col = new_col, title = "Raw Reads.")

  ## Potentially useful for proteomics data and subtracted data.
  if (!is.null(color_na)) {
    na_style <- openxlsx::createStyle(fontColour = color_na)
    nas <- se[["na_values"]]
    for (col in colnames(nas)) {
      row_definition <- which(nas[[col]]) + 2
      col_idx <- colnames(read_info) == col
      col_definition <- which(col_idx)
      colored <- openxlsx::addStyle(wb, sheet = sheet, na_style,
                                    rows = row_definition, cols = col_definition)
    }
  }

  ## Write some graphs for the raw data
  mesg("Graphing the raw reads.")
  sheet <- "raw_graphs"
  newsheet <- try(openxlsx::addWorksheet(wb, sheetName = sheet))
  if (class(newsheet) == "try-error") {
    warning("Failed to add the sheet: ", sheet)
  }

  ## Start with library sizes.
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Legend.",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Raw library sizes.",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Non-zero genes.",
                                      startRow = new_row, startCol = new_col)
  new_row <- new_row + 1
  new_col <- 1
  legend_plot <- metrics[["legend"]]
  try_result <- xlsx_insert_png(legend_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "01_legend", savedir = excel_basename, fancy_type = "svg")
  if ("try-error" %in% class(try_result)) {
    warning("Failed to add the raw legend.")
  } else {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  libsize_plot <- metrics[["libsize"]]
  try_result <- xlsx_insert_png(libsize_plot, wb = wb, sheet = sheet,
                                width = plot_dim, height = plot_dim,
                                start_col = new_col, start_row = new_row,
                                plotname = "02_libsize", savedir = excel_basename)
  if ("try-error" %in% class(try_result)) {
    warning("Failed to add the raw library sizes.")
  } else {
    image_files <- c(image_files, try_result[["filename"]])
  }

  ## Same row, non-zero plot
  new_col <- new_col + plot_cols + 1
  nonzero_plot <- metrics[["nonzero"]]
  try_result <- xlsx_insert_png(nonzero_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "03_nonzero", savedir = excel_basename)
  if ("try-error" %in% class(try_result)) {
    warning("Failed to add the raw nonzero plot.")
  } else {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1

  ## Visualize distributions
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Raw data density plot.",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Raw Boxplot.",
                                      startRow = new_row, startCol = new_col)
  new_col <- 1
  density_plot <- metrics[["density"]]
  new_row <- new_row + 1
  try_result <- xlsx_insert_png(density_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "04_density", savedir = excel_basename)
  if ("try-error" %in% class(try_result)) {
    warning("Failed to add the raw density plot.")
  } else {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  boxplot_plot <- metrics[["boxplot"]]
  try_result <- xlsx_insert_png(boxplot_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "05_boxplot", savedir = excel_basename)
  if ("try-error" %in% class(try_result)) {
    warning("Failed to add the raw boxplot.")
  } else {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  topn_plot <- metrics[["topnplot"]]
  try_result <- xlsx_insert_png(topn_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "06_topnplot", savedir = excel_basename)
  if ("try-error" %in% class(try_result)) {
    warning("Failed to add the raw top-n plot.")
  } else {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  cv_plot <- metrics[["cvplot"]]
  try_result <- xlsx_insert_png(cv_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "07_cvplot", savedir = excel_basename)
  if ("try-error" %in% class(try_result)) {
    warning("Failed to insert the raw cv plot.")
  } else {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- 1

  ## Move down next set of rows, heatmaps
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Raw correlation heatmap.",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Raw distance heatmap.",
                                      startRow = new_row, startCol = new_col)
  new_col <- 1
  new_row <- new_row + 1
  corheat_plot <- metrics[["corheat"]]
  try_result <- xlsx_insert_png(corheat_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "08_corheat", savedir = excel_basename)
  if ("try-error" %in% class(try_result)) {
    warning("Failed to add the raw correlation heatmap.")
  } else {
    image_files <- c(image_files, try_result[["filename"]])
  }
  disheat_plot <- metrics[["disheat"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(disheat_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "09_disheat", savedir = excel_basename)
  if ("try-error" %in% class(try_result)) {
    warning("Failed to add the raw distance heatmap.")
  } else {
    image_files <- c(image_files, try_result[["filename"]])
  }
  if (isTRUE(sample_heat)) {
    tmp_se <- sm(normalize(se, transform = "log2", filter = TRUE))
    sampleheat_plot <- plot_sample_heatmap(tmp_se)
    new_col <- new_col + plot_cols + 1
    try_result <- xlsx_insert_png(sampleheat_plot[["plot"]], wb = wb,
                                  sheet = sheet, width = plot_dim, height = plot_dim,
                                  start_col = new_col, start_row = new_row,
                                  plotname = "09a_sampleheat", savedir = excel_basename)
    if ("try-error" %in% class(try_result)) {
      warning("Failed to add the sample heatmap.")
    } else {
      image_files <- c(image_files, try_result[["filename"]])
    }
  }
  new_col <- 1

  ## SM plots
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Raw standard median correlation.",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Raw standard distance correlation.",
                      startRow = new_row, startCol = new_col)
  new_col <- 1
  new_row <- new_row + 1
  smc_plot <- metrics[["smc"]]
  try_result <- xlsx_insert_png(smc_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "10_smc", savedir = excel_basename, fancy_type = "svg")
  if ("try-error" %in% class(try_result)) {
    warning("Failed to add the raw correlation standard median plot.")
  } else {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  smd_plot <- metrics[["smd"]]
  try_result <- xlsx_insert_png(smd_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "11_smd", savedir = excel_basename, fancy_type = "svg")
  if ("try-error" %in% class(try_result)) {
    warning("Failed to add the raw distance standard median plot.")
  } else {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- 1

  ## PCA, PCA(l2cpm) and qq_log
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Raw PCA.",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "top 40 PC1 loadings.",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "PCA(log2(cpm())).",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Raw TSNE.",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "TSNE(log2(cpm())).",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Raw QQ, log scale.",
                                      startRow = new_row, startCol = new_col)
  new_col <- 1
  new_row <- new_row + 1
  pca_plot <- metrics[["pc_plot"]]
  pca_topn <- metrics[["pc_loadplot"]]
  pca_table <- metrics[["pc_table"]]
  tsne_plot <- metrics[["tsne_plot"]]
  tsne_table <- metrics[["tsne_table"]]
  try_result <- xlsx_insert_png(pca_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "12_pcaplot", savedir = excel_basename, fancy_type = "svg")
  if ("try-error" %in% class(try_result)) {
    warning("Failed to add the initial PCA plot.")
  } else {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(pca_topn[["plot"]], wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "13_pctopn", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  tmp_data <- sm(normalize(se, transform = "log2", convert = "cpm", filter = filter,
                                ...))
  rpca <- plot_pca(tmp_data,
                   ...)
  rtsne <- plot_tsne(tmp_data,
                     ...)
  rspca_plot <- rpca[["plot"]]
  rtsne_plot <- rtsne[["plot"]]
  rpca_table <- rpca[["residual_df"]]
  rtsne_table <- rtsne[["residual_df"]]
  rm(tmp_data)
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(rspca_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "14_norm_pcaplot", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(tsne_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "15_tsneplot", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(rtsne_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "16_rtsneplot", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  qq_plot <- metrics[["qqlog"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(qq_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "17_qqlog", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- 1

  violin_plot <- NULL
  pct_plot <- NULL
  ## Violin plots
  if (isTRUE(violin)) {
    filt <- sm(normalize(se, filter = "simple"))
    do_varpart <- TRUE
    full_model <- as.formula("~ condition + batch")
    reduced_model <- as.formula("~ condition")
    data_full_model <- try(stats::model.matrix.default(full_model, data = pData(filt)), silent = TRUE)
    data_reduced_model <- NULL
    full_model_columns <- 0
    reduced_model_columns <- 0
    full_model_rank <- 0
    reduced_model_rank <- 0
    varpart_factors <- c("condition")
    fstring <- "~ condition + batch"
    if ("try-error" %in% class(data_full_model)) {
      do_varpart <- FALSE
      message("The expressionset has a minimal or missing set of conditions/batches.")
    } else {
      data_reduced_model <-  stats::model.matrix.default(reduced_model, data = pData(se))
      full_model_columns <- ncol(data_full_model)
      reduced_model_columns <- ncol(data_reduced_model)
      full_model_rank <- qr(data_full_model)[["rank"]]
      reduced_model_rank <- qr(data_reduced_model)[["rank"]]
      varpart_factors <- c("condition", "batch")
    }

    varpart_raw <- NULL
    if (full_model_rank < full_model_columns) {
      message("This expressionset does not support lmer with condition+batch")
      if (reduced_model_rank < reduced_model_columns) {
        message("This expressionset also does not support lmer with just condition!")
        do_varpart <- FALSE
      } else {
        varpart_factors <- "condition"
        fstring <- "~ condition"
      }
    }
    if (isTRUE(do_varpart)) {
      varpart_raw <- sm(suppressWarnings(try(simple_varpart(filt, fstring = fstring),
                                             silent = TRUE)))
    }
    if (! "try-error" %in% class(varpart_raw)) {
      varpart_raw <- NULL
      do_varpart <- FALSE
    }
    if (!is.null(varpart_raw)) {
      violin_plot <- varpart_raw[["partition_plot"]]
      new_row <- new_row + plot_rows + 2
      new_col <- 1
      try_result <- xlsx_insert_png(violin_plot, wb = wb, sheet = sheet, width = plot_dim,
                                    height = plot_dim, start_col = new_col, start_row = new_row,
                                    plotname = "18_violin", savedir = excel_basename)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
      new_col <- new_col + plot_cols + 1

      pct_plot <- varpart_raw[["percent_plot"]]
      try_result <- xlsx_insert_png(pct_plot, wb = wb, sheet = sheet, width = plot_dim,
                                    height = plot_dim, start_col = new_col, start_row = new_row,
                                    plotname = "19_pctvar", savedir = excel_basename)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
    }
  }

  ## PCA table
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Raw PCA res.",
                                      startRow = new_row, startCol = new_col)
  new_row <- new_row + 1
  xls_result <- write_xlsx(data = metrics[["pc_summary"]], wb = wb, rownames = FALSE,
                           sheet = sheet, start_col = new_col, start_row = new_row)
  new_col <- xls_result[["end_col"]] + 6
  new_row <- new_row - 1
  written_name <- openxlsx::writeData(wb, sheet, "Raw PCA table.",
                                      startRow = new_row, startCol = new_col)
  new_row <- new_row + 1
  xls_result <- write_xlsx(data = metrics[["pc_table"]], wb = wb, rownames = FALSE,
                           sheet = sheet, start_row = new_row, start_col = new_col)

  ## Move on to the next sheet, normalized data
  mesg("Writing the normalized reads.")
  sheet <- "norm_data"
  new_col <- 1
  new_row <- 1
  ## Perform a quick query to see if sva will explode on this data.
  test_norm <- normalize(se = se, transform = transform,
                              convert = convert, filter = filter)
  test_zeros <- sum(rowSums(exprs(test_norm)) == 0)
  if (test_zeros > 0) {
    actual_filter <- "simple"
  } else {
    actual_filter <- filter
  }
  norm_data <- sm(normalize(se = se, transform = transform,
                                 convert = convert, batch = batch,
                                 filter = actual_filter,
                                 ...))
  norm_reads <- exprs(norm_data)
  info <- fData(norm_data)
  read_info <- merge(norm_reads, info, by = "row.names")
  title <- what_happened(norm_data)
  xls_result <- write_xlsx(wb = wb, data = read_info, rownames = FALSE,
                           start_row = new_row, start_col = new_col, sheet = sheet, title = title)

  ## Potentially useful for proteomics data and subtracted data.
  if (!is.null(color_na)) {
    na_style <- openxlsx::createStyle(fontColour = color_na)
    nas <- se[["na_values"]]
    for (col in colnames(nas)) {
      row_definition <- which(nas[[col]]) + 2
      col_idx <- colnames(read_info) == col
      col_definition <- which(col_idx)
      colored <- openxlsx::addStyle(wb, sheet = sheet, na_style,
                                    rows = row_definition, cols = col_definition)
    }
  }

  ## Graphs of the normalized data
  mesg("Graphing the normalized reads.")
  sheet <- "norm_graphs"
  newsheet <- try(openxlsx::addWorksheet(wb, sheetName = sheet))
  norm_metrics <- sm(graph_metrics(norm_data, qq = do_qq, gene_heat = do_sample_heat,
                                   ...))
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Raw PCA res.",
                                      startRow = new_row, startCol = new_col)
  ## Start with library sizes.
  written_name <- openxlsx::writeData(wb, sheet, "Legend.",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet, "Normalized library sizes.",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Non-zero genes.",
                                      startRow = new_row, startCol = new_col)
  new_col <- 1
  new_row <- new_row + 1
  new_plot <- norm_metrics[["legend"]]
  try_result <- xlsx_insert_png(new_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  nlibsize_plot <- norm_metrics[["libsize"]]
  try_result <- xlsx_insert_png(nlibsize_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "20_nlibsize", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  ## Same row, non-zero plot
  new_col <- new_col + plot_cols + 1
  nnzero_plot <- norm_metrics[["nonzero"]]
  try_result <- xlsx_insert_png(nnzero_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "21_nnzero", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1

  ## Visualize distributions
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Normalized data density plot.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Normalized Boxplot.",
                      startRow = new_row, startCol = new_col)
  new_col <- 1
  ndensity_plot <- norm_metrics[["density"]]
  new_row <- new_row + 1
  try_result <- xlsx_insert_png(ndensity_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "22_ndensity", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  nboxplot_plot <- norm_metrics[["boxplot"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(nboxplot_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "23_nboxplot", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  ntopn_plot <- norm_metrics[["topnplot"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(ntopn_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "24_nboxplot", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- 1

  ## Move down next set of rows, heatmaps
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Normalized correlation heatmap.",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Normalized distance heatmap.",
                                      startRow = new_row, startCol = new_col)
  new_col <- 1
  ncorheat_plot <- norm_metrics[["corheat"]]
  new_row <- new_row + 1
  try_result <- xlsx_insert_png(ncorheat_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "25_ncorheat", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  ndisheat_plot <- norm_metrics[["disheat"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(ndisheat_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "26_ndisheat", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  if (isTRUE(sample_heat)) {
    sampleheat_plot <- plot_sample_heatmap(norm_data)
    new_col <- new_col + plot_cols + 1
    try_result <- xlsx_insert_png(sampleheat_plot[["plot"]], wb = wb, sheet = sheet, width = plot_dim,
                                  height = plot_dim, start_col = new_col, start_row = new_row,
                                  plotname = "26a_sampleheat", savedir = excel_basename)
    if (! "try-error" %in% class(try_result)) {
      image_files <- c(image_files, try_result[["filename"]])
    }
  }
  new_col <- 1

  ## SM plots
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  openxlsx::writeData(wb, sheet = sheet, x = "Normalized standard median correlation.",
                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  openxlsx::writeData(wb, sheet = sheet, x = "Normalized standard distance correlation.",
                      startRow = new_row, startCol = new_col)
  new_col <- 1
  nsmc_plot <- norm_metrics[["smc"]]
  new_row <- new_row + 1
  try_result <- xlsx_insert_png(nsmc_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "27_nsmc", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  nsmd_plot <- norm_metrics[["smd"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(nsmd_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "28_nsmd", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- 1

  ## PCA and qq_log
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Normalized PCA.",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Top 40 PC1 loadings.",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Normalized TSNE.",
                                      startRow = new_row, startCol = new_col)
  new_col <- new_col + plot_cols + 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Normalized QQ, log scale.",
                                      startRow = new_row, startCol = new_col)
  new_col <- 1
  npca_plot <- norm_metrics[["pc_plot"]]
  npc_topnplot <- norm_metrics[["pc_loadplot"]]
  ntsne_plot <- norm_metrics[["tsne_plot"]]
  npca_table <- norm_metrics[["pc_table"]]
  ntsne_table <- norm_metrics[["tsne_table"]]
  new_row <- new_row + 1
  try_result <- xlsx_insert_png(npca_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "29_npcaplot", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(npc_topnplot[["plot"]], wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "30_npcloadplot", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(ntsne_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "31_ntsneplot", savedir = excel_basename, fancy_type = "svg")
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }
  nqq_plot <- norm_metrics[["qqlog"]]
  new_col <- new_col + plot_cols + 1
  try_result <- xlsx_insert_png(nqq_plot, wb = wb, sheet = sheet, width = plot_dim,
                                height = plot_dim, start_col = new_col, start_row = new_row,
                                plotname = "32_nqqplot", savedir = excel_basename)
  if (! "try-error" %in% class(try_result)) {
    image_files <- c(image_files, try_result[["filename"]])
  }

  new_col <- 1

  ## Violin plots
  nvarpart_plot <- NULL
  npct_plot <- NULL
  if (isTRUE(violin) && isTRUE(do_varpart)) {
    varpart_norm <- suppressWarnings(try(simple_varpart(norm_data, fstring = fstring),
                                         silent = TRUE))
    if (! "try-error" %in% class(varpart_norm)) {
      nvarpart_plot <- varpart_norm[["partition_plot"]]
      new_row <- new_row + plot_rows + 2
      new_col <- 1
      try_result <- xlsx_insert_png(nvarpart_plot, wb = wb, sheet = sheet, width = plot_dim,
                                    height = plot_dim, start_col = new_col, start_row = new_row,
                                    plotname = "33_nviolin", savedir = excel_basename)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
      new_col <- new_col + plot_cols + 1
      npct_plot <- varpart_norm[["percent_plot"]]
      try_result <- xlsx_insert_png(npct_plot, wb = wb, sheet = sheet, width = plot_dim,
                                    height = plot_dim, start_col = new_col, start_row = new_row,
                                    plotname = "34_npctplot", savedir = excel_basename)
      if (! "try-error" %in% class(try_result)) {
        image_files <- c(image_files, try_result[["filename"]])
      }
    }
  }

  ## PCA table
  new_row <- new_row + plot_rows + 2
  new_col <- 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Normalized PCA res.",
                                      startRow = new_row, startCol = new_col)
  new_row <- new_row + 1
  xls_result <- write_xlsx(data = norm_metrics[["pc_summary"]], wb = wb, rownames = FALSE,
                           sheet = sheet, start_col = new_col, start_row = new_row)
  new_col <- xls_result[["end_col"]] + 6
  new_row <- new_row - 1
  written_name <- openxlsx::writeData(wb, sheet = sheet, x = "Normalized PCA table.",
                      startRow = new_row, startCol = new_col)
  new_row <- new_row + 1
  xls_result <- write_xlsx(data = norm_metrics[["pc_table"]], wb = wb, sheet = sheet,
                           rownames = FALSE, start_col = new_col, start_row = new_row)

  ## Give a median-by-factor accounting of the data
  mesg("Writing the median reads by factor.")
  sheet <- "median_data"
  new_col <- 1
  new_row <- 1
  median_data <- sm(median_by_factor(norm_data, fun = med_or_mean,
                                     fact = norm_data[["conditions"]]))
  med <- median_data[["medians"]]
  colnames(med) <- paste0(med_or_mean, "_", colnames(med))
  cv <- median_data[["cvs"]]
  colnames(cv) <- paste0("cv_", colnames(cv))
  median_data <- merge(med, cv, by = "row.names")
  rownames(median_data) <- median_data[["Row.names"]]
  median_data[["Row.names"]] <- NULL
  median_data_merged <- data.frame()
  if (merge_order == "annot_first") {
    median_data_merged <- merge(info, median_data, by.x = "row.names", by.y = "row.names")
  } else {
    median_data_merged <- merge(median_data, info, by.x = "row.names", by.y = "row.names")
  }
  rownames(median_data_merged) <- median_data_merged[["Row.names"]]
  median_data_merged[["Row.names"]] <- NULL
  xls_result <- write_xlsx(wb, data = median_data_merged, start_row = new_row, start_col = new_col,
                           rownames = TRUE, sheet = sheet, title = "Median Reads by factor.")

  ## Save the result
  save_result <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  retlist <- list(
    "excel" = excel,
    "save" = save_result,
    "legend" = legend,
    "annotations" = info,
    "raw_reads" = reads,
    "design" = annot,
    "legend_plot" = legend_plot,
    "raw_libsize" = libsize_plot,
    "raw_nonzero" = nonzero_plot,
    "raw_density" = density_plot,
    "raw_cv" = cv_plot,
    "raw_boxplot" = boxplot_plot,
    "raw_corheat" = corheat_plot,
    "raw_disheat" = disheat_plot,
    "raw_smc" = smc_plot,
    "raw_smd" = smd_plot,
    "raw_pctopn" = pca_topn,
    "raw_pca" = pca_plot,
    "raw_pca_table" = pca_table,
    "raw_tsne" = tsne_plot,
    "raw_tsne_table" = tsne_table,
    "raw_scaled_pca" = rspca_plot,
    "raw_scaled_pca_table" = rpca_table,
    "raw_scaled_tsne" = rtsne_plot,
    "raw_scaled_tsne_table" = rtsne_table,
    "raw_qq" = qq_plot,
    "raw_violin" = violin_plot,
    "raw_percent" = pct_plot,
    "norm_reads" = norm_reads,
    "norm_libsize" = nlibsize_plot,
    "norm_nonzero" = nnzero_plot,
    "norm_density" = ndensity_plot,
    "norm_boxplot" = nboxplot_plot,
    "norm_corheat" = ncorheat_plot,
    "norm_disheat" = ndisheat_plot,
    "norm_smc" = nsmc_plot,
    "norm_smd" = nsmd_plot,
    "norm_pca" = npca_plot,
    "norm_pctopn" = npc_topnplot,
    "norm_pca_table" = npca_table,
    "norm_tsne" = ntsne_plot,
    "norm_tsne_table" = ntsne_table,
    "norm_qq" = nqq_plot,
    "norm_violin" = nvarpart_plot,
    "norm_pct" = npct_plot,
    "medians" = median_data,
    "saved" = save_result
  )
  for (img in image_files) {
    removed <- try(suppressWarnings(file.remove(img)), silent = TRUE)
  }
  class(retlist) <- "written_se"
  return(retlist)
}

#' Print the result from write_se.
#'
#' @param x List containing all the many plots, the dataframes, etc.
#' @param ... Other args to match the generic.
#' @export
print.written_se <- function(x, ...) {
  result_string <- glue("The result from write_se() sent to:
{x[['excel']]}")
  message(result_string)
  return(invisible(x))
}



## EOF
