# expt.r: Create and play with expressionsets.  This actually is an
## implementation of a non-S4 superclass of the expressionSet.  I did this
## because I was having annoying problems creating expressionSets.  Thus, >90%
## of the logic here is intended to simplify and standardize that process.

## Ensure imports are loaded from 01_hpgltools.R
#' @include 01_hpgltools.R
NULL

#' Take two expressionsets and smoosh them together.
#'
#' Because of the extra sugar I added to expressionSets, the combine() function
#' needs a little help when combining expts.  Notably, the information from
#' tximport needs some help.
#'
#' @param expt1 First expt object.
#' @param expt2 Second expt object.
#' @param condition Column with which to reset the conditions.
#' @param all_x Keep all of the first expt's annotations/counts if there are mismatches?
#' @param all_y Keep all the second expt's annotations/counts if there are mismatches?
#' @param batch Column with which to reset the batches.
#' @param merge_meta Merge the metadata when they mismatch?  This should perhaps default to TRUE.
#' @return Larger expt.
#' @seealso [set_expt_batches()] [set_expt_conditions()] [set_expt_colors()]
#'  [set_expt_genenames()] [set_expt_samplenames()] [subset_expt()] [create_expt()]
#' @examples
#'  \dontrun{
#'   ## I am trying to get rid of all my dontrun sections, but I don't have two
#'   ## expressionsets to combine.
#'   expt1 <- create_expt(first_meta)
#'   expt2 <- create_expt(second_meta)
#'   combined <- combine_expts(expt1, expt2, merge_meta = TRUE)
#' }
#' @export
combine_expts <- function(expt1, expt2, condition = "condition", all_x = TRUE, all_y = TRUE,
                          batch = "batch", merge_meta = TRUE) {
  exp1 <- expt1[["expressionset"]]
  exp2 <- expt2[["expressionset"]]
  fData(exp2) <- fData(exp1)
  ##testthat::expect_equal(rownames(exprs(exp1)), rownames(exprs(exp2)))
  if (isTRUE(merge_meta)) {
    design1 <- pData(exp1)
    d1_rows <- seq_len(nrow(design1))
    design2 <- pData(exp2)
    both <- as.data.frame(rbindlist(list(design1, design2), fill = TRUE))
    d2_rows <- (nrow(design1) + 1):nrow(both)
    new_design1 <- both[d1_rows, ]
    rownames(new_design1) <- rownames(design1)
    new_design2 <- both[d2_rows, ]
    rownames(new_design2) <- rownames(design2)
    pData(exp1) <- new_design1
    pData(exp2) <- new_design2
  }

  ## new <- a4Base::combineTwoExpressionSet(exp1, exp2)
  new <- Biobase::combine(exp1, exp2)
  expt1[["expressionset"]] <- new
  expt1[["conditions"]] <- pData(expt1)[["condition"]]
  names(expt1[["conditions"]]) <- rownames(pData(expt1))
  expt1[["batches"]] <- pData(expt1)[["batch"]]
  names(expt1[["batches"]]) <- rownames(pData(expt1))
  expt1[["colors"]] <- c(expt1[["colors"]], expt2[["colors"]])
  expt1 <- set_expt_conditions(expt1, fact = condition)

  if (!is.null(expt1[["tximport"]])) {
    raw1 <- expt1[["tximport"]][["raw"]]
    raw2 <- expt2[["tximport"]][["raw"]]
    merged <- merge(raw1[["abundance"]], raw2[["abundance"]], by = "row.names",
                    all.x = all_x, all.y = all_y)
    rownames(merged) <- merged[["Row.names"]]
    merged[["Row.names"]] <- NULL
    expt1[["tximport"]][["raw"]][["abundance"]] <- merged
    merged <- merge(raw1[["counts"]], raw2[["counts"]], by = "row.names",
                    all.x = all_x, all.y = all_y)
    rownames(merged) <- merged[["Row.names"]]
    merged[["Row.names"]] <- NULL
    expt1[["tximport"]][["raw"]][["counts"]] <- merged
    merged <- merge(raw1[["length"]], raw2[["length"]], by = "row.names",
                    all.x = all_x, all.y = all_y)
    rownames(merged) <- merged[["Row.names"]]
    merged[["Row.names"]] <- NULL
    expt1[["tximport"]][["raw"]][["length"]] <- merged
    scaled1 <- expt1[["tximport"]][["scaled"]]
    scaled2 <- expt2[["tximport"]][["scaled"]]
    merged <- merge(scaled1[["abundance"]], scaled2[["abundance"]], by = "row.names",
                    all.x = all_x, all.y = all_y)
    rownames(merged) <- merged[["Row.names"]]
    merged[["Row.names"]] <- NULL
    expt1[["tximport"]][["scaled"]][["abundance"]] <- merged
    merged <- merge(scaled1[["counts"]], scaled2[["counts"]], by = "row.names",
                    all.x = all_x, all.y = all_y)
    rownames(merged) <- merged[["Row.names"]]
    merged[["Row.names"]] <- NULL
    expt1[["tximport"]][["scaled"]][["counts"]] <- merged
    merged <- merge(scaled1[["length"]], scaled2[["length"]], by = "row.names")
    rownames(merged) <- merged[["Row.names"]]
    merged[["Row.names"]] <- NULL
    expt1[["tximport"]][["scaled"]][["length"]] <- merged
  }

  ## It appears combining expressionsets can lead to NAs?
  ## I ran into a situation where the expression data of the expressionset somehow got
  ## cast as dataframe.  Let us make sure that doesn't happen here, though I think
  ## there are better places for this type of check.  I did add some checks in my
  ## setMethod calls for exprs<-
  if (class(exprs(expt1))[1] == "data.frame") {
    exprs(expt1) <- as.matrix(exprs(expt1))
  }
  na_idx <- is.na(exprs(expt1))
  if (sum(na_idx) > 0) {
    warning("There appear to be NA values in the new expressionset.  Beware.")
    exprs(expt1)[na_idx] <- 0
  }

  if (isFALSE(all_x) || isFALSE(all_y)) {
    found_first <- rownames(exprs(expt1)) %in% rownames(exprs(exp1))
    shared_first_ids <- rownames(exprs(expt1))[found_first]
    expt1 <- subset_genes(expt1, ids = shared_first_ids, method = "keep")

    found_second <- rownames(exprs(expt1)) %in% rownames(exprs(exp2))
    shared_second_ids <- rownames(exprs(expt1))[found_second]
    expt1 <- subset_genes(expt1, ids = shared_second_ids, method = "keep")
  }
  return(expt1)
}

#' Wrap bioconductor's expressionset to include some extra information.
#'
#' Note: You should just be using create_se().  It does everything the
#' expt does, but better.
#'
#' The primary innovation of this function is that it will check the metadata
#' for columns containing filenames for the count tables, thus hopefully making
#' the collation and care of metadata/counts easier.  For example, I have some
#' data which has been mapped against multiple species.  I can use this function
#' and just change the file_column argument to pick up each species' tables.
#'
#' @param metadata Comma separated file (or excel) describing the samples with
#'  information like condition, batch, count_filename, etc.
#' @param gene_info Annotation information describing the rows of the data set,
#'  this often comes from a call to import.gff() or biomart or organismdbi.
#' @param count_dataframe If one does not wish to read the count tables from the
#'  filesystem, they may instead be fed as a data frame here.
#' @param sanitize_rownames Clean up weirdly written gene IDs?
#' @param sample_colors List of colors by condition, if not provided it will
#'  generate its own colors using colorBrewer.
#' @param title Provide a title for the expt?
#' @param notes Additional notes?
#' @param include_type I have usually assumed that all gff annotations should be
#'  used, but that is not always true, this allows one to limit to a specific
#'  annotation type.
#' @param countdir Directory containing count tables.
#' @param include_gff Gff file to help in sorting which features to keep.
#' @param file_column Column to use in a gene information dataframe for
#' @param id_column Column which contains the sample IDs.
#' @param savefile Rdata filename prefix for saving the data of the resulting
#'  expt.
#' @param low_files Explicitly lowercase the filenames when searching the
#'  filesystem?
#' @param handle_na How does one wish to deal with NA values in the data?
#' @param researcher Used to make the creation of gene sets easier,
#'  set the researcher tag.
#' @param study_name Ibid, but set the study tag.
#' @param file_type Explicitly state the type of files containing the
#'  count data.  I have code which autodetects the method used to
#'  import count data, this short-circuits it.
#' @param annotation_name Ibid, but set the orgdb (or other annotation) instance.
#' @param tx_gene_map Dataframe of transcripts to genes, primarily for
#'  tools like salmon.
#' @param feature_type Make explicit the type of feature used so it
#'  may be printed later.
#' @param ignore_tx_version When using tximport, one may strictly
#'  match the transcript versions, or not.
#' @param keep_underscore Sanitize out underscores from the columns?
#' @param ... More parameters are fun!
#' @return experiment an expressionset
#' @seealso [Biobase] [cdm_expt_rda] [example_gff] [sb_annot] [sb_data] [extract_metadata()]
#'  [set_expt_conditions()] [set_expt_batches()] [set_expt_samplenames()] [subset_expt()]
#'  [set_expt_colors()] [set_expt_genenames()] [tximport] [load_annotations()]
#' @examples
#'  ## There are lots of other ways to use this, for example:
#'  \dontrun{
#'   new_experiment <- create_expt(metadata = "some_csv_file.csv", gene_info = gene_df)
#'   ## Remember that this depends on an existing data structure of gene annotations.
#'   meta <- extract_metadata("some_supplementary_materials_xls_file_I_downloaded.xls")
#'   another_expt <- create_expt(metadata = meta, gene_info = annotations, count_dataframe = df_I_downloaded)
#'  }
#' @import Biobase
#' @exportClass expt
#' @export
create_expt <- function(metadata = NULL, gene_info = NULL, count_dataframe = NULL,
                        sanitize_rownames = TRUE, sample_colors = NULL, title = NULL,
                        notes = NULL, include_type = "all", countdir = NULL,
                        include_gff = NULL, file_column = "file", id_column = NULL,
                        savefile = NULL, low_files = FALSE, handle_na = "drop",
                        researcher = "elsayed", study_name = NULL, file_type = NULL,
                        annotation_name = "org.Hs.eg.db", tx_gene_map = NULL,
                        feature_type = "gene", ignore_tx_version = TRUE,
                        keep_underscore = TRUE, ...) {
  arglist <- list(...)  ## pass stuff like sep=, header=, etc here

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
    title <- "This is an expt class."
  }
  if (is.null(notes)) {
    notes <- glue("Created on {date()}.
")
  }

  ## Palette for colors when auto-chosen
  chosen_palette <- "Dark2"
  if (!is.null(arglist[["palette"]])) {
    chosen_palette <- arglist[["palette"]]
  }
  file_suffix <- ".count.gz"
  if (!is.null(arglist[["file_suffix"]])) {
    file_suffix <- arglist[["file_suffix"]]
  }
  file_prefix <- ""
  if (!is.null(arglist[["file_prefix"]])) {
    file_prefix <- arglist[["file_prefix"]]
  }
  gff_type <- "all"
  if (!is.null(arglist[["include_type"]])) {
    gff_type <- arglist[["include_type"]]
  }

  if (is.null(id_column)) {
    id_column <- "sampleid"
  } else {
    id_column <- tolower(id_column)
    if (isTRUE(keep_underscore)) {
      id_column <- gsub(pattern = "[^_[:^punct:]]", replacement = "", x = id_column, perl = TRUE)
    } else {
      id_column <- gsub(pattern = "[[:punct:]]", replacement = "", x = id_column)
    }
  }
  file_column <- tolower(file_column)
  if (isTRUE(keep_underscore)) {
    file_column <- gsub(pattern = "[^_[:^punct:]]", replacement = "", x = file_column, perl = TRUE)
  } else {
    file_column <- gsub(pattern = "[[:punct:]]", replacement = "", x = file_column)
  }

  round <- FALSE
  if (!is.null(arglist[["round"]])) {
    round <- arglist[["round"]]
  }

  ## Read in the metadata from the provided data frame, csv, or xlsx.
  message("Reading the sample metadata.")
  sample_definitions <- extract_metadata(metadata, id_column = id_column,
                                         keep_underscore = keep_underscore,
                                         ...)
  ## Add an explicit removal of the file column if the option file_column is NULL.
  ## This is a just in case measure to avoid conflicts.
  if (is.null(file_column)) {
    if (!is.null(metadata[["file"]])) {
      message("file_column is NULL, removing the column named 'file'.")
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
      message("The metadata row names are: ",
              toString(sort(rownames(sample_definitions))))
      stop("The count table column names are not the same as the sample definition row names.")
    }
    all_count_tables <- as.data.table(count_dataframe, keep.rownames = "rownames")
    ## If neither of these cases is true, start looking for the files in the
    ## processed_data/ directory
  } else if (is.null(sample_definitions[[file_column]])) {
    message("The file column: ", file_column, " is not in the set of sample columns.")
    print(colnames(sample_definitions))
    stop("This requires either a count dataframe/matrix or column containing the filenames.")
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
                              tx_gene_map = tx_gene_map, file_type = file_type,
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
  ## While we are removing stuff...
  ## I have had a couple data sets with incomplete counts, get rid of those rows
  ## before moving on.

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
  ## Features like exon:alicethegene-1 are annoying and entirely too common in TriTrypDB data
  if (isTRUE(sanitize_rownames)) {
    mesg("Sanitizing rownames, if something unexpected happens later, look here first.")
    all_count_tables[["rownames"]] <- gsub(
      x = all_count_tables[["rownames"]],
      pattern = "^exon:", replacement = "") %>%
      gsub(pattern = "^gene:", replacement = "")

    all_count_tables[["rownames"]] <- make.names(gsub(
      x = all_count_tables[["rownames"]],
      pattern = ":\\d+", replacement = ""), unique = TRUE)
  }

  ## There is an important caveat here!!
  ## as.data.table(stuff, keep.rownames='column') will change the
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

  ## Try a couple different ways of getting gene-level annotations into the expressionset.
  annotation <- NULL
  if (is.null(gene_info)) {
    ## Including, if all else fails, just grabbing the gene names from the count tables.
    if (is.null(include_gff)) {
      gene_info <- as.data.table(all_count_tables[["rownames"]],
                                             keep.rownames = "rownames")
      names(gene_info) <- "rownames"
    } else {
      ## Or reading a gff file.
      message("create_expt(): Reading annotation gff, this is slow.")
      annotation <- load_gff_annotations(gff = include_gff, type = gff_type)
      gene_info <- as.data.table(annotation, keep.rownames = "rownames")
    }
  } else if (class(gene_info)[[1]] == "list" && !is.null(gene_info[["genes"]])) {
    ## In this case, it is using the output of reading a OrgDB instance
    gene_info <- as.data.table(gene_info[["genes"]], keep.rownames = "rownames")
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
    gene_info <- as.data.table(gene_info, keep.rownames = "rownames")
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
  if (!is.null(arglist[["tx_gene_map"]])) {
    tx_gene_map <- arglist[["tx_gene_map"]]
    ## This is wrong, we recast gene_info as a data.table and thus it does not have rownames.
    ## matched_rows <- sum(rownames(gene_info) %in% tx_gene_map[[2]])
    ## Because of the previous line, the following test always sees the match between the tx_gene_map
    ## and gene annotations as 0, and therefore it messed up the rownames of the annotations.
    matched_rows <- gene_info[["rownames"]] %in% tx_gene_map[[2]]
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

  ## There are some shenanigans, Maddy is getting an error on countsdt...
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
  chosen_colors <- generate_se_colors(sample_definitions, sample_colors = sample_colors,
                                      chosen_palette = chosen_palette)

  ## Fill in incomplete tables.
  if (is.null(sample_definitions[["condition"]])) {
    sample_definitions[["condition"]] <- "unknown"
  }
  if (is.null(sample_definitions[["batch"]])) {
    sample_definitions[["batch"]] <- "unknown"
  }
  if (is.null(sample_definitions[["file"]])) {
    sample_definitions[["file"]] <- "null"
  }

  ## Finally, create the ExpressionSet using the counts, annotations, and metadata.
  requireNamespace("Biobase")  ## AnnotatedDataFrame is from Biobase
  metadata <- methods::new("AnnotatedDataFrame", sample_definitions)
  new_samplenames <- try(Biobase::sampleNames(metadata) <- colnames(final_counts))
  if (class(new_samplenames) == "try-error") {
    message("Something is wrong with the sample names for the experimental metadata.")
    message("They must be set to the column names of the count table, which are: ",
            toString(colnames(final_counts)))
    message("If the above column names have .x/.y in them, they may be duplicated,",
            " check your sample sheet.")
    message("The sample definitions are:")
    print(knitr::kable(sample_definitions))
    stop()
  }

  feature_data <- methods::new("AnnotatedDataFrame", final_annotations)
  Biobase::featureNames(feature_data) <- rownames(final_counts)

  ## An expressionset needs to have a Biobase::annotation() in order for
  ## GSEABase to work with it. Reading the documentation, these are primarily
  ## used for naming the type of microarray chip used.
  experiment <- methods::new("ExpressionSet",
                             exprs = final_counts,
                             phenoData = metadata,
                             annotation = annotation_name,
                             featureData = feature_data)

  Biobase::notes(experiment) <- toString(notes)
  ## These entries in new_expt are intended to maintain a record of
  ## the transformation status of the data, thus if we now call
  ## normalize_expt() it should change these.
  ## Therefore, if we call a function like DESeq() which requires
  ## non-log2 counts, we can check these values and convert accordingly

  ## expt <- sm(subset_expt(experiment)) ## I think this is spurious now.

  if (is.null(study_name)) {
    study_name <- basename(getwd())
  }
  expt <- list(
    "expressionset" = experiment,
    "annotation" = annotation,
    "colors" = chosen_colors,
    ## "design" = sample_definitions,
    "feature_type" = feature_type,
    "gff_file" = include_gff,
    "libsize" = colSums(exprs(experiment)),
    "notes" = toString(notes),
    "researcher" = researcher,
    "study" = study_name,
    "title" = title,
    "tximport" = tximport_data)
  ## the 'state' slot in the expt is used to keep track of how the data is modified over time.
  starting_state <- list(
    "batch" = "raw",
    "conversion" = "raw",
    "filter" = "raw",
    "normalization" = "raw",
    "transform" = "raw")
  state(expt) <- starting_state
  ## Just in case there are condition names which are not used.
  ## Ditto here, this should not be needed.
  ## expt[["conditions"]] <- droplevels(as.factor(sample_definitions[, "condition"]))
  ## This might be redundant, but it ensures that no all-numeric conditions exist.
  ## expt[["conditions"]] <- gsub(pattern = "^(\\d+)$", replacement = "c\\1", x = expt[["conditions"]])
  ##names(expt[["conditions"]]) <- rownames(sample_definitions)
  ## Ditto for batches
  ## expt[["batches"]] <- droplevels(as.factor(sample_definitions[, "batch"]))
  ## expt[["batches"]] <- gsub(pattern = "^(\\d+)$", replacement = "b\\1", x = expt[["batches"]])
  ## names(expt[["batches"]]) <- rownames(sample_definitions)
  ## Keep a backup of the library sizes for limma.
  if (sum(expt[["libsize"]] == 0) > 0) {
    zero_idx <- expt[["libsize"]] == 0
    zero_samples <- names(expt[["libsize"]])[zero_idx]
    warning("The following samples have no counts! ", zero_samples)
  }
  names(expt[["libsize"]]) <- rownames(sample_definitions)
  ## Save the chosen colors
  names(expt[["colors"]]) <- rownames(sample_definitions)
  class(expt) <- "expt"
  ## Save an rdata file of the expressionset.
  if (is.null(savefile)) {
    if ("character" %in% class(metadata)) {
      name <- paste0(gsub(x = basename(metadata), pattern = "^(.*)\\..*",
                          replacement = "\\1"), ".rda")
    } else {
      message("Saving the expressionset to 'expt.rda'.")
      savefile <- "expt.rda"
    }
  }
  if (!dir.exists(dirname(savefile))) {
    created <- dir.create(dirname(savefile), recursive = TRUE)
  }
  save_result <- try(save(expt, file = savefile), silent = TRUE)
  if (class(save_result) == "try-error") {
    warning("Saving the expt object failed, perhaps you do not have permissions?")
  }
  message("The final expressionset has ", nrow(exprs(expt)),
          " features and ", ncol(exprs(expt)), " samples.")
  return(expt)
}
setOldClass("expt")

#' Modified print function for an expt.
#'
#' I am trying to understand how R collates functions.
#' @param x List from create_expt containing the expressionSet,
#'  annotation data, batches, conditions, colors, libsizes, etc.
#' @param ... Other args to match the generic.
#' @export
print.expt <- function(x, ...) {
  feature_type <- "genes"
  if (is.null(x[["feature_type"]])) {
    feature_type <- x[["feature_type"]]
  }
  num_annotations <- ncol(fData(x))
  num_meta <- ncol(pData(x))
  num_samples <- ncol(exprs(x))
  num_rows <- nrow(exprs(x))
  state <- what_happened(x)
  condition_set <- toString(levels(as.factor(pData(x)[["condition"]])))
  ## I have a typeographic error somewhere in the following multi-line statement:
  ## and I cannot seem to find it.
  ##summary_string <- glue("An expressionSet containing experiment with {num_rows}
##{feature_type} and {num_samples} samples. There are {num_meta} metadata columns and
##{num_annotations} annotation columns; the primary condition is comprised of:
##{condition_set}.
##Its current state is: {state}.")
  message("A modified expressionSet containing ", num_rows, " ",
          feature_type, " and ", num_samples, " sample. There are ",
          num_meta, " metadata columns and ", num_annotations, " annotation columns.
The primary condition is comprised of:
", condition_set, ".
Its current state is: ", state, ".")
  return(invisible(x))
}

#' Synchronize the extra elements of an expt with a new expressionset.
#'
#' @param expt Modified/new expt
#' @param previous Optional previous state to use as a template.
#' @param ... Parameters used to fill in other optional slots.
synchronize_expt <- function(expt, previous = NULL, ...) {
  arglist <- list(...)
  expt[["batches"]] <- pData(expt)[["batch"]]
  expt[["conditions"]] <- pData(expt)[["condition"]]
  expt[["libsize"]] <- colSums(exprs(expt))
  if (is.null(previous)) {
    expt[["colors"]] <- generate_se_colors(pData(expt), cond_column = "condition")
    expt[["state"]] <- list(
      "batch" = "raw",
      "conversion" = "raw",
      "filter" = "raw",
      "normalization" = "raw",
      "transform" = "raw")
  } else {
    expt[["researcher"]] <- previous[["researcher"]]
    expt[["notes"]] <- previous[["notes"]]
    expt[["feature_type"]] <- previous[["feature_type"]]
    expt[["gff_file"]] <- previous[["gff_file"]]
    expt[["study"]] <- previous[["study"]]
    expt[["title"]] <- previous[["title"]]
    expt[["state"]] <- previous[["state"]]

    if (nrow(pData(expt)) == length(previous[["conditions"]])) {
      expt[["colors"]] <- previous[["colors"]]
    } else {
      old_colors <- previous[["colors"]]
      new_ids <- rownames(pData(expt))
      new_colors <- old_colors[new_ids]
      expt[["colors"]] <- new_colors

      df <- expt[["tximport"]][["raw"]][["abundance"]][new_ids, ]
      expt[["tximport"]][["raw"]][["abundance"]] <- df
      df <- expt[["tximport"]][["raw"]][["counts"]][new_ids, ]
      expt[["tximport"]][["raw"]][["counts"]] <- df
      df <- expt[["tximport"]][["raw"]][["length"]][new_ids, ]
      expt[["tximport"]][["raw"]][["length"]] <- df

      df <- expt[["tximport"]][["scaled"]][["abundance"]][new_ids, ]
      expt[["tximport"]][["scaled"]][["abundance"]] <- df
      df <- expt[["tximport"]][["scaled"]][["counts"]][new_ids, ]
      expt[["tximport"]][["scaled"]][["counts"]] <- df
      df <- expt[["tximport"]][["scaled"]][["length"]][new_ids, ]
      expt[["tximport"]][["scaled"]][["length"]] <- df
    }
  }
  if (!is.null(arglist[["researcher"]])) {
    expt[["researcher"]] <- arglist[["researcher"]]
  }
  if (!is.null(arglist[["notes"]])) {
    expt[["notes"]] <- paste(expt[["notes"]], arglist[["notes"]])
  }
  if (!is.null(arglist[["feature_type"]])) {
    expt[["feature_type"]] <- arglist[["feature_type"]]
  }
  if (!is.null(arglist[["gff_file"]])) {
    expt[["gff_file"]] <- arglist[["gff_file"]]
  }
  if (!is.null(arglist[["study"]])) {
    expt[["study"]] <- arglist[["study"]]
  }
  return(expt)
}
setOldClass("expt")

#' Do features_greater_than() inverted!
#'
#' @param  ... Arguments passed to features_greather_than()
#' @return The set of features less than whatever you would have done with
#'   features_greater_than().
#' @seealso [features_greater_than()]
#' @export
features_less_than <- function(...) {
  features_greater_than(..., inverse = TRUE)
}

#' Count the number of features(genes) greater than x in a data set.
#'
#' Sometimes I am asked how many genes have >= x counts.  Well, here you go.
#'
#' Untested as of 2016-12-01 but used with Lucia.  I think it would be interesting to iterate
#' this function from small to large cutoffs and plot how the number of kept genes decreases.
#'
#' @param data Dataframe/exprs/matrix/whatever of counts.
#' @param cutoff Minimum number of counts.
#' @param hard Greater-than is hard, greater-than-equals is not.
#' @param inverse when inverted, this provides features less than the cutoff.
#' @return A list of two elements, the first comprised of the number of genes
#'   greater than the cutoff, the second with the identities of said genes.
#' @seealso [Biobase]
#' @export
features_greater_than <- function(data, cutoff = 1, hard = TRUE, inverse = FALSE) {
  if ("expt" %in% class(data) ||
        "ExpressionSet" %in% class(data) ||
          "SummarizedExperiment" %in% class(data)) {
    data <- as.data.frame(exprs(data))
  } else {
    data <- as.data.frame(data)
  }
  number_table <- numeric(length = ncol(data))
  names(number_table) <- colnames(data)
  feature_tables <- list()
  for (col in seq_along(colnames(data))) {
    column_name <- colnames(data)[col]
    column_data <- data[[column_name]]
    num_features <- NULL
    if (isTRUE(hard)) {
      if (isTRUE(inverse)) {
        feature_idx <- column_data < cutoff
      } else {
        feature_idx <- column_data > cutoff
      }
    } else {
      if (isTRUE(inverse)) {
        feature_idx <- column_data <= cutoff
      } else {
        feature_idx <- column_data >= cutoff
      }
    }
    num_features <- sum(feature_idx)
    number_table[[column_name]] <- num_features
    passed_filter <- rownames(data)[feature_idx]
    feature_tables[[column_name]] <- passed_filter
  }
  result <- list(
    "number" = number_table,
    "features" = feature_tables)
  return(result)
}

#' I want an easy way to answer the question: what features are in only condition x?
#'
#' The answer to this lies in a combination of subset_expt() and
#' features_greater_than().
#'
#' @param expt An experiment to query.
#' @param cutoff What is the minimum number of counts required to define
#'  'included.'
#' @param factor What metadata factor to query?
#' @param chosen Either choose a subset or all conditions to query.
#' @return A set of features.
#' @seealso [subset_expt()]
#' @export
features_in_single_condition <- function(expt, cutoff = 2, factor = "condition", chosen = NULL) {
  condition_set <- levels(as.factor(pData(expt)[[factor]]))
  solo_this_list <- list()
  solo_other_list <- list()
  shared_list <- list()
  neither_list <- list()
  for (cond in condition_set) {
    if (!is.null(chosen)) {
      if (! cond %in% chosen) {
        next
      }
    }
    extract_string <- glue("{factor} == '{cond}'")
    single_expt <- sm(subset_expt(expt = expt, subset = extract_string))
    extract_string <- glue("{factor} != '{cond}'")
    others_expt <- sm(subset_expt(expt = expt, subset = extract_string))
    single_data <- exprs(single_expt)
    others_data <- exprs(others_expt)

    single_positive_idx <- rowSums(single_data) >= cutoff
    single_negative_idx <- rowSums(single_data) < cutoff
    single_positive <- single_data[single_positive_idx, ]
    single_negative <- single_data[single_negative_idx, ]

    others_positive_idx <- rowSums(others_data) >= cutoff
    others_negative_idx <- rowSums(others_data) < cutoff
    others_positive <- others_data[others_positive_idx, ]
    others_negative <- others_data[others_negative_idx, ]

    in_both_sets_idx <- rownames(single_positive) %in% rownames(others_positive)
    in_both_sets <- rownames(single_positive)[in_both_sets_idx]

    only_this_set_idx <- rownames(single_positive) %in% rownames(others_negative)
    only_this_set <- rownames(single_positive)[only_this_set_idx]

    only_other_set_idx <- rownames(others_positive) %in% rownames(single_negative)
    only_other_set <- rownames(others_positive)[only_other_set_idx]

    neither_set_this_idx <- rownames(single_negative) %in% rownames(others_negative)
    neither_set_this <- rownames(single_negative)[neither_set_this_idx]
    neither_set_that_idx <- rownames(others_negative) %in% rownames(single_negative)
    neither_set_that <- rownames(others_negative)[neither_set_that_idx]

    shared_list[[cond]] <- in_both_sets
    solo_this_list[[cond]] <- only_this_set
    solo_other_list[[cond]] <- only_other_set
    neither_list[[cond]] <- neither_set_this
  }
  retlist <- list(
    "shared" = shared_list,
    "solo_this" = solo_this_list,
    "solo_other" = solo_other_list,
    "neither" = neither_list
  )
  return(retlist)
}

#' Runs median_by_factor with fun set to 'mean'.
#'
#' @param data Input expt
#' @param fact Metadata factor over which to perform mean().
#' @export
mean_by_factor <- function(data, fact = "condition") {
  median_by_factor(data, fact = fact, fun = "mean")
}

#' Create a data frame of the medians of rows by a given factor in the data.
#'
#' This assumes of course that (like expressionsets) there are separate columns
#' for each replicate of the conditions.  This will just iterate through the
#' levels of a factor describing the columns, extract them, calculate the
#' median, and add that as a new column in a separate data frame.
#'
#' Used in write_expt() as well as a few random collaborations.
#'
#' @param data Data frame, presumably of counts.
#' @param fact Factor describing the columns in the data.
#' @param fun Optionally choose mean or another function.
#' @return Data frame of the medians.
#' @seealso [Biobase] [matrixStats]
#' @export
median_by_factor <- function(data, fact = "condition", fun = "median") {
  if (length(fact) == 1) {
    design <- pData(data)
    fact <- design[[fact]]
    na_idx <- is.na(fact)
    if (sum(na_idx) > 0) {
      fact[na_idx] <- "unknown"
    }
    names(fact) <- rownames(design)
  }
  if ("expt" %in% class(data) ||
        "ExpressionSet" %in% class(data) ||
          "SummarizedExperiment" %in% class(data)) {
    data <- exprs(data)
  }

  medians <- data.frame("ID" = rownames(data), stringsAsFactors = FALSE)
  cvs <- data.frame("ID" = rownames(data), stringsAsFactors = FALSE)
  mins <- data.frame("ID" = rownames(data), stringsAsFactors = FALSE)
  maxs <- data.frame("ID" = rownames(data), stringsAsFactors = FALSE)
  sums <- data.frame("ID" = rownames(data), stringsAsFactors = FALSE)
  data <- as.matrix(data)
  rownames(medians) <- rownames(data)
  rownames(cvs) <- rownames(data)
  rownames(mins) <- rownames(data)
  rownames(maxs) <- rownames(data)
  rownames(sums) <- rownames(data)
  fact <- as.factor(fact)
  used_columns <- c()
  group_indexes <- list()
  samples_per_condition <- c()
  condition_names <- c()
  for (type in levels(fact)) {
    ## columns <- grep(pattern = type, x = fact)
    columns <- as.character(fact) == type
    group_indexes[[type]] <- columns
    samples_per_condition <- c(sum(columns), samples_per_condition)
    condition_names <- c(type, condition_names)
    med <- NULL
    if (sum(columns) < 1) {
      warning("The level ", type, " of the factor has no columns.")
      next
    }
    used_columns <- c(used_columns, type)
    if (sum(columns) == 1) {
      message("The factor ", type, " has only 1 row.")
      med <- as.data.frame(data[, columns], stringsAsFactors = FALSE)
      cv <- as.data.frame(data[, columns], stringsAsFactors = FALSE)
      min <- as.data.frame(data[, columns], stringsAsFactors = FALSE)
      max <- as.data.frame(data[, columns], stringsAsFactors = FALSE)
      sum <- as.data.frame(data[, columns], stringsAsFactors = FALSE)
    } else {
      min <- MatrixGenerics::rowMins(data[, columns], na.rm = TRUE)
      max <- MatrixGenerics::rowMaxs(data[, columns], na.rm = TRUE)
      sum <- rowSums(data[, columns], na.rm = TRUE)
      if (fun == "median") {
        message("The factor ", type, " has ", sum(columns), " rows.")
        med <- matrixStats::rowMedians(data[, columns], na.rm = TRUE)
        cv <- matrixStats::rowMads(data[, columns], na.rm = TRUE)
        cv <- cv / med
        ## I am not really sure if this is appropriate.
        nan_idx <- is.nan(cv)
        cv[nan_idx] <- 0
      } else if (fun == "mean") {
        message("The factor ", type, " has ", sum(columns), " rows.")
        ## Strangely, recently R says BiocGenerics does not export rowMeans.
        ## I stil see its S4 dispatch; but I will assume that R will figure this out.
        ## med <- BiocGenerics::rowMeans(data[, columns], na.rm = TRUE)
        med <- rowMeans(data[, columns], na.rm = TRUE)
        cv <- matrixStats::rowSds(data[, columns], na.rm = TRUE)
        cv <- cv / med
        nan_idx <- is.nan(cv)
        cv[nan_idx] <- 0
      } else {
        stop("I do not understand that function.")
      }
    }
    medians <- cbind(medians, med)
    cvs <- cbind(cvs, cv)
    mins <- cbind(mins, min)
    maxs <- cbind(maxs, max)
    sums <- cbind(sums, sum)
  }
  names(samples_per_condition) <- condition_names
  medians <- medians[, -1, drop = FALSE]
  cvs <- cvs[, -1, drop = FALSE]
  mins <- mins[, -1, drop = FALSE]
  maxs <- maxs[, -1, drop = FALSE]
  sums <- sums[, -1, drop = FALSE]
  ## Sometimes not all levels of the original experimental design are used.
  ## Thus lets make sure to use only those which appeared.
  colnames(medians) <- used_columns
  colnames(cvs) <- used_columns
  colnames(mins) <- used_columns
  colnames(maxs) <- used_columns
  colnames(sums) <- used_columns
  retlist <- list(
    "samples_per_condition" = samples_per_condition,
    "method" = fun,
    "medians" = medians,
    "cvs" = cvs,
    "mins" = mins,
    "maxs" = maxs,
    "sums" = sums,
    "indexes" = group_indexes)
  return(retlist)
}

#' Create a Schizosaccharomyces cerevisiae expt.
#'
#' This just saves some annoying typing if one wishes to make a standard
#' expressionset superclass out of the publicly available fission data set.
#'
#' @param annotation Add annotation data?
#' @param host ensembl host to query.
#' @return Expressionset/expt of fission.
#' @seealso [fission] [create_expt()]
#' @example inst/examples/create_pombe.R
#' @export
make_pombe_expt <- function(annotation = TRUE, host = "nov2020-fungi.ensembl.org") {
  fission <- new.env()
  tt <- sm(requireNamespace("fission"))
  tt <- sm(try(attachNamespace("fission"), silent = TRUE))
  tt <- data(fission, envir = fission)
  ## some minor shenanigans to get around the oddities of loading from data()
  fission <- fission[["fission"]]
  meta <- as.data.frame(fission@colData)
  meta[["condition"]] <- glue("{meta[['strain']]}.{meta[['minute']]}")
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
                        "external_gene_name", "description", "gene_biotype"),
      species = "spombe", overwrite = TRUE))
    if ("try-error" %in% class(pombe_annotations)) {
      warning("There was an error downloading the pombe annotations, this will still return.")
    } else {
      pombe_mart <- pombe_annotations[["mart"]]
      annotations <- pombe_annotations[["annotation"]]
      ## I think ensembl changed the IDs to match and the following line is no longer needed.
      ## rownames(annotations) <- make.names(gsub(pattern = "\\.\\d+$",
      ##                                         replacement = "",
      ##                                         x = rownames(annotations)), unique = TRUE)
    }
  }
  pombe_expt <- create_expt(metadata = meta,
                            count_dataframe = fission_data,
                            gene_info = annotations, annotation_name = "org.Sp.eg.db")
  detach("package:fission")
  return(pombe_expt)
}

#' Remove/keep specifically named genes from an expt.
#'
#' I find subsetting weirdly confusing.  Hopefully this function will allow one
#' to include/exclude specific genes/families based on string comparisons.
#'
#' @param input Expt to filter.
#' @param invert The default is to remove the genes with the semantic strings.
#'  Keep them when inverted.
#' @param topn Take the topn most abundant genes rather than a text based heuristic.
#' @param semantic Character list of strings to search for in the annotation
#'  data.
#' @param semantic_column Column in the annotations to search.
#' @return A presumably smaller expt.
#' @seealso [Biobase]
#' @export
semantic_expt_filter <- function(input, invert = FALSE, topn = NULL,
                                 semantic = c("mucin", "sialidase", "RHS", "MASP", "DGF", "GP63"),
                                 semantic_column = "description") {
  mtrx <- exprs(input)
  annots <- fData(input)
  if (isTRUE(invert)) {
    new_annots <- data.frame()
    new_mtrx <- data.frame()
  } else {
    new_annots <- annots
    new_mtrx <- mtrx
  }
  start_rows <- nrow(mtrx)
  numbers_removed <- 0
  if (is.null(topn)) {
    for (string in semantic) {
      idx <- NULL
      if (isTRUE(invert)) {
        ## Keep the rows which match the ~7 strings above.
        ## For these, we will re-grep the full table each time and just add the matches.
        type <- "Kept"
        if (semantic_column == "rownames") {
          idx <- grepl(pattern = string, x = rownames(annots))
        } else {
          idx <- grepl(pattern = string, x = annots[, semantic_column])
        }
        message("Hit ", sum(idx), " genes for term ", string, ".")
        ## Then, after grepping, just append the matched rows to the new annotations and matrix.
        tmp_annots <- annots[idx, ]
        tmp_mtrx <- mtrx[idx, ]
        new_annots <- rbind(new_annots, tmp_annots)
        new_mtrx <- rbind(new_mtrx, tmp_mtrx)
      } else {
        type <- "Removed"
        ## In the case of removals, I need to only grep what is left after each iteration.
        if (semantic_column == "rownames") {
          idx <- grepl(pattern = string, x = rownames(new_annots))
        } else {
          idx <- grepl(pattern = string, x = new_annots[, semantic_column])
        }
        mesg("Hit ", sum(idx), " genes for term ", string, ".")
        idx <- ! idx
        ## So, we take the index of stuff to keep, and just subset on that index.
        new_annots <- new_annots[idx, ]
        new_mtrx <- new_mtrx[idx, ]
      }
    } ## End for loop
    end_rows <- nrow(new_mtrx)
    lost_rows <- start_rows - end_rows
    message("semantic_expt_filter(): Removed ", lost_rows, " genes.")
  } else {
    ## Instead of a string based sematic filter, take the topn most abundant
    medians <- rowMedians(mtrx)
    new_order <- order(medians, decreasing = TRUE)
    reordered <- mtrx[new_order, ]
    subset <- rownames(head(reordered, n = topn))
    new_annots <- annots[subset, ]
  }

  expressionset <- input[["expressionset"]]
  keepers <- rownames(new_annots)
  new_expressionset <- expressionset[keepers, ]
  new_libsizes <- colSums(exprs(new_expressionset))
  input[["expressionset"]] <- new_expressionset
  input[["libsize"]] <- new_libsizes
  return(input)
}

#' Extract a subset of samples following some rule(s) from an
#' experiment class.
#'
#' Sometimes an experiment has too many parts to work with conveniently, this
#' operation allows one to break it into smaller pieces.
#'
#' @param expt Expt chosen to extract a subset of data.
#' @param subset Valid R expression which defines a subset of the design to keep.
#' @param nonzero Look for a minimal number of nonzero genes.
#' @param ids List of sample IDs to extract.
#' @param coverage Request a minimum coverage/sample rather than text-based subset.
#' @param print_excluded Print out the samples which are removed via this filter?
#' @return metadata Expt class which contains the smaller set of data.
#' @seealso [Biobase] [pData()] [exprs()] [fData()]
#' @export
subset_expt <- function(expt, subset = NULL, ids = NULL,
                        nonzero = NULL, coverage = NULL,
                        print_excluded = FALSE) {
  starting_expressionset <- NULL
  starting_metadata <- NULL
  starting_samples <- sampleNames(expt)

  if ("ExpressionSet" %in% class(expt)) {
    starting_expressionset <- expt
    starting_metadata <- pData(starting_expressionset)
  } else if ("expt" %in% class(expt)) {
    starting_expressionset <- expt[["expressionset"]]
    starting_metadata <- pData(expt)
  } else {
    stop("expt is neither an expt nor ExpressionSet")
  }

  if (class(subset)[1] == "logical") {
    ids <- rownames(pData(expt))[subset]
    subset <- NULL
  }

  if (is.null(starting_metadata[["sampleid"]])) {
    starting_metadata[["sampleid"]] <- rownames(starting_metadata)
  }

  excluded <- c()
  if (!is.null(ids)) {
    if (is.numeric(ids)) {
      ids <- rownames(starting_metadata)[ids]
    } else if (is.logical(ids)) {
      ids <- rownames(starting_metadata)[ids]
    }

    string <- ""
    for (id in ids) {
      string <- glue("{string}|sampleid=='{id}'")
    }
    ## Remove the leading |
    subset <- substring(string, 2)
  }
  note_appended <- NULL
  subset_design <- NULL
  if (is.null(coverage) && is.null(nonzero)) {
    if (is.null(subset)) {
      subset_design <- starting_metadata
    } else {
      mesg("Using a subset expression.")
      r_expression <- paste("subset(starting_metadata,", subset, ")")
      subset_design <- eval(parse(text = r_expression))
      note_appended <- glue("Subsetted with {subset} on {date()}.
")
    }
    if (nrow(subset_design) == 0) {
      stop("When the subset was taken, the resulting design has 0 members.")
    }
    subset_design <- as.data.frame(subset_design, stringsAsFactors = FALSE)
    excluded_idx <- ! rownames(starting_metadata) %in% rownames(subset_design)
    excluded <- starting_samples[excluded_idx]
    if (isTRUE(print_excluded)) {
      message("The samples excluded are: ", toString(excluded), ".")
    }
  } else if (is.null(nonzero)) {
    ## If coverage is defined, then use it to subset based on the minimal desired coverage
    ## Perhaps in a minute I will make this work for strings like '1z' to get the lowest
    ## standard deviation or somesuch...
    mesg("Subsetting given a minimal number of counts/sample.")
    coverages <- colSums(exprs(expt))

    if (is.null(pData(expt)[["sample_coverage"]])) {
      pData(expt)[["sample_coverage"]] <- coverages
    }
    subset_idx <- coverages >= as.numeric(coverage) ## In case I quote it on accident.
    subset_design <- starting_metadata[subset_idx, ]
    subset_design <- as.data.frame(subset_design, stringsAsFactors = FALSE)
    if (isTRUE(print_excluded)) {
      message("The samples removed (and read coverage) when filtering samples with less than ",
              coverage, " reads are: ")
    }
    print(colSums(exprs(expt))[!subset_idx])
  } else if (is.null(coverage)) {
    ## Remove samples with less than this number of non-zero genes.
    mesg("Subsetting given a minimal number genes observed/sample.")
    nonzero_idx <- exprs(expt) != 0
    num_nonzero <- colSums(nonzero_idx)
    if (is.null(pData(expt)[["num_nonzero"]])) {
      pData(expt)[["num_nonzero"]] <- num_nonzero
    }
    remove_idx <- num_nonzero < nonzero
    if (sum(remove_idx) == 0) {
      message("No samples have fewer than ", nonzero, " observed genes.")
      return(expt)
    }
    samples_dropped <- num_nonzero[remove_idx]
    subset_design <- starting_metadata[!remove_idx, ]
    subset_design <- as.data.frame(subset_design, stringsAsFactors = FALSE)
    message("The samples (and read coverage) removed when filtering ",
            nonzero, " non-zero genes are: ")
    dropped <- exprs(expt)[, remove_idx]
    dropped_names <- names(num_nonzero[remove_idx])
    num_genes <- num_nonzero[!remove_idx]
    if (class(dropped)[1] == "numeric") {
      print(num_nonzero[remove_idx])
    } else if (class(dropped)[1] == "matrix") {
      print(colSums(dropped))
    }
    message("by number of genes.")
  } else {
    stop("Unable to determine what is being subset.")
  }
  ## This is to get around stupidity with respect to needing all factors to be
  ## in a DESeqDataSet
  starting_ids <- rownames(starting_metadata)
  subset_ids <- rownames(subset_design)
  subset_positions <- starting_ids %in% subset_ids
  starting_colors <- expt[["colors"]]
  subset_colors <- starting_colors[subset_positions, drop = TRUE]
  starting_conditions <- expt[["conditions"]]
  subset_conditions <- starting_conditions[subset_positions, drop = TRUE]
  starting_batches <- expt[["batches"]]
  subset_batches <- starting_batches[subset_positions, drop = TRUE]
  current_libsize <- expt[["libsize"]]
  subset_current_libsize <- current_libsize[subset_positions, drop = TRUE]
  subset_expressionset <- starting_expressionset[, subset_positions]
  notes <- expt[["notes"]]
  if (!is.null(note_appended)) {
    notes <- glue("{notes}{note_appended}")
  }
  for (col in seq_len(ncol(subset_design))) {
    if (class(subset_design[[col]])[1] == "factor") {
      subset_design[[col]] <- droplevels(subset_design[[col]])
    }
  }
  pData(subset_expressionset) <- subset_design
  ## Ensure that the tximport information is maintained!

  new_expt <- list(
    "title" = expt[["title"]],
    "notes" = toString(notes),
    "initial_metadata" = subset_design,
    "expressionset" = subset_expressionset,
    "design" = subset_design,
    "conditions" = subset_conditions,
    "batches" = subset_batches,
    "samplenames" = subset_ids,
    "colors" = subset_colors,
    "state" = expt[["state"]],
    "libsize" = subset_current_libsize)

  if (!is.null(expt[["tximport"]])) {
    raw_counts <- expt[["tximport"]][["raw"]][["counts"]]
    raw_abundance <- expt[["tximport"]][["raw"]][["abundance"]]
    raw_length <- expt[["tximport"]][["raw"]][["length"]]
    raw_subset_counts <- raw_counts[, subset_positions]
    raw_subset_abundance <- raw_abundance[, subset_positions]
    raw_subset_length <- raw_length[, subset_positions]
    scaled_counts <- expt[["tximport"]][["scaled"]][["counts"]]
    scaled_abundance <- expt[["tximport"]][["scaled"]][["abundance"]]
    scaled_length <- expt[["tximport"]][["scaled"]][["length"]]
    scaled_subset_counts <- scaled_counts[, subset_positions]
    scaled_subset_abundance <- scaled_abundance[, subset_positions]
    scaled_subset_length <- scaled_length[, subset_positions]
    new_tximport <- list(
      "raw" = list(
        "counts" = raw_subset_counts,
        "abundance" = raw_subset_abundance,
        "length" = raw_subset_length,
        "countsFromAbundance" = expt[["tximport"]][["raw"]][["countsFromAbundance"]]),
      "scaled" = list(
        "counts" = scaled_subset_counts,
        "abundance" = scaled_subset_abundance,
        "length" = scaled_subset_length,
        "countsFromAbundance" = expt[["tximport"]][["scaled"]][["countsFromAbundance"]]))
    new_expt[["tximport"]] <- new_tximport
  }

  class(new_expt) <- "expt"
  final_samples <- sampleNames(new_expt)
  message("subset_expt(): There were ", length(starting_samples), ", now there are ",
          length(final_samples), " samples.")
  return(new_expt)
}
setGeneric("subset_expt")

#' Subset a SummarizedExperiment with some extra syntax.
#'
#' @param expt Expt chosen to extract a subset of data.
#' @param subset Valid R expression which defines a subset of the design to keep.
#' @param nonzero Look for a minimal number of nonzero genes.
#' @param ids List of sample IDs to extract.
#' @param coverage Request a minimum coverage/sample rather than text-based subset.
#' @param print_excluded Print out the samples which are removed via this filter?
#' @return metadata Expt class which contains the smaller set of data.
#' @seealso [Biobase] [pData()] [exprs()] [fData()]
#' @export
setMethod(
  "subset_expt", signature = signature(expt = "SummarizedExperiment"),
  definition = function(expt, subset = NULL, ids = NULL,
                        nonzero = NULL, coverage = NULL,
                        print_excluded = TRUE) {
    subset_se(expt, subset = subset, ids = ids,
              nonzero = nonzero, coverage = coverage,
              print_excluded = print_excluded)
  })

#' I want an easy way to sum counts in eupathdb-derived data sets.
#' These have a few things which should make this relatively easy.
#' Notably: The gene IDs look like: "exon_ID-1 exon_ID-2 exon_ID-3"
#' Therefore we should be able to quickly merge these.
#'
#' @param counts Matrix/df/dt of count data.
#' @return The same data type but with the exons summed.
sum_eupath_exon_counts <- function(counts) {
  rownames(counts) <- gsub(pattern = "^exon_", replacement = "", x = rownames(counts))
  rownames(counts) <- gsub(pattern = "\\-1$", replacement = "", x = rownames(counts))
  multi_exon_idx <- grep(pattern = "\\-\\d+$", x = rownames(counts))
  for (idx in multi_exon_idx) {
    gene <- gsub(pattern = "\\-\\d+$", replacement = "", x = rownames(counts)[idx])
    exon_counts <- counts[idx, ]
    gene_counts <- counts[gene, ]
    total_counts <- gene_counts + exon_counts
    counts[gene, ] <- total_counts
  }
  counts <- counts[-multi_exon_idx, ]
  counts <- as.matrix(counts)
  return(counts)
}

#' Add some gene annotations based on the mean/variance in the data.
#'
#' Why?  Maria Adelaida is interested in pulling the least-variant
#' genes in our data, this seems like it might be generally
#' applicable.  Note, I made this slightly redundant by doing a cpm on
#' the data; as a result the proportion and mean values are
#' effectively identical.
#'
#' @param expt Expressionset to which to add this information.
#' @param convert Use this conversion,
#' @param transform and transformation,
#' @param norm and normalization.
#' @return Slightly modified gene annotations including the mean/variance.
#' @export
variance_expt <- function(expt, convert = "cpm", transform = "raw", norm = "raw") {
  start <- normalize_expt(expt, convert = convert,
                          transform = transform, norm = norm)
  df <- exprs(start)
  na_idx <- is.na(df)
  if (sum(na_idx) > 0) {
    warning("There are ", sum(na_idx), " NAs in this data, removing them.")
    message("There are ", sum(na_idx), " NAs in this data, removing them.")
    df[na_idx] <- 0
  }
  raw_sum <- rowSums(exprs(expt))
  raw_total <- sum(raw_sum)
  sums <- rowSums(df)
  total <- sum(sums)
  vars <- matrixStats::rowVars(df)
  sds <- matrixStats::rowSds(df)
  meds <- matrixStats::rowMedians(df)
  iqrs <- matrixStats::rowIQRs(df)
  mean <- rowMeans(df)
  fData(expt)[["exprs_gene_prop"]] <- sums / total
  fData(expt)[["exprs_gene_rawprop"]] <- raw_sum / raw_total
  fData(expt)[["exprs_gene_variance"]] <- vars
  fData(expt)[["exprs_gene_stdev"]] <- sds
  fData(expt)[["exprs_gene_mean"]] <- mean
  fData(expt)[["exprs_gene_median"]] <- meds
  fData(expt)[["exprs_gene_iqrs"]] <- iqrs
  fData(expt)[["exprs_cv"]] <- sds / mean
  return(expt)
}

#' An expt is an ExpressionSet superclass with a shorter name.
#'
#' It is also a simple list so that one may summarize it more simply,
#' provides colors and some slots to make one's life easier.
#' It is created via the function create_expt() which perhaps should be changed.
#'
#' Another important caveat: expressionSets and their methods are all S4; but I
#' did not want to write S4 methods, so I made my expt a S3 class.  As a result,
#' in order to make use of exprs, notes, pData, fData, and friends, I made use
#' of setMethod() to set up calls for the expressionSet portion of the expt
#' objects.
#'
#' @param ... Parameters for create_expt()
#' @slot title Title for the expressionSet.
#' @slot notes Notes for the expressionSet (redundant with S4 notes()).
#' @slot design Copy of the experimental metadata (redundant with pData()).
#' @slot annotation Gene annotations (redundant with fData()).
#' @slot gff_file filename of a gff file which feeds this data.
#' @slot state What is the state of the data vis a vis normalization,
#'  conversion, etc.
#' @slot conditions Usually the condition column from pData.
#' @slot batches Usually the batch column from pData.
#' @slot libsize Library sizes of the data in its current state.
#' @slot colors Chosen colors for plotting the data.
#' @slot tximport Data provided by tximport() to create the exprs() data.
#' @export expt
expt <- function(...) {
  create_expt(...)
}
setClass("expt")

#' @export
write_expt <- function(...) {
  write_se(...)
}

## EOF
