#' @include 01_hpgltools.R
NULL

#' Figure out when mappings were performed by their timestamp
#'
#' I got bit in the butt by mismatching ensembl IDs from some older
#' count tables and newer annotations.  Happily my biomart annotation
#' gatherer is smart enough to collect from the archive servers, so it
#' should not be difficult for me to ensure that they match in the
#' future.
#'
#' With that in mind, provide this function with the filename of some
#' metadata and the file column in it, and it will look at the first
#' file and return the year and month it was created.  Therefore, you
#' may ask ensembl for the appropriately dated gene annotations.
#'
#' @param metadata File containing the metadata for this experiment.
#'  If none is provided, this function will just give the current
#'  year, which is only what you want if this is brand new data.
#' @param column Sanitized column name in the metadata containing the
#'  count tables of interest.  If this is not provided, it will return
#'  the month/year of the timestamp for the metadata.  This has a
#'  reasonable chance of giving correct information.
#' @export
check_metadata_year <- function(metadata = NULL, column = NULL) {
  retlist <- list(
    "year" = format(Sys.time(), "%Y"),
    "month" = format(Sys.time(), "%m"))
  if (is.null(metadata)) {
    message("No metadata provided, assuming now is sufficient.")
  } else if (is.null(column)) {
    message("No file column was requested, assuming the metadata creation time is sufficient.")
    retlist[["year"]] <- format(info[["ctime"]], "%Y")
    retlist[["month"]] <- format(info[["ctime"]], "%m")
  } else {
    message("Checking the creation time on the first count table.")
    meta <- extract_metadata(metadata)
    files <- meta[[column]]
    first_file <- files[1]
    info <- file.info(first_file)
    retlist[["year"]] <- format(info[["ctime"]], "%Y")
    retlist[["month"]] <- format(info[["ctime"]], "%m")
  }
  return(retlist)
}

#' A silly function to guesstimate factor columns in metadata.
#'
#' Use the heuristic that any column with a number of different
#' elements which is <= (samples / ratio) has a reasonable chance of
#' being usable as a categorical.
#'
#' @param meta_df Input metadata.
#' @param ratio Heuristic of levels/ratio as the maximum number of
#'  allowed levels in the new factor.  If the putative number of
#'  levels is more than this, then assume it is not a factor.
guess_factors <- function(meta_df, ratio = 3) {
  guesses <- c()
  num_samples <- nrow(meta_df)
  max_levels <- num_samples / ratio
  for (f in seq_len(ncol(meta_df))) {
    name <- colnames(meta_df)[f]
    num_levels <- length(levels(as.factor(meta_df[[f]])))
    if (num_levels <= max_levels) {
      guesses <- c(guesses, name)
    }
  }
  return(guesses)
}

#' Pull metadata from a table (xlsx/xls/csv/whatever)
#'
#' I find that when I acquire metadata from a paper or collaborator, annoyingly
#' often there are many special characters or other shenanigans in the column
#' names.  This function performs some simple sanitizations.  In addition, if I
#' give it a filename it calls my generic 'read_metadata()' function before
#' sanitizing.
#'
#' @param metadata file or df of metadata
#' @param id_column Column in the metadat containing the sample names.
#' @param fill Fill missing data with this.
#' @param fill_condition Add a condition column if there is not one?
#' @param fill_batch Add a batch column if there is not one?
#' @param keep_underscore Sanitize underscores?
#' @param sanitize Perform my various sanitizers on the data?
#' @param condition_column conditions metadata column to be protected.
#' @param batch_column batch metadata column to be protected.
#' @param ... Arguments to pass to the child functions (read_csv etc).
#' @return Metadata dataframe hopefully cleaned up to not be obnoxious.
#' @examples
#'  \dontrun{
#'   sanitized <- extract_metadata("some_random_supplemental.xls")
#'   saniclean <- extract_metadata(some_goofy_df)
#' }
#' @export
extract_metadata <- function(metadata, id_column = "sampleid", fill = NULL,
                             fill_condition = TRUE, fill_batch = TRUE,
                             keep_underscore = TRUE, sanitize = TRUE,
                             condition_column = NULL, batch_column = NULL, ...) {
  ## FIXME: Now that this has been yanked into its own function,
  ## Make sure it sets good, standard rownames.
  file <- NULL

  meta_dataframe <- NULL
  meta_file <- NULL
  ## Handle this via dispatch
  if ("character" %in% class(metadata)) {
    ## This is a filename containing the metadata
    meta_file <- metadata
  } else if ("data.frame" %in% class(metadata) ||
               "DFrame" %in% class(metadata) ||
                 "data.table" %in% class(metadata)) {
    ## A data frame of metadata was passed.
    meta_dataframe <- metadata
  } else {
    stop("This requires either a file or meta data.frame.")
  }

  ## The two primary inputs for metadata are a csv/xlsx file or a dataframe, check for them here.
  if (is.null(meta_dataframe) && is.null(meta_file)) {
    stop("This requires either a csv file or dataframe of metadata describing the samples.")
  } else if (is.null(meta_file)) {
    ## punctuation is the devil
    sample_definitions <- meta_dataframe
  }  else {
    sample_definitions <- read_metadata(meta_file, fill_condition = fill_condition,
                                        fill_batch = fill_batch,
                                        ...)
    ## sample_definitions <- read_metadata(meta_file)
  }

  sample_definitions <- as.data.frame(sample_definitions)
  if (isTRUE(sanitize)) {
    if (isTRUE(keep_underscore)) {
      colnames(sample_definitions) <- gsub(pattern = "[^_[:^punct:]]", replacement = "",
                                           x = colnames(sample_definitions), perl = TRUE)
    } else {
      colnames(sample_definitions) <- gsub(pattern = "[[:punct:]]",
                                           replacement = "",
                                           x = colnames(sample_definitions))
    }
  }

  ## Try to ensure that we have a useful ID column by:
  ## 1. Look for data in the id_column column.
  ##  a.  If it is null, look at the rownames
  ##    i.  If they are 1...n, arbitrarily grab the first column.
  ##    ii. If not, use the rownames.
  ## Get appropriate row and column names.
  ## NOTE: the use of identical() vs. numbers fails when a sample is
  ## removed when reading the metadata, because the samples will have a break.
  current_rownames <- rownames(sample_definitions)
  numeric_test <- suppressWarnings(as.numeric(current_rownames))
  non_numeric <- sum(is.na(numeric_test))  ## If this is > 0, then these are not just numbers.
  if (is.null(sample_definitions[[id_column]])) {
    message("Did not find the column: ", id_column, ".")
    if (non_numeric > 0) {
      message("The rownames do not appear numeric, using them.")
      sample_definitions[[id_column]] <- rownames(sample_definitions)
    } else {
      message("Setting the ID column to the first column.")
      id_column <- colnames(sample_definitions)[1]
      rownames(sample_definitions) <- make.names(sample_definitions[[id_column]], unique = TRUE)
    }
  } else {
    ## 202311: I am not completely certain this logic change is what I want.
    rownames(sample_definitions) <- make.names(sample_definitions[[id_column]], unique = TRUE)
  }

  ## Drop empty rows in the sample sheet
  empty_samples <- which(sample_definitions[, id_column] == "" |
                           grepl(x = sample_definitions[, id_column], pattern = "^undef") |
                           is.na(sample_definitions[, id_column]) |
                           grepl(pattern = "^#", x = sample_definitions[, id_column]))
  if (length(empty_samples) > 0) {
    message("Dropped ", length(empty_samples),
            " rows from the sample metadata because the sample ID is blank.")
    sample_definitions <- sample_definitions[-empty_samples, ]
  }

  ## Drop duplicated elements.
  num_duplicated <- sum(duplicated(sample_definitions[[id_column]]))
  if (num_duplicated > 0) {
    message("There are ", num_duplicated,
            " duplicate rows in the sample ID column.")
    sample_definitions[[id_column]] <- make.names(sample_definitions[[id_column]],
                                                  unique = TRUE)
  }

  ## Now we should have consistent sample IDs, set the rownames.
  rownames(sample_definitions) <- sample_definitions[[id_column]]
  ## Check that condition and batch have been filled in.
  sample_columns <- colnames(sample_definitions)

  ## The various proteomics data I am looking at annoyingly starts with a number
  ## So make.names() prefixes it with X which is ok as far as it goes, but
  ## since it is a 's'amplename, I prefer an 's'.
  if (isTRUE(sanitize)) {
    rownames(sample_definitions) <- gsub(pattern = "^X([[:digit:]])",
                                         replacement = "s\\1",
                                         x = rownames(sample_definitions))
  }

  sample_columns_to_remove <- NULL
  for (col in seq_along(colnames(sample_definitions))) {
    sum_na <- sum(is.na(sample_definitions[[col]]))
    sum_null <- sum(is.null(sample_definitions[[col]]))
    sum_empty <- sum_na + sum_null
    if (sum_empty ==  nrow(sample_definitions)) {
      ## This column is empty.
      sample_columns_to_remove <- append(sample_columns_to_remove, col)
    }
  }
  if (length(sample_columns_to_remove) > 0) {
    sample_definitions <- sample_definitions[-sample_columns_to_remove]
  }

  sample_columns <- colnames(sample_definitions)
  if (!is.null(condition_column) && condition_column %in% sample_columns) {
    found_condition <- "condition" %in% sample_columns
    if (isTRUE(found_condition)) {
      message("There is already a condition column, backing it up to condition_bak")
      sample_definitions[["condition_backup"]] <- sample_definitions[["condition"]]
    }
    sample_definitions[["condition"]] <- sample_definitions[[condition_column]]
  }
  if (!is.null(batch_column) && batch_column %in% sample_columns) {
    found_batch <- "batch" %in% sample_columns
    if (!isTRUE(found_batch)) {
      message("There is already a batch column, backing it up to batch_bak")
      sample_definitions[["batch_backup"]] <- sample_definitions[["batch"]]
    }
    sample_definitions[["batch"]] <- sample_definitions[[batch_column]]
  }

  ## Now check for columns named condition and batch
  if (isTRUE(fill_condition)) {
    found_condition <- "condition" %in% sample_columns
    if (!isTRUE(found_condition)) {
      message("Did not find the condition column in the sample sheet.")
      message("Filling it in as undefined.")
      sample_definitions[["condition"]] <- "undefined"
    }
  }

  if (isTRUE(fill_batch)) {
    found_batch <- "batch" %in% sample_columns
    if (!isTRUE(found_batch)) {
      message("Did not find the batch column in the sample sheet.")
      message("Filling it in as undefined.")
      sample_definitions[["batch"]] <- "undefined"
    }
  }

  message("Checking the state of the condition column.")
  ## Check the condition/batch columns (if they exist) for NA values.
  if (!is.null(sample_definitions[["condition"]])) {
    ## Adding this casting of the condition/batch columns
    ## in case they are something weird like 'Date'
    sample_definitions[["condition"]] <- as.character(sample_definitions[["condition"]])
    na_idx <- is.na(sample_definitions[["condition"]])
    if (sum(na_idx) > 0) {
      warning("There were NA values in the condition column, setting them to 'undefined'.")
      sample_definitions[na_idx, "condition"] <- "undefined"
    }
    sample_definitions[["condition"]] <- as.factor(sample_definitions[["condition"]])
  }

  message("Checking the state of the batch column.")
  if (!is.null(sample_definitions[["batch"]])) {
    sample_definitions[["batch"]] <- as.character(sample_definitions[["batch"]])
    na_idx <- is.na(sample_definitions[["batch"]])
    if (sum(na_idx) > 0) {
      warning("There were NA values in the condition column, setting them to 'undefined'.")
      sample_definitions[na_idx, "batch"] <- "undefined"
    }
    sample_definitions[["batch"]] <- as.factor(sample_definitions[["batch"]])
  }

  message("Checking the condition factor.")
  ## Extract out the condition names as a factor
  condition_names <- unique(sample_definitions[["condition"]])
  if (is.null(condition_names)) {
    warning("There is no 'condition' field in the definitions, this will make many
analyses more difficult/impossible.")
  }
  ## Condition and Batch are not allowed to be numeric, so if they are just numbers,
  ## prefix them with 'c' and 'b' respectively.
  pre_condition <- unique(sample_definitions[["condition"]])
  pre_batch <- unique(sample_definitions[["batch"]])
  sample_definitions[["condition"]] <- gsub(pattern = "^(\\d+)$", replacement = "c\\1",
                                            x = sample_definitions[["condition"]])
  sample_definitions[["batch"]] <- gsub(pattern = "^(\\d+)$", replacement = "b\\1",
                                        x = sample_definitions[["batch"]])
  sample_definitions[["condition"]] <- factor(sample_definitions[["condition"]],
                                              levels = unique(sample_definitions[["condition"]]),
                                              labels = pre_condition)
  sample_definitions[["batch"]] <- factor(sample_definitions[["batch"]],
                                          levels = unique(sample_definitions[["batch"]]),
                                          labels = pre_batch)

  if (!is.null(fill)) {
    na_idx <- is.na(sample_definitions)
    sample_definitions[na_idx] <- fill
  }

  return(sample_definitions)
}

#' Automagically fill in a sample sheet with the results of the
#' various preprocessing tools.
#'
#' I am hoping to fill this little function out with a bunch of useful
#' file specifications and regular expressions.  If I do a good job,
#' then it should become trivial to fill in a sample sheet with lots
#' of fun useful numbers in preparations for creating a nice table
#' S1.  I am thinking to split this up into sections for
#' trimming/mapping/etc.  But for the moment I just want to add some
#' specifications/regexes and see if it proves itself robust.  If
#' Theresa reads this, I think this is another good candidate for a
#' true OO implmentation.  E.g. make a base-class for the metadata and
#' use S4 multi-dispatch to pick up different log files.  I wrote the
#' downstream functions with this in mind already, but I am too
#' stupid/lazy to do the full implementation until I am confident that
#' these functions/ideas actually have merit.
#'
#' @param starting_metadata Existing sample sheet or NULL.  When NULL
#'  it will look in basedir for subdirectories not named 'test' and
#'  ontaining subdirectories named 'scripts' and use them to create an
#'  empty sample sheet.
#' @param specification List containing one element for each new
#'  column to append to the sample sheet.  Each element in turn is a
#'  list containing column names and/or input filenames (and
#'  presumably other stuff as I think of it).
#' @param basedir Root directory containing the files/logs of metadata.
#' @param new_metadata Filename to which to write the new metadata
#' @param species Define a desired species when file hunting.
#' @param type Define a feature type when file hunting.
#' @param subtype More specific type in the file spec.
#' @param tag Filename tag.
#' @param verbose Currently just used to debug the regexes.
#' @param id_column Column number or name containing the gene IDs.
#' @param r1_input_column Column containing the required names.
#' @param r2_input_column Ibid.
#' @param md5 Seek out and/or calculate md5 checksums for this column?
#' @param ... This is one of the few instances where I used
#' ... intelligently.  Pass extra variables to the file specification
#' and glue will pick them up (note the species entries in the
#' example specifications.
#' @return For the moment it just returns the modified metadata, I
#'  suspect there is something more useful it should do.
#' @export
gather_preprocessing_metadata <- function(starting_metadata = NULL, specification = NULL,
                                          basedir = "preprocessing", new_metadata = NULL,
                                          species = "*", type = "genome", subtype = "gene",
                                          tag = "ID", verbose = FALSE, id_column = 1,
                                          r1_input_column = NULL, r2_input_column = NULL,
                                          md5 = FALSE, ...) {
  ## I want to create a set of specifications for different tasks:
  ## tnseq/rnaseq/assembly/phage/phylogenetics/etc.
  ## For the moment it assumes phage assembly.
  if (is.null(specification)) {
    specification <- make_rnaseq_spec()
  } else if (class(specification)[1] == "list") {
    mesg("Using provided specification")
  } else if (specification == "rnaseq") {
    specification <- make_rnaseq_spec()
  } else if (specification == "dnaseq") {
    specification <- make_dnaseq_spec()
  } else if (specification == "phage") {
    specification <- make_assembly_spec()
  }
  if (isTRUE(md5)) {
    specification[["md5_r1_raw"]] <- list(
      "column" = "input_r1")
    specification[["md5_r2_raw"]] <- list(
      "column" = "input_r2")
#    specification[["md5_r1_trimmed"]] <- list(
#      "column" = "trimmed_r1")
#    specification[["md5_r2_trimmed"]] <- list(
#      "column" = "trimmed_r2")
    if (!is.null(r1_input_column)) {
      specification[["md5_r1_raw"]] <- list(
        "column" = r1_input_column)
    }
    if (!is.null(r2_input_column)) {
      specification[["md5_r2_raw"]] <- list(
        "column" = r2_input_column)
    }
  }

  meta <- NULL
  if (is.null(starting_metadata)) {
    sample_directories <- list.dirs(path = basedir, recursive = FALSE)
    sample_names <- basename(sample_directories)
    meta <- data.frame(sampleid = sample_names)
    rownames(meta) <- meta[["sampleid"]]
    meta[["condition"]] <- "undefined"
    meta[["batch"]] <- "undefined"
    include_idx <- !grepl(pattern = "^test$", x = sample_names)
    meta <- meta[include_idx, ]
    script_test <- list.dirs(path = file.path(basedir, rownames(meta)))
    script_idx <- grepl(pattern = "scripts$", x = script_test)
    include_idx <- sort(basename(dirname(script_test[script_idx])))
    meta <- meta[include_idx, ]
    starting_metadata <- "sample_sheets/autodetected_samples.xlsx"
  } else if (class(starting_metadata)[1] == "data.frame") {
    meta <- starting_metadata
  } else {
    meta <- extract_metadata(starting_metadata, id_column = id_column,
                             ...)
  }

  if (class(starting_metadata)[1] == "data.frame" && is.null(new_metadata)) {
    new_metadata <- NULL
  } else if (is.null(new_metadata)) {
    new_metadata <- gsub(x = starting_metadata, pattern = "\\.xlsx$",
                         replacement = "_modified.xlsx")
  }

  old_meta <- meta
  ## Perhaps use sanitize instead?
  if (!is.null(meta[["sampleid"]])) {
    meta[["sampleid_backup"]] <- meta[["sampleid"]]
  }
  meta[["sampleid"]] <- gsub(pattern = "\\s+", replacement = "", x = meta[[id_column]])
  new_columns <- c()

  for (entry in seq_along(specification)) {
    entry_type <- names(specification[entry])
    if (verbose) {
      message("Starting ", entry_type, ": ", entry, ".")
    }
    new_column <- entry_type
    if (!is.null(specification[[entry_type]][["column"]])) {
      new_column <- specification[[entry_type]][["column"]]
    }
    if (new_column %in% colnames(meta)) {
      warning("Column: ", new_column, " already exists, replacing it.")
    }
    input_file_spec <- specification[[entry_type]][["file"]]
    new_spec <- NULL
    if (!is.null(specification[[entry_type]][["file_new"]])) {
      new_spec <- specification[[entry_type]][["file_new"]]
    }
    if (verbose) {
      message("Checking input_file_spec: ", input_file_spec, ".")
      message("basedir is ", basedir, ", species is ", species, " and type is ", type, ".")
    }
    new_entries <- dispatch_metadata_extract(
      meta, entry_type, input_file_spec, specification,
      basedir = basedir, verbose = verbose, species = species, type = type,
      subtype = subtype, tag = tag, new_spec = new_spec) ##, ...)
    test_entries <- new_entries
    na_idx <- is.na(test_entries)
    na_entries <- sum(na_idx)
    test_entries <- test_entries[!na_idx]
    empty_entries <- sum(test_entries == "")
    zero_entries <- sum(test_entries == "0")
    uninformative_entries <- na_entries + empty_entries + zero_entries
    all_empty <- uninformative_entries == length(new_entries)
    if (is.null(new_entries) || isTRUE(all_empty)) {
      if (isTRUE(verbose)) {
        message("Not including new entries for: ", new_column, ", it is empty.")
      }
    } else if ("data.frame" %in% class(new_entries)) {
      for (sp in colnames(new_entries)) {
        new_column_name <- glue("{entry_type}_{sp}")
        new_columns <- c(new_columns, new_column_name)
        meta[[new_column_name]] <- new_entries[[sp]]
      }
    } else {
      new_columns <- c(new_columns, new_column)
      meta[[new_column]] <- new_entries
    }
  } ## End iterating over every metadatum

  is_numberp <- function(x) all(varhandle::check.numeric(as.character(x)))
  for (colname in colnames(meta)) {
    test_number <- is_numberp(meta[[colname]])
    if (isTRUE(test_number)) {
      meta[[colname]] <- as.numeric(meta[[colname]])
    }
  }

  written <- NULL
  if (is.null(new_metadata)) {
    message("Not writing new metadata.")
  } else {
    message("Writing new metadata to: ", new_metadata)
    written <- write_xlsx(data = meta, excel = new_metadata)
  }
  ret <- list(
    "xlsx_result" = written,
    "new_file" = new_metadata,
    "new_columns" = new_column,
    "old_meta" = old_meta,
    "new_meta" = meta)
  class(ret) <- "preprocessing_metadata"
  return(ret)
}

#' Print the result from gather_preprocessing_metadata().
#'
#' @param x List composed of the xlsx output file, new columns added
#'  to it, and copies of the metadata before/after modification.
#' @param ... Other args to match the generic.
#' @export
print.preprocessing_metadata <- function(x, ...) {
  message("The result of performing extract_preprocessing_metadata().")
  columns <- toString(colnames(x[["new_meta"]]))
  message("The columns now are: ", columns, ".")
  return(invisible(x))
}

#' This is basically just a switch and set of regexes for finding the
#' numbers of interest in the various log files.
#'
#' When I initially wrote this, it made sense to me to have it
#' separate from the top-level function.  I am not sure that is true
#' now, having slept on it.
#'
#' @param meta Starting metadata
#' @param entry_type String which defines the type of log entry to
#'  hunt down.  If the specification does not include a column, this
#'  will be used as the column name to write to the metadata.
#' @param input_file_spec Glue specification defining the log file for
#'  each sample to hunt down.
#' @param specification This is the reason I am thinking having this
#'  as a separate function might be stupid.  I added it to make it
#'  easier to calculate ratios of column_x/column_y; but it is a
#'  def-facto argument to either get rid of input_file_spec as an arg
#'  or to just get rid of this function.
#' @param basedir Root directory containing the files/logs of metadata.
#' @param verbose used for testing regexes.
#' @param new_spec New specification resulting from this invocation.
#' @param species Choose a specific species for which to search (for filenames generally).
#' @param type Set the type of file to search.
#' @param subtype Specific annotation type to seek.
#' @param tag GFF tag to use when seeking IDs.
#' @param ... passed to glue to add more variables to the file spec.
#' @return Vector of entries which will be used to populate the new
#'  column in the metadata.
dispatch_metadata_extract <- function(meta, entry_type, input_file_spec,
                                      specification, basedir = "preprocessing", verbose = FALSE,
                                      new_spec = NULL, species = "*", type = "genome",
                                      subtype = "gene", tag = "ID", ...) {
  switchret <- switch(
    entry_type,
    "aragorn_tRNAs" = {
      search <- "^\\d+ genes found"
      replace <- "^(\\d+) genes found"
      mesg("Searching aragorn outputs for tRNA-like genes.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, basedir = basedir,
        verbose = verbose, new_spec = new_spec, as = "numeric", ...)
    },
    "assembly_fasta_nt" = {
      mesg("Searching for assembly fasta files.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, new_spec = new_spec, basedir = basedir)
    },
    "assembly_genbank_annotated" = {
      mesg("Searching for assembly genbank files.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, new_spec = new_spec, basedir = basedir)
    },
    "assembly_genbank_stripped" = {
      mesg("Searching for stripped down assembly genbank files.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, new_spec = new_spec, basedir = basedir)
    },
    "assembly_cds_amino_acids" = {
      mesg("Searching for assembly amino acid fasta files.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, new_spec = new_spec, basedir = basedir)
    },
    "assembly_cds_nucleotides" = {
      mesg("Searching for assembly CDS fasta files.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, new_spec = new_spec, basedir = basedir)
    },
    "assembly_gff" = {
      mesg("Searching for assembly gff files.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, new_spec = new_spec, basedir = basedir)
    },
    "assembly_tsv" = {
      mesg("Searching for assembly tsv feature files.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, new_spec = new_spec, basedir = basedir)
    },
    "assembly_xls" = {
      mesg("Searching for assembly xlsx feature files.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, new_spec = new_spec, basedir = basedir)
    },
    "bedtools_coverage_file" = {
      mesg("Searching for files containing coverage from bedtools")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, new_spec = new_spec,
        species = species, basedir = basedir)
    },
    "bbmap_coverage_stats" = {
      mesg("Searching for coverage statistics from bbmap.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, new_spec = new_spec, basedir = basedir)
    },
    "bbmap_coverage_per_nt" = {
      mesg("Searching for coverage/nucleotide from bbmap.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, new_spec = new_spec, basedir = basedir)
    },
    "deduplication_stats" = {
      mesg("Searching for the GATK/picard deduplication stats file.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, new_spec = new_spec,
        species = species, type = type, basedir = basedir)
    },
    "deeptools_coverage" = {
      mesg("Searching for deeptools coverage values.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, new_spec = new_spec,
        species = species, type = type, basedir = basedir)
    },
    "fastqc_pct_gc" = {
      ## %GC     62
      search <- "^%GC	\\d+$"
      replace <- "^%GC	(\\d+)$"
      mesg("Searching for the percent GC content from fastqc.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, new_spec = new_spec,
        as = "numeric", basedir = basedir, ...)
    },
    "fastqc_most_overrepresented" = {
      ## Two lines after:
      ## >>Overrepresented sequences     fail
      ## #Sequence       Count   Percentage      Possible Source
      ## CTCCGCTATCGGTTCTACATGCTTAGCCAGCTCTACTGAGTTAACTCCGCGCCGCCCGA     68757   1.9183387064598885      No Hit
      mesg("Skipping fastqc_most_overrepresented.")
      ##entries <- dispatch_regex_search(meta, search, replace,
      ##                                 input_file_spec, verbose = verbose, basedir = basedir,
      ##                                 ...)
      entries <- NULL
    },
    "filtered_relative_coverage" = {
      ## >1 length=40747 coverage=1.00x circular=true
      search <- "^>\\d+ length=\\d+ coverage=.*x.*$"
      replace <- "^>\\d+ length=\\d+ coverage=(.*)x.*$"
      mesg("Searching for relative coverage observed by unicycler/spades.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, new_spec = new_spec,
        basedir = basedir, which = "all", ...)
    },
    "final_gc_content" = {
      mesg("Searching for the GC content of a finalized assembly.")
      entries <- dispatch_gc(
        meta, input_file_spec, basedir = basedir, verbose = verbose, new_spec = new_spec)
    },
    "gatk_unpaired" = {
      search <- "^.+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t0\\.\\d+\t\\d+$"
      replace <- "^.+\t(\\d+)\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t0\\.\\d+\t\\d+$"
      mesg("Searching the number of unpaired reads observed by GATK.")
      entries <- dispatch_regex_search(
        meta, search, replace, species = species, input_file_spec, new_spec = new_spec,
        verbose = verbose, basedir = basedir, which = "first", ...)
    },
    "gatk_paired" = {
      search <- "^.+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t0\\.\\d+\t\\d+$"
      replace <- "^.+\t\\d+\t(\\d+)\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t0\\.\\d+\t\\d+$"
      mesg("Searching for the number of paired reads observed by GATK.")
      entries <- dispatch_regex_search(
        meta, search, replace, species = species, input_file_spec, new_spec = new_spec,
        verbose = verbose, basedir = basedir, which = "first", ...)
    },
    "gatk_supplementary" = {
      search <- "^.+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t0\\.\\d+\t\\d+$"
      replace <- "^.+\t\\d+\t\\d+\t(\\d+)\t\\d+\t\\d+\t\\d+\t\\d+\t0\\.\\d+\t\\d+$"
      mesg("Searching for the number of supplementary reads observed by GATK.")
      entries <- dispatch_regex_search(
        meta, search, replace, species = species, input_file_spec, verbose = verbose,
        new_spec = new_spec, basedir = basedir, which = "first", ...)
    },
    "gatk_unmapped" = {
      search <- "^.+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t0\\.\\d+\t\\d+$"
      replace <- "^.+\t\\d+\t\\d+\t\\d+\t(\\d+)\t\\d+\t\\d+\t\\d+\t0\\.\\d+\t\\d+$"
      mesg("Searching for the number of unmapped reads by GATK.")
      entries <- dispatch_regex_search(
        meta, search, replace, species = species, input_file_spec, verbose = verbose,
        new_spec = new_spec, basedir = basedir, which = "first", ...)
    },
    "gatk_unpaired_duplicates" = {
      search <- "^.+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t0\\.\\d+\t\\d+$"
      replace <- "^.+\t\\d+\t\\d+\t\\d+\t\\d+\t(\\d+)\t\\d+\t\\d+\t0\\.\\d+\t\\d+$"
      mesg("Searching for the number of unpaired duplicate reads observed by GATK.")
      entries <- dispatch_regex_search(
        meta, search, replace, species = species, input_file_spec, verbose = verbose,
        new_spec = new_spec, basedir = basedir, which = "first", ...)
    },
    "gatk_paired_duplicates" = {
      search <- "^.+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t0\\.\\d+\t\\d+$"
      replace <- "^.+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t(\\d+)\t\\d+\t0\\.\\d+\t\\d+$"
      mesg("Searching for the number of paired duplicate reads observed by GATK.")
      entries <- dispatch_regex_search(
        meta, search, replace, species = species, input_file_spec, verbose = verbose,
        new_spec = new_spec, basedir = basedir, which = "first", ...)
    },
    "gatk_paired_opt_duplicates" = {
      mesg("Searching for the number of optical paired duplicate reads observed by GATK.")
      search <- "^.+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t0\\.\\d+\t\\d+$"
      replace <- "^.+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t(\\d+)\t0\\.\\d+\t\\d+$"
      entries <- dispatch_regex_search(
        meta, search, replace, species = species, input_file_spec, verbose = verbose,
        new_spec = new_spec, basedir = basedir, which = "first", ...)
    },
    "gatk_duplicate_pct" = {
      mesg("Searching for the percentage duplication observed by GATK.")
      search <- "^.+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t0\\.\\d+\t\\d+$"
      replace <- "^.+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t(0\\.\\d+)\t\\d+$"
      entries <- dispatch_regex_search(
        meta, search, replace, species = species, input_file_spec, verbose = verbose,
        new_spec = new_spec, basedir = basedir, which = "first", ...)
    },
    "gatk_libsize" = {
      mesg("Searching for the library size by GATK.")
      search <- "^.+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t0\\.\\d+\t\\d+$"
      replace <- "^.+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t\\d+\t0\\.\\d+\t(\\d+)$"
      entries <- dispatch_regex_search(
        meta, search, replace, species = species, input_file_spec, verbose = verbose,
        new_spec = new_spec, basedir = basedir, which = "first", ...)
    },
    "glimmer_positive_strand" = {
      search <- "\\s+\\+\\d+\\s+"
      mesg("Searching for the number of + strand entries observed by glimmer.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir, new_spec = new_spec)
    },
    "glimmer_negative_strand" = {
      search <- "\\s+-\\d+\\s+"
      mesg("Searching for the number of - strand entries observed by glimmer.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir, new_spec = new_spec)
    },
    "hisat_alignment" = {
      mesg("Searching for the hisat2 bam output.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species, type = type,
        subtype = subtype, tag = tag, basedir = basedir, new_spec = new_spec)
    },
    "hisat_count_table" = {
      mesg("Searching for the count tables from hisat->htseq.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species, type = type,
        subtype = subtype, tag = tag, basedir = basedir, new_spec = new_spec)
    },
    "hisat_genome_single_concordant" = {
      search <- "^\\s+\\d+ \\(.+\\) aligned concordantly exactly 1 time"
      replace <- "^\\s+(\\d+) \\(.+\\) aligned concordantly exactly 1 time"
      mesg("Searching for number of concordant single-mapped reads to CDS by hisat.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, basedir = basedir, as = "numeric",
        verbose = verbose, type = type, species = species, new_spec = new_spec, ...)
    },
    "hisat_genome_multi_concordant" = {
      search <- "^\\s+\\d+ \\(.+\\) aligned concordantly >1 times"
      replace <- "^\\s+(\\d+) \\(.+\\) aligned concordantly >1 times"
      mesg("Searching for number of concordant multi-mapped reads to RNA by hisat.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, basedir = basedir, as = "numeric",
        verbose = verbose, type = type, species = species, new_spec = new_spec, ...)
    },
    "hisat_genome_single_all" = {
      search <- "^\\s+\\d+ \\(.+\\) aligned exactly 1 time"
      replace <- "^\\s+(\\d+) \\(.+\\) aligned exactly 1 time"
      mesg("Searching for number of all reads mapped to CDS by hisat.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, as = "numeric",
        basedir = basedir, type = type, species = species, new_spec = new_spec,
        ...)
    },
    "hisat_genome_multi_all" = {
      search <- "^\\s+\\d+ \\(.+\\) aligned >1 times"
      replace <- "^\\s+(\\d+) \\(.+\\) aligned >1 times"
      mesg("Searching for number of all reads multi-mapped to CDS by hisat.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, as = "numeric",
        basedir = basedir, type = type, species = species, new_spec = new_spec, ...)
    },
    "hisat_genome_percent_log" = {
      search <- "^.+% overall alignment rate"
      replace <- "^(.+)% overall alignment rate"
      mesg("Searching for the percent mapped in the hisat log.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir,
        type = type, species = species, new_spec = new_spec, ...)
    },
    "hisat_genome_pct_vs_trimmed" = {
      numerator_column <- specification[["hisat_genome_single_concordant"]][["column"]]
      if (is.null(numerator_column)) {
        numerator_column <- "hisat_genome_single_concordant"
      }
      denominator_column <- specification[["trimomatic_output"]][["column"]]
      if (is.null(denominator_column)) {
        denominator_column <- "trimomatic_output"
      }
      mesg("Searching for percent reads mapped by hisat.")
      entries <- dispatch_metadata_ratio(meta, numerator_column, denominator_column)
    },
    "hisat_genome_input_reads" = {
      search <- "^[[:digit:]]+ reads; of these"
      replace <- "^([[:digit:]]+) reads.+$"
      mesg("Searching for number of concordant single-mapped reads to rRNA by hisat.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, basedir = basedir, as = "numeric",
        verbose = verbose, species = species, new_spec = new_spec,
        ...)
    },
    "hisat_observed_genes" = {
      search <- "^.*\\s+0$"
      ## In featureCounts, the column of interest is the 7th
      backup_search <- ""
      mesg("Searching for the number of genes observed by hisat.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir, inverse = TRUE,
        new_spec = new_spec, type = type, subtype = subtype, tag = tag, species = species)
    },
    "hisat_sum_genes" = {
      mesg("Searching for the number of reads counted to genes.")
      entries <- dispatch_sum_column(
        meta, input_file_spec, verbose = verbose, basedir = basedir, new_spec = new_spec,
        type = type, subtype = subtype, tag = tag, species = species)

    },
    "hisat_sum_genes_pct" = {
      numerator_column <- specification[["hisat_sum_genes"]][["column"]]
      if (is.null(numerator_column)) {
        numerator_column <- "hisat_sum_genes"
      }
      denominator_column <- specification[["trimomatic_output"]][["column"]]
      if (is.null(denominator_column)) {
        denominator_column <- "trimomatic_output"
      }
      mesg("Searching for percent reads counted by hisat.")
      entries <- dispatch_metadata_ratio(meta, numerator_column, denominator_column,
                                         species = species, as = "percent")

    },
    "hisat_observed_mean_exprs" = {
      mesg("Searching for the mean expression observed by hisat.")
      entries <- dispatch_csv_search(
        meta, "last", input_file_spec, file_type = "tsv", which = "function", chosen_func = "mean",
        verbose = verbose, basedir = basedir, type = type, subtype = subtype, tag = tag,
        species = species, as = "numeric")
    },
    "hisat_observed_median_exprs" = {
      mesg("Searching for the median expression observed by hisat.")
      entries <- dispatch_csv_search(
        meta, "last", input_file_spec, file_type = "tsv", which = "function", chosen_func = "median",
        verbose = verbose, basedir = basedir, type = type, subtype = subtype, tag = tag,
        species = species, as = "numeric")
    },
    "hisat_rrna_input_reads" = {
      search <- "^[[:digit:]]+ reads; of these"
      replace <- "^([[:digit:]]+) reads.+"
      mesg("Searching for number of concordant single-mapped reads to rRNA by hisat.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, basedir = basedir, as = "numeric",
        verbose = verbose, species = species, type = "rRNA", new_spec = new_spec,
        ...)
    },
    "hisat_rrna_single_concordant" = {
      search <- "^\\s+\\d+ \\(.+\\) aligned concordantly exactly 1 time"
      replace <- "^\\s+(\\d+) \\(.+\\) aligned concordantly exactly 1 time"
      mesg("Searching for number of concordant single-mapped reads to rRNA by hisat.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, basedir = basedir, as = "numeric",
        verbose = verbose, species = species, type = "rRNA", new_spec = new_spec, ...)
    },
    "hisat_rrna_multi_concordant" = {
      search <- "^\\s+\\d+ \\(.+\\) aligned concordantly >1 times"
      replace <- "^\\s+(\\d+) \\(.+\\) aligned concordantly >1 times"
      mesg("Searching for number of concordant multi-mapped reads to rRNA by hisat.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, basedir = basedir, as = "numeric",
        verbose = verbose, type = "rRNA", species = species, new_spec = new_spec, ...)
    },
    "hisat_rrna_pct_vs_trimmed" = {
      numerator_column <- specification[["hisat_rrna_multi_concordant"]][["column"]]
      if (is.null(numerator_column)) {
        numerator_column <- "hisat_rrna_multi_concordant"
      }
      numerator_add <- specification[["hisat_rrna_single_concordant"]][["column"]]
      if (is.null(numerator_add)) {
        numerator_add <- "hisat_rrna_single_concordant"
      }
      denominator_column <- specification[["trimomatic_output"]][["column"]]
      if (is.null(denominator_column)) {
        denominator_column <- "trimomatic_output"
      }
      mesg("Searching for percent rRNA mapped by hisat.")
      entries <- dispatch_metadata_ratio(
        meta, numerator_column, denominator_column, numerator_add = numerator_add)
    },
    "hisat_rrna_percent_log" = {
      search <- "^.+% overall alignment rate"
      replace <- "^(.+)% overall alignment rate"
      mesg("Searching for the percent mapped in the hisat log.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir,
        type = "rRNA", species = species, new_spec = new_spec, ...)
    },
    "hisat_unmapped" = {
      search <- "^\\s+\\d+ \\(.+\\) aligned 0 times"
      replace <- "^\\s+(\\d+) \\(.+\\) aligned 0 times"
      mesg("Searching for number of all reads that did not map to CDS by hisat.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, as = "numeric",
        basedir = basedir, type = type, species = species, new_spec = new_spec, ...)
    },
    "host_filter_species" = {
      search <- "^.*$"
      replace <- "^(.*)$"
      mesg("Searching for the most likely host species for reads according to kraken.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, which = "first",
        verbose = verbose, basedir = basedir)
    },
    "ictv_taxonomy" = {
      ## column <- "taxon"
      column <- "name"
      mesg("Searching for the likely ICTV taxonomy of an assembly.")
      entries <- dispatch_csv_search(
        meta, column, input_file_spec, file_type = "tsv", which = "all",
        verbose = verbose, basedir = basedir, ...)
    },
    "ictv_accession" = {
      column <- "hit_accession"
      mesg("Searching for the likely ICTV accession of an assembly.")
      entries <- dispatch_csv_search(
        meta, column, input_file_spec, file_type = "tsv", which = "all", verbose = verbose,
        basedir = basedir, ...)
    },
    "ictv_family" = {
      ## column <- "taxon"
      column <- "hit_family"
      mesg("Searching for the likely ICTV virus family of an assembly.")
      entries <- dispatch_csv_search(
        meta, column, input_file_spec, file_type = "tsv", which = "first",
        verbose = verbose, basedir = basedir, ...)
    },
    "ictv_genus" = {
      ## column <- "taxon"
      column <- "hit_genus"
      mesg("Searching for the likely ICTV virus genus of an assembly.")
      entries <- dispatch_csv_search(
        meta, column, input_file_spec, file_type = "tsv", which = "first",
        verbose = verbose, basedir = basedir, ...)
    },
    "input_r1" = {
      search <- "^\\s+<\\(\\w* .+\\).*$"
      replace <- "^\\s+<\\(\\w* (.+?)\\).*$"
      backup_search <- "^\\s*.*\\.fastq\\.gz\\s+.*\\.fastq\\.gz.*$"
      backup_replace <- "^\\s+(.*?\\.fastq\\.gz).*"
      mesg("Searching for the input R1.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir,
        backup_search = backup_search, backup_replace = backup_replace,
        ...)
    },
    "input_r2" = {
      search <- "^\\s+<\\(\\w* .+\\) <\\(\\w* .+\\).*$"
      replace <- "^\\s+<\\(\\w* .+\\) <\\(\\w* (.+)\\).*$"
      backup_search <- "^\\s*.*\\.fastq\\.gz\\s+.*\\.fastq\\.gz.*$"
      backup_replace <- "^\\s+.*?\\.fastq\\.gz\\s+(.*\\.fastq\\.gz).*$"
      mesg("Searching for the input R2.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir,
        backup_search = backup_search, backup_replace = backup_replace,
        ...)
    },
    "interpro_signalp_hits" = {
      search <- "\\tSignalP.*\\t"
      mesg("Searching for interpro signalP hits.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "interpro_phobius_hits" = {
      search <- "\\tPhobius\\t"
      mesg("Searching for interpro phobius hits.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "interpro_pfam_hits" = {
      search <- "\\tPfam\\t"
      mesg("Searching for interpro pfam hits.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "interpro_tmhmm_hits" = {
      search <- "\\tTMHMM\\t"
      mesg("Searching for interpro TMHMM hits.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "interpro_cdd_hits" = {
      search <- "\\tCDD\\t"
      mesg("Searching for interpro CDD hits.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "interpro_smart_hits" = {
      search <- "\\tSMART\\t"
      mesg("Searching for interpro SMART hits.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "interpro_gene3d_hits" = {
      search <- "\\tGene3D\\t"
      mesg("Searching for interpro gene3D hits.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "interpro_superfamily_hits" = {
      search <- "\\tSUPERFAMILY\\t"
      mesg("Searching for interpro superfamily hits.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "jellyfish_count_table" = {
      mesg("Searching for hits/kmer matrices by jellyfish.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "jellyfish_observed" = {
      search <- "^.*$"
      mesg("Searching for the number of separate kmers observed by jellyfish.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "kraken_viral_classified" = {
      search <- "^\\s+\\d+ sequences classified.*$"
      replace <- "^\\s+(\\d+) sequences classified.*$"
      mesg("Searching for reads classified by the kraken virus database.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir, ...)
    },
    "kraken_viral_unclassified" = {
      search <- "^\\s+\\d+ sequences unclassified.*$"
      replace <- "^\\s+(\\d+) sequences unclassified.*$"
      mesg("Searching for reads not classified by the kraken virus database.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir, ...)
    },
    "kraken_first_viral_species" = {
      search <- "^.*s__.*\\t\\d+$"
      replace <- "^.*s__(.*)\\t\\d+$"
      mesg("Searching for the most represented species by the kraken viral database.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir, ...)
    },
    "kraken_first_viral_species_reads" = {
      search <- "^.*s__.*\\t\\d+$"
      replace <- "^.*s__.*\\t(\\d+)$"
      mesg("Searching for the number of reads most represented by the kraken viral database.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir, ...)
    },
    "kraken_matrix" = {
      mesg("Searching for a matrix of reads/taxonomy from kraken.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "kraken_standard_classified" = {
      search <- "^\\s+\\d+ sequences classified.*$"
      replace <- "^\\s+(\\d+) sequences classified.*$"
      mesg("Searching for reads classified by the kraken standard database.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir, ...)
    },
    "kraken_standard_unclassified" = {
      search <- "^\\s+\\d+ sequences unclassified.*$"
      replace <- "^\\s+(\\d+) sequences unclassified.*$"
      mesg("Searching for reads not classified by the kraken standard database.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir, ...)
    },
    "kraken_first_standard_species" = {
      search <- "^.*s__.*\\t\\d+$"
      replace <- "^.*s__(.*)\\t\\d+$"
      mesg("Searching for the most represented species by the kraken standard database.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir, ...)
    },
    "kraken_first_standard_species_reads" = {
      search <- "^.*s__.*\\t\\d+$"
      replace <- "^.*s__.*\\t(\\d+)$"
      mesg("Searching for the number of reads most represented by the kraken standard database.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir, ...)
    },
    "md5_r1_raw" = {
      mesg("Gathering MD5 Sums for r1.")
      entries <- dispatch_md5(meta, entry_type, input_file_spec, specification,
                              basedir = basedir, species = species, type = type,
                              subtype = subtype, tag = tag, ...)
    },
    "md5_r2_raw" = {
      mesg("Gathering MD5 Sums for r2.")
      entries <- dispatch_md5(meta, entry_type, input_file_spec, specification,
                              basedir = basedir, species = species, type = type,
                              subtype = subtype, tag = tag, ...)
    },
    #"md5_r1_trimmed" = {
    #  mesg("Gathering MD5 Sums for r1 trimmed.")
    #  entries <- dispatch_md5(meta, entry_type, input_file_spec, specification,
    #                          basedir = basedir, species = species, type = type,
    #                          subtype = subtype, tag = tag, ...)
    #},
    #"md5_r2_trimmed" = {
    #  mesg("Gathering MD5 Sums for r2 trimmed.")
    #  entries <- dispatch_md5(meta, entry_type, input_file_spec, specification,
    #                          basedir = basedir, species = species, type = type,
    #                          subtype = subtype, tag = tag, ...)
    #},

    "notes" = {
      search <- "^.*$"
      replace <- "^(.*)$"
      mesg("Searching for the text in a notes.txt file.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir,
        which = "all")
    },
    "pernt_mean_coverage" = {
      column <- "Coverage"
      mesg("Searching for mean coverage/nucleotide from bbmap.")
      entries <- dispatch_csv_search(
        meta, column, input_file_spec, file_type = "tsv", which = "function",
        chosen_func = "mean", verbose = verbose, basedir = basedir, as = "numeric")
    },
    "pernt_median_coverage" = {
      column <- "Coverage"
      mesg("Searching for median coverage/nucleotide from bbmap.")
      entries <- dispatch_csv_search(
        meta, column, input_file_spec, file_type = "tsv", chosen_func = "median",
        basedir = basedir, which = "function", as = "numeric", verbose = verbose)
    },
    "pernt_max_coverage" = {
      column <- "Coverage"
      mesg("Searching for maximum coverage/nucleotide from bbmap.")
      entries <- dispatch_csv_search(
        meta, column, input_file_spec, file_type = "tsv", which = "function", chosen_func = "max",
        verbose = verbose, basedir = basedir)
    },
    "pernt_min_coverage" = {
      column <- "Coverage"
      mesg("Searching for minimum coverage/nucleotide from bbmap.")
      entries <- dispatch_csv_search(
        meta, column, input_file_spec, file_type = "tsv", which = "function", chosen_func = "min",
        verbose = verbose, basedir = basedir)
    },
    "phageterm_dtr_length" = {
      mesg("Searching for the DTR length according to phageterm.")
      entries <- dispatch_fasta_lengths(
        meta, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "phanotate_positive_strand" = {
      search <- "\\t\\+\\t"
      mesg("Searching for the number of features on the + strand observed by phanotate.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "phanotate_negative_strand" = {
      search <- "\\t-\\t"
      mesg("Searching for the number of features on the - strand observed by phanotate.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "phastaf_num_hits" = {
      search <- "^.*$"
      mesg("Searching for the number of potential hits observed by phastaf.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "prodigal_positive_strand" = {
      search <- "\\t\\+\\t"
      mesg("Searching for the number of features on the + strand observed by prodigal.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "prodigal_negative_strand" = {
      search <- "\\t-\\t"
      mesg("Searching for the number of features on the - strand observed by prodigal.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "possible_host_species" = {
      search <- ".*times\\.$"
      mesg("Searching for the species of a likely host according to kraken.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "racer_changed" = {
      search <- "^Number of changed positions"
      replace <- "^Number of changed positions\\s+(\\d+)$"
      mesg("Searching for the number of bases corrected by RACER.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir,
        which = "first", ...)
    },
    "salmon_stranded" = {
      search <- "Automatically detected most likely library type as .+$"
      replace <- "^.*Automatically detected most likely library type as (\\w+)$"
      mesg("Searching the likely library type observed by salmon.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose,
        basedir = basedir, which = "first", ...)
    },
    "salmon_mapped" = {
      search <- "Counted .+ total reads in the equivalence classes"
      replace <- "^.*Counted (.+)? total reads in the equivalence classes"
      mesg("Searching the number of reads quantified by salmon.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose,
        basedir = basedir, ...)
    },
    "salmon_percent" = {
      search <- "^.* Mapping rate = \\d+\\.\\d+%"
      replace <- "^.* Mapping rate = (\\d+\\.\\d+)%"
      mesg("Searching the percentage reads quantified by salmon.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose,
        basedir = basedir, ...)
    },
    "salmon_count_table" = {
      mesg("Searching for the salmon quantitation file.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species,
        type = "genome", basedir = basedir)
    },
    "salmon_observed_genes" = {
      search <- "^.*\t0\\.0+$"
      mesg("Searching for the number of genes observed by salmon.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, basedir = basedir, inverse = TRUE,
        type = type, species = species)
    },
    "shovill_contigs" = {
      ## [shovill] It contains 1 (min=131) contigs totalling 40874 bp.
      search <- "^\\[shovill\\] It contains \\d+ .*$"
      replace <- "^\\[shovill\\] It contains (\\d+) .*$"
      mesg("Searching for the number of contigs observed by shovill.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir,
        as = "numeric", which = "last", ...)
    },
    "shovill_length" = {
      ## [shovill] It contains 1 (min=131) contigs totalling 40874 bp.
      search <- "^\\[shovill\\] It contains \\d+ .* contigs totalling \\d+ bp\\.$"
      replace <- "^\\[shovill\\] It contains \\d+ .* contigs totalling (\\d+) bp\\.$"
      mesg("Searching for the assembly length observed by shovill.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir,
        which = "last", as = "numeric", ...)
    },
    "shovill_estlength" = {
      ## [shovill] Assembly is 109877, estimated genome size was 117464 (-6.46%)
      search <- "^\\[shovill\\] Assembly is \\d+\\, estimated genome size was \\d+.*$"
      replace <- "^\\[shovill\\] Assembly is \\d+\\, estimated genome size was (\\d+).*$"
      mesg("Searching for the estimated genome size observed by shovill.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir,
        which = "last", as = "numeric", ...)
    },
    "shovill_minlength" = {
      ## [shovill] It contains 1 (min=145) contigs totalling 109877 bp.
      search <- "^\\[shovill\\] It contains \\d+ \\(min=\\d+\\).*$"
      replace <- "^\\[shovill\\] It contains \\d+ \\(min=(\\d+)\\).*$"
      mesg("Searching for the minimum contig size observed by shovill.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir,
        which = "last", as = "numeric", ...)
    },
    "trimomatic_input" = {
      search <- "^Input (Read Pairs|Reads): \\d+ .*$"
      replace <- "^Input (Read Pairs|Reads): (\\d+) .*$"
      mesg("Searching for the number of trimomatic input reads.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir,
        extraction = "\\2", ...)
    },
    "trimomatic_output" = {
      search <- "^Input (Read Pairs|Reads): \\d+ (Both Surviving|Surviving): \\d+ .*$"
      replace <- "^Input (Read Pairs|Reads): \\d+ (Both Surviving|Surviving): (\\d+) .*$"
      mesg("Searching for the number of trimomatic output reads.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir,
        extraction = "\\3", ...)
    },
    "trimomatic_ratio" = {
      ## I think we can assume that the trimomatic ratio will come immediately after input/output
      numerator_column <- "trimomatic_output"
      if (!is.null(specification[["trimomatic_output"]][["column"]])) {
        numerator_column <- specification[["trimomatic_output"]][["column"]]
      }
      denominator_column <- "trimomatic_input"
      if (!is.null(specification[["trimomatic_input"]][["column"]])) {
        denominator_column <- specification[["trimomatic_input"]][["column"]]
      }
      mesg("Searching for the proportion of output/input from trimomatic.")
      entries <- dispatch_metadata_ratio(
        meta, numerator_column, denominator_column, verbose = verbose)
    },
    "trimmed_r1" = {
      search <- "^output_r1=.+$"
      replace <- "^output_r1=(.+)$"
      mesg("Searching for the trimmed R1.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir, ...)
    },
    "trimmed_r2" = {
      search <- "^output_r2=.+$"
      replace <- "^output_r2=(.+)$"
      mesg("Searching for the trimmed R2.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir, ...)
    },
    "tRNA_hits" = {
      ## >1 length=40747 depth=1.00x circular=true
      search <- "^.*Found \\d+ tRNAs"
      replace <- "^.*Found (\\d+) tRNAs"
      mesg("Searching for the number of putative tRNAs observed by prokka.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "umi_dedup_output_bam" = {
      mesg("Searching for umi-filtered alignments.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species,
        type = "genome", basedir = basedir)
    },
    "umi_dedup_output_count" = {
      mesg("Searching for umi-filtered counts.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species,
        type = "genome", basedir = basedir)
    },
    "umi_dedup_chimeric" = {
      ## Example line:
      ## 2024-06-07 20:26:12,771 INFO Reads: Input Reads: 18614709, Chimeric read pair: 1452937, Single end unmapped: 834185
      search <- "^.* Chimeric read pair: \\d+\\,.*$"
      replace <- "^.* Chimeric read pair: (\\d+)\\,.*$"
      mesg("Searching for the umi_tools chimeric reads.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "umi_dedup_num_reads_in" = {
      ## Example line:
      ## 2024-06-07 20:26:12,771 INFO Reads: Input Reads: 18614709, Chimeric read pair: 1452937, Single end unmapped: 834185
      search <- "^.* INFO Reads: Input Reads: \\d+\\,.*$"
      replace <- "^.* INFO Reads: Input Reads: (\\d+)\\,.*$"
      mesg("Searching for the umi_tools number of input reads.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "umi_dedup_num_reads_out" = {
      ## Example line:
      ## 2024-06-07 20:26:12,771 INFO Number of reads out: 8648625
      search <- "^.* Number of reads out: \\d+$"
      replace <- "^.* Number of reads out: (\\d+)$"
      mesg("Searching for the umi_tools number of output reads.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "umi_dedup_pct_reads" = {
      mesg("Getting the percent of deduplication reads output/input.")
      numerator_column <- "umi_dedup_num_reads_out"
      denominator_column <- "umi_dedup_num_reads_in"
      entries <- dispatch_metadata_ratio(meta, numerator_column, denominator_column)
    },
    "umi_dedup_deduplicated_positions" = {
      ## Example line:
      ## 2024-06-07 20:26:12,771 INFO Total number of positions deduplicated: 1753865
      search <- "^.* number of positions deduplicated: \\d+$"
      replace <- "^.* number of positions deduplicated: (\\d+)$"
      mesg("Searching for the umi_tools number of deduplicated positions.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "umi_dedup_mean_umi_per_pos" = {
      ## Example line:
      ## 2024-06-07 20:26:12,772 INFO Mean number of unique UMIs per position: 8.60
      search <- "^.* Mean number of unique UMIs per position: \\d+\\.*\\d*$"
      replace <- "^.* Mean number of unique UMIs per position: (\\d+\\.*\\d*)$"
      mesg("Searching for the umi_tools mean number of umis per position.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "umi_dedup_max_umi_per_pos" = {
      ## Example line:
      ## 2024-06-07 20:26:12,772 INFO Max. number of unique UMIs per position: 38924
      search <- "^.* Max\\. number of unique UMIs per position: \\d+$"
      replace <- "^.* Max\\. number of unique UMIs per position: (\\d+)$"
      mesg("Searching for the umi_tools maximum number of umis per position.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "umi_dedup_edit_distance" = {
      mesg("Searching for umi-dedup distance statistics.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species,
        type = "genome", basedir = basedir)
    },
    "umi_dedup_umi_positions" = {
      mesg("Searching for umi-dedup position statistics.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species,
        type = "genome", basedir = basedir)
    },
    "umi_dedup_umi_stats" = {
      mesg("Searching for umi-dedup statistics.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species,
        type = "genome", basedir = basedir)
    },
    "umi_extract_r1_output" = {
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species,
        type = "genome", basedir = basedir)
    },
    "umi_extract_r2_output" = {
      ## Example line:
      ## # read2_out                               : outputs/01umi_tools/r2_extracted.fastq.gz
      search <- "^# read2_out\\s+:.*$"
      replace <- "^# read2_out\\s+: (.+)$"
      mesg("Searching for the umi_tools R2 extract output.")
      entries <- dispatch_filename_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir)
    },
    "umi_extract_matches" = {
      ## Example line:
      ## 2024-06-06 11:27:37,727 INFO regex matches read2: 21454044
      search <- "^.+ regex matches read2: \\d+$"
      replace <- "^.+ regex matches read2: (\\d+)$"
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose,
        basedir = basedir, as = "numeric", ...)
    },
    "unicycler_lengths" = {
      ## >1 length=40747 depth=1.00x circular=true
      search <- "^>\\d+ length=\\d+ depth.*$"
      replace <- "^>\\d+ length=(\\d+) depth.*$"
      mesg("Searching for contig lengths observed by unicycler/spades.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir,
        which = "all", ...)
    },
    "unicycler_relative_coverage" = {
      ## >1 length=40747 depth=1.00x circular=true
      search <- "^>\\d+ length=\\d+ depth=.*x.*$"
      replace <- "^>\\d+ length=\\d+ depth=(.*)x.*$"
      mesg("Searching for relative coverage observed by unicycler.")
      entries <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose, basedir = basedir,
        which = "all", ...)
    },
    "freebayes_observed" = {
      search <- ".*"
      mesg("Searching for the number of variants observed by freebayes.")
      entries <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose, species = species, basedir = basedir)
    },
    "freebayes_observed_file" = {
      mesg("Searching for the matrix of all freebayes observations.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species, basedir = basedir)
    },
    "freebayes_variants_by_gene" = {
      mesg("Searching for variants/gene observed by freebayes.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species, basedir = basedir)
    },
    "freebayes_variants_table" = {
      mesg("REMOVE ME? Redundant with freebayes_observed_file, REMOVE ME?")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species, basedir = basedir)
    },
    "freebayes_modified_genome" = {
      mesg("Searching for the modified genome from freebayes.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species, basedir = basedir)
    },
    "freebayes_bcf_file" = {
      mesg("Searching for the freebayes bcf output.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species, basedir = basedir)
    },
    "freebayes_penetrance_file" = {
      mesg("Searching for a table of variant penetrance according to freebayes.")
      entries <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species, basedir = basedir)
    },
    {
      stop("I do not know this spec: ", entry_type)
    })
  if ("data.frame" %in% class(entries)) {
    return(entries)
  }
  if (!is.null(entries)) {
    entries <- gsub(pattern = ",", replacement = "", x = entries)
  }
  return(entries)
}

#' Count the number of lines in an input file spec and add it to the metadata.
#'
#' Sometimes the number of lines of a file is a good proxy for some
#' aspect of a sample. For example, jellyfish provides 1 line for
#' every kmer observed in a sample.  This function extracts that
#' number and puts it into each cell of a sample sheet.
#'
#' @param meta Input metadata
#' @param search Pattern to count
#' @param input_file_spec Input file specification to hunt down the
#'  file of interest.
#' @param verbose Print diagnostic information while running?
#' @param species Specify a species to search for.
#' @param basedir Root directory containing the files/logs of metadata.
#' @param type Add columns for only the genome mapping and/or rRNA by default.
#' @param subtype A specific annotation type to seek.
#' @param tag GFF tag to use for IDs.
#' @param new_spec New specification resulting from this function.
#' @param inverse Count the lines that do _not_ match the pattern.
dispatch_count_lines <- function(meta, search, input_file_spec, verbose = verbose,
                                 species = "*", basedir = "preprocessing",
                                 type = "genome", subtype = "gene", tag = "ID",
                                 new_spec = NULL, inverse = FALSE) {

  if (length(species) > 1) {
    output_entries <- data.frame(row.names = seq_len(nrow(meta)))
    for (s in seq_len(length(species))) {
      species_name <- species[s]
      output_entries[[species_name]] <- dispatch_count_lines(
        meta, search, input_file_spec, verbose = verbose,
        species = species_name, basedir = basedir, inverse = inverse)
    }
    return(output_entries)
  }

  filenames_with_wildcards <- glue(input_file_spec)
  if (isTRUE(verbose)) {
    message("Example count filename: ", filenames_with_wildcards[1], ".")
  }
  output_entries <- rep(0, length(filenames_with_wildcards))
  for (row in seq_len(nrow(meta))) {
    found <- 0
    ## Just in case there are multiple matches
    input_file <- Sys.glob(filenames_with_wildcards[row])[1]
    if (length(input_file) == 0) {
      mesg("There is no file matching: ", filenames_with_wildcards[row],
           ".")
      output_entries[row] <- ''
      next
    }
    if (is.na(input_file)) {
      mesg("The input file is NA for: ", filenames_with_wildcards[row], ".")
      output_entries[row] <- ''
      next
    }

    input_handle <- file(input_file, "r", blocking = FALSE)
    input_vector <- readLines(input_handle)
    last_found <- NULL
    this_found <- NULL
    all_found <- c()
    found_idx <- grepl(x = input_vector, pattern = search)
    num_hits <- 0
    if (inverse) {
      num_hits <- sum(!found_idx)
    } else {
      num_hits <- sum(found_idx)
    }
    close(input_handle)
    output_entries[row] <- num_hits
  } ## End looking at every row of the metadata
  return(output_entries)
}

#' Given an arbitrary input column, get the sum of its values.
#'
#' Note: This is not taking some important things into account yet.
#'
#' @param meta Input metadata
#' @param input_file_spec Filename specification to seek out the input data.
#' @param verbose Be verbose about what is happening?
#' @param species Choose a specific species, or all.
#' @param basedir Start searching here.
#' @param type Analysis type to search.
#' @param subtype Specific annotation type to search.
#' @param new_spec New specification resulting from this invocation.
#' @param tag GFF tag used to extract IDs.
dispatch_sum_column <- function(meta, input_file_spec, verbose = verbose,
                                species = "*", basedir = "preprocessing",
                                type = "genome", subtype = "gene",
                                new_spec = NULL, tag = "ID") {

  if (length(species) > 1) {
    output_entries <- data.frame(row.names = seq_len(nrow(meta)))
    for (s in seq_len(length(species))) {
      species_name <- species[s]
      output_entries[[species_name]] <- dispatch_sum_column(
        meta, input_file_spec, verbose = verbose,
        species = species_name, basedir = basedir, type = type)
    }
    return(output_entries)
  }

  filenames_with_wildcards <- glue(input_file_spec)
  if (isTRUE(verbose)) {
    message("Example count filename: ", filenames_with_wildcards[1], ".")
  }
  output_entries <- rep(0, length(filenames_with_wildcards))
  for (row in seq_len(nrow(meta))) {
    found <- 0
    ## Just in case there are multiple matches
    input_file <- Sys.glob(filenames_with_wildcards[row])[1]
    if (length(input_file) == 0) {
      mesg("There is no file matching: ", filenames_with_wildcards[row],
           ".")
      output_entries[row] <- ''
      next
    }
    if (is.na(input_file)) {
      mesg("The input file is NA for: ", filenames_with_wildcards[row], ".")
      output_entries[row] <- ''
      next
    }
    input_start <- try(read.table(input_file))
    if (class(input_start)[1] == "try-error") {
      warning("There was a problem reading: <", input_file, ">.")
      return(NULL)
    }

    input_df <- as.data.frame(read.table(input_file))
    if (rownames(input_df)[1] == "1") {
      rownames(input_df) <- input_df[["V1"]]
      input_df[["V1"]] <- NULL
    }
    ## FIXME: This is only useful for count tables.
    if (class(input_df) == "data.frame") {
      drop_idx <- grepl(x = rownames(input_df), pattern = "^_")
      ## I am not sure why, but this is getting coerced to a vector.
      input_df <- input_df[!drop_idx, ]
      num_hits <- sum(input_df)
    } else {
      drop_idx <- grepl(x = names(input_df), pattern = "^_")
      input_df <- input_df[!drop_idx]
      num_hits <- sum(input_df)
    }
    output_entries[row] <- num_hits
  } ## End looking at every row of the metadata
  return(output_entries)
}

#' Get the lengths of sequences from a fasta file.
#'
#' @param meta Input metadata
#' @param input_file_spec Input file specification to hunt down the
#'  file of interest.
#' @param verbose Print diagnostic information while running?
#' @param basedir Root directory containing the files/logs of metadata.
dispatch_fasta_lengths <- function(meta, input_file_spec, verbose = verbose,
                                   basedir = "preprocessing") {
  filenames_with_wildcards <- glue(input_file_spec)
  if (isTRUE(verbose)) {
    message("Example length filename: ", filenames_with_wildcards[1], ".")
  }
  output_entries <- rep(0, length(filenames_with_wildcards))
  for (row in seq_len(nrow(meta))) {
    found <- 0
    ## Just in case there are multiple matches
    input_file <- Sys.glob(filenames_with_wildcards[row])[1]
    if (length(input_file) == 0) {
      warning("There is no file matching: ", filenames_with_wildcards[row],
              ".")
      output_entries[row] <- ''
      next
    }
    if (is.na(input_file)) {
      warning("The input file is NA for: ", filenames_with_wildcards[row], ".")
      output_entries[row] <- ''
      next
    }

    dtrs <- Biostrings::readBStringSet(input_file)
    output_entries[row] <- Biostrings::width(dtrs[1])
  } ## End looking at every row of the metadata
  return(output_entries)
}

#' Pull out the filename matching an input spec
#'
#' This is useful for putting the count table name into a metadata file.
#' @param meta Input metadata
#' @param input_file_spec Input file specification to hunt down the
#'  file of interest.
#' @param verbose Print diagnostic information while running?
#' @param species Specify a species to search for, or '*' for anything.
#' @param type Some likely filename searches may be for genome
#'  vs. rRNA vs other feature types.
#' @param subtype Specific annotation type to seek.
#' @param tag GFF tag from which to seek IDs.
#' @param basedir Root directory containing the files/logs of metadata.
#' @param new_spec New specification resulting from this.
dispatch_filename_search <- function(meta, input_file_spec, verbose = verbose,
                                     species = "*", type = "genome", subtype = "gene",
                                     tag = "ID", basedir = "preprocessing", new_spec = NULL) {
  if (length(species) > 1) {
    output_entries <- data.frame(row.names = seq_len(nrow(meta)))
    for (s in seq_len(length(species))) {
      species_name <- species[s]
      output_entries[[species_name]] <- dispatch_filename_search(
        meta, input_file_spec, verbose = verbose, species = species_name,
        type = type, subtype = subtype, tag = tag, basedir = basedir, new_spec = new_spec)
    }
    return(output_entries)
  }

  filenames_with_wildcards <- glue(input_file_spec)
  new_filenames <- NULL
  if (!is.null(new_spec)) {
    new_filenames <- glue(new_spec)
    print(new_filenames)
  }
  if (isTRUE(verbose)) {
    message("Example filename: ", filenames_with_wildcards[1], ".")
  }
  output_entries <- rep(0, length(filenames_with_wildcards))
  for (row in seq_len(nrow(meta))) {
    found <- 0
    ## Just in case there are multiple matches
    input_file <- Sys.glob(filenames_with_wildcards[row])[1]
    input_new <- NULL
    if (!is.null(new_spec)) {
      input_new <- Sys.glob(new_filenames[row])[1]
    }
    if (length(input_file) == 0) {
      mesg("There is no file matching: ", filenames_with_wildcards[row], ".")
      if (!is.null(input_new)) {
        if (length(input_new) > 0) {
          input_file <- input_new
        } else {
          input_file <- ''
        }
      } else {
        input_file <- ''
      }
    }
    if (is.na(input_file)) {
      mesg("The input file is NA for: ", filenames_with_wildcards[row], ".")
      if (!is.null(input_new)) {
        if (!is.na(input_new)) {
          input_file <- input_new
        } else {
          input_file <- ''
        }
      } else {
        input_file <- ''
      }
    }
    output_entries[row] <- input_file
  }

  return(output_entries)
}

#' Pull GC content into the metadata sheet.
#'
#' As the name suggests, this only works for fasta files.
#'
#' @param meta Input metadata
#' @param input_file_spec Input file specification to hunt down the
#'  file of interest.
#' @param verbose Print diagnostic information while running?
#' @param new_spec New specification from this.
#' @param basedir Root directory containing the files/logs of metadata.
dispatch_gc <- function(meta, input_file_spec, verbose = FALSE,
                        new_spec = NULL, basedir = "preprocessing") {
  filenames_with_wildcards <- glue(input_file_spec)
  if (isTRUE(verbose)) {
    message("Example gc filename: ", filenames_with_wildcards[1], ".")
  }
  output_entries <- rep(0, length(filenames_with_wildcards))
  for (row in seq_len(nrow(meta))) {
    found <- 0
    ## Just in case there are multiple matches
    input_file <- Sys.glob(filenames_with_wildcards[row])[1]
    if (length(input_file) == 0) {
      warning("There is no file matching: ", filenames_with_wildcards[row],
              ".")
      next
    }
    if (is.na(input_file)) {
      warning("The input file is NA for: ", filenames_with_wildcards[row], ".")
      next
    }

    attribs <- sequence_attributes(input_file)
    output_entries[row] <- signif(x = attribs[["gc"]], digits = 3)
  } ## End looking at every row of the metadata
  return(output_entries)
}

#' Pull raw MD5 sums into the metadata sheet.
#'
#' @param meta Input metadata
#' @param entry_type Perhaps not used, I copy/pasted the function declaration.
#' @param input_file_spec Input file specification to hunt down the
#'  file of interest.
#' @param specification Full specification list
#' @param basedir Used to find the path to the input to md5sum
#' @param species Probably not needed
#' @param type Almost definitely not needed.
#' @param subtype Clean this too
#' @param tag And this
#' @param verbose Print diagnostic information while running?
#' @param ... Extra arguments passed through, not currently used.
#' @export
dispatch_md5 <- function(meta, entry_type, input_file_spec,
                         specification, basedir = "preprocessing",
                         species = "*", type = "genome", subtype = "gene",
                         tag = "ID", verbose = TRUE, ...) {
  input_column_name <- specification[[entry_type]][["column"]]
  input_filenames <- meta[[input_column_name]]
  output_entries <- rep("", nrow(meta))
  for (row in seq_len(nrow(meta))) {
    found <- 0
    ## Just in case there are multiple matches
    input_base_file <- input_filenames[row]
    sample_id <- rownames(meta)[row]
    input_file <- file.path(basedir, sample_id, input_base_file)
    if (isFALSE(file.exists(input_file))) {
      warning("There is no file matching: ", input_file,
              ".")
      next
    }
    if (is.na(input_file)) {
      warning("The input file is NA for: ", input_file, ".")
      next
    }

    cmd <- glue("md5sum {input_file}")
    md5 <- try(system(cmd, intern = TRUE))
    a_sum <- gsub(md5, pattern = "^([^[:space:]]+)[[:space:]]+([^[:space:]]+)$",
                  replacement = "\\1")
    mesg(cmd, " produced ", a_sum, ".")
    output_entries[row] <- a_sum
  } ## End looking at every row of the metadata
  return(output_entries)
}

#' Given two metadata columns, print a ratio.
#'
#' @param meta metadata, contains the column names!
#' @param numerator_column what it says on the tin.
#' @param denominator_column what it says on the tin.
#' @param digits Number of significant digits to keep in the output.
#' @param numerator_add Add this column to the numerator in case one
#'  needs multiple columns.
#' @param verbose unsed for the moment.
#' @param as Recast the data to this type.
#' @param species Search a specific species, or all.
dispatch_metadata_ratio <- function(meta, numerator_column = NULL,
                                    denominator_column = NULL, digits = 3,
                                    numerator_add = FALSE, verbose = FALSE,
                                    as = "numeric", species = "*") {
  if (length(species) > 1) {
    output_entries <- data.frame(row.names = seq_len(nrow(meta)))
    for (s in seq_len(length(species))) {
      species_name <- species[s]
      output_entries[[species_name]] <- dispatch_metadata_ratio(
        meta, numerator_column = numerator_column,
        denominator_column = denominator_column,
        digits = digits, numerator_add = numerator_add,
        verbose = verbose, species = species_name)
    }
    return(output_entries)
  }
  column_number <- ncol(meta)
  if (is.null(numerator_column)) {
    numerator_column <- colnames(meta)[ncol(meta)]
  }
  if (is.null(denominator_column)) {
    denominator_column <- colnames(meta)[ncol(meta) - 1]
  }
  if (!is.null(species)) {
    test_numerator <- paste0(numerator_column, "_", species)
    if (!is.null(meta[[test_numerator]])) {
      numerator_column <- test_numerator
    }
    test_denominator <- paste0(denominator_column, "_", species)
    if (!is.null(meta[[test_denominator]])) {
      denominator_column <- test_denominator
    }
    ##test_numerator_regex <- grep(x = colnames(meta), pattern = paste0("^", numerator_column))
  }
  entries <- NULL
  if (is.null(meta[[numerator_column]]) | is.null(meta[[denominator_column]])) {
    if (isTRUE(verbose)) {
      message("Missing data to calculate the ratio between: ", numerator_column,
              " and ", denominator_column, ".")
    }
  } else {
    if (isTRUE(verbose)) {
      message("The numerator column is: ", numerator_column, ".")
      message("The denominator column is: ", denominator_column, ".")
    }
    if (isTRUE(numerator_add)) {
      entries <- (as.numeric(meta[[numerator_column]]) + as.numeric(meta[[numerator_add]])) /
        as.numeric(meta[[denominator_column]])
    } else {
      entries <- as.numeric(meta[[numerator_column]]) / as.numeric(meta[[denominator_column]])
    }
    if (as == "percent") {
      entries <- entries * 100
    }
    if (!is.null(digits)) {
      entries <- signif(entries, digits)
    }
  }
  return(entries)
}

#' Generic dispatcher to hunt down useful information from logs.
#'
#' Given the metadata, a couple of regular expressions, and a filename
#' specification, this should be able to pull out the interesting
#' number(s) from one logfile per sample from the metadata.
#'
#' @param meta Input metadata.
#' @param search regex used to go hunting for the line of interest.
#' @param replace probably the same regex with parentheses in place
#'  for gsub().
#' @param input_file_spec filename extractor expression.
#' @param species Specify a species or glob it.
#' @param basedir Root directory containing the files/logs of metadata.
#' @param extraction the replacement portion of gsub(). I am thinking
#'  to make it possible to have this function return more interesting
#'  outputs if this changes, but for the moment I am sort of assuming
#'  \\1 will always suffice.
#' @param which Usually 'first', which means grab the first match and get out.
#' @param new_spec New specification resulting from this invocation.
#' @param as Coerce the output to a specific data type (numeric/character/etc).
#' @param verbose For testing regexes.
#' @param type Make explicit the type of data (genome/rRNA/Tx/etc).
#' @param backup_search A secondary search for when I change stuff.
#' @param backup_replace A secondary replace for when I change stuff.
#' @param ... Used to pass extra variables to glue for finding files.
dispatch_regex_search <- function(meta, search, replace, input_file_spec,
                                  species = "*", basedir = "preprocessing",
                                  extraction = "\\1", which = "first", new_spec = NULL,
                                  as = NULL, verbose = FALSE, type = "genome",
                                  backup_search = NULL, backup_replace = NULL,
                                  ...) {
  arglist <- list(...)
  ## if (length(arglist) > 0) {
  ##
  ## }
  if (length(species) > 1) {
    output_entries <- data.frame(row.names = seq_len(nrow(meta)))
    for (s in seq_len(length(species))) {
      species_name <- species[s]
      output_entries[[species_name]] <- dispatch_regex_search(
        meta, search, replace, input_file_spec, verbose = verbose,
        species = species_name, basedir = basedir, type = type,
        extraction = extraction, which = which, as = as, new_spec = new_spec,
        backup_search = backup_search, backup_replace = backup_replace,
        ...)
    }
    return(output_entries)
  }
  filenames_with_wildcards <- glue(input_file_spec,
                                   ...)
  new_filenames <- NULL
  if (!is.null(new_spec)) {
    new_filenames <- glue(new_spec, ...)
  }
  if (isTRUE(verbose)) {
    message("Example regex filename: ", filenames_with_wildcards[1], ".")
  }
  output_entries <- rep(0, length(filenames_with_wildcards))
  for (row in seq_len(nrow(meta))) {
    found <- 0
    ## Just in case there are multiple matches
    input_file <- Sys.glob(filenames_with_wildcards[row])[1]
    if (is.na(input_file)) {
      new_file <- NA
      if (!is.null(new_filenames[row])) {
        new_file <- Sys.glob(new_filenames[row])[1]
      }
      if (is.na(new_file)) {
        output_entries[row] <- ''
        ## The file did not exist.
        next
      } else {
        message("Setting input file to new_file: ", new_file, ".")
        input_file <- new_file
      }
    }
    if (length(input_file) == 0) {
      warning("There is no file matching: ", filenames_with_wildcards[row], ".")
      output_entries[row] <- ''
      next
    }

    input_handle <- file(input_file, "r") ## , blocking = FALSE)
    input_vector <- readLines(input_handle)
    if (length(input_vector) == 0) {
      ## Empty file, move on.
      output_entries[row] <- ''
      next
    }
    last_found <- NULL
    this_found <- NULL
    all_found <- c()
    for (i in seq_along(input_vector)) {
      if (which == "first" && found == 1) {
        output_entries[row] <- last_found
        next
      }
      input_line <- input_vector[i]
      found_boolean <- grepl(x = input_line, pattern = search)
      backup_found_boolean <- FALSE
      if (!is.null(backup_search)) {
        backup_found_boolean <- grepl(x = input_line, pattern = backup_search)
      }
      if (found_boolean || backup_found_boolean) {
        if (isTRUE(verbose)) {
          message("Found the correct line: ")
          message(input_line)
        }
        this_found <- NULL
        if (isTRUE(found_boolean)) {
          this_found <- gsub(x = input_line,
                             pattern = replace,
                             replacement = extraction)
        } else if (isTRUE(backup_found_boolean)) {
          this_found <- gsub(x = input_line,
                             pattern = backup_replace,
                             replacement = extraction)
        }
        ## Drop leading or terminal spaces.
        this_found <- gsub(x = this_found,
                           pattern = "^ +| +$", replacement = "")
        found <- found + 1
      } else {
        next
      }
      last_found <- this_found
      all_found <- c(all_found, this_found)
      output_entries[row] <- last_found
    } ## End looking at every line of the log file specified by the input file spec for this row
    close(input_handle)

    ## Handle cases where one might want to pull only the last entry in a log, or all of them.
    if (which == "last") {
      output_entries[row] <- last_found
    } else if (which == "all") {
      initial_string <- toString(all_found)
      output_entries[row] <- gsub(x = initial_string, pattern = ",", replacement = ";")
    }
  } ## End looking at every row of the metadata
  if (!is.null(as)) {
    if (as == "numeric") {
      output_entries <- as.numeric(output_entries)
    } else if (as == "character") {
      output_entries <- as.character(output_entries)
    }
  }
  return(output_entries)
}

#' Pull some information from a csv/tsv file.
#'
#' This function is a bit more generic than the others, but it grabs from a
#' column of a csv/tsv file.
#'
#' @param meta Input metadata
#' @param column Column to yank from
#' @param input_file_spec Input file specification to hunt down the
#'  file of interest.
#' @param file_type csv or tsv?
#' @param chosen_func If set, use this function to summarize the result.
#' @param species Specify a species, or glob it.
#' @param type Specify a type of search, usually genome and/or rRNA.
#' @param subtype Specific annotation type to seek.
#' @param tag GFF tag to use for IDs.
#' @param basedir Root directory containing the files/logs of metadata.
#' @param which Take the first entry, or some subset.
#' @param verbose Print diagnostic information while running?
#' @param ... Other arguments for glue.
dispatch_csv_search <- function(meta, column, input_file_spec, file_type = "csv",
                                chosen_func = NULL, species = "*", type = "genome",
                                subtype = "gene", tag = "ID",
                                basedir = "preprocessing", which = "first",
                                verbose = FALSE, ...) {
  arglist <- list(...)

  if (length(species) > 1) {
    output_entries <- data.frame(row.names = seq_len(nrow(meta)))
    for (s in seq_len(length(species))) {
      species_name <- species[s]
      output_entries[[species_name]] <- dispatch_csv_search(
        meta, column, input_file_spec, file_type = file_type,
        chosen_func = chosen_func, species = species_name,
        type = type, subtype = subtype, tag = tag,
        basedir = basedir, which = which, verbose = verbose, ...)
    }
    return(output_entries)
  }

  filenames_with_wildcards <- glue(input_file_spec, ...)
  mesg("Example csv filename: ", filenames_with_wildcards[1], ".")

  output_entries <- rep(0, length(filenames_with_wildcards))
  for (row in seq_len(nrow(meta))) {
    found <- 0
    ## Just in case there are multiple matches
    input_file <- Sys.glob(filenames_with_wildcards[row])[1]
    if (is.na(input_file)) {
      ## The file did not exist.
      output_entries[row] <- ''
      next
    }
    if (length(input_file) == 0) {
      warning("There is no file matching: ", filenames_with_wildcards[row], ".")
      output_entries[row] <- ''
      next
    }

    if (file_type == "csv") {
      input_df <- sm(readr::read_csv(input_file))
    } else if (file_type == "tsv") {
      input_df <- sm(readr::read_tsv(input_file))
    } else {
      if (isTRUE(verbose)) {
        message("Assuming csv input.")
      }
      input_df <- sm(readr::read_csv(input_file))
    }
    if (column == "last") {
      column <- ncol(input_df)
    }

    if (which == "first") {
      output_entries[row] <- input_df[1, column]
    } else if (which == "all") {
      stringified <- gsub(x = toString(input_df[[column]]), pattern = ",", replacement = ";")
      output_entries[row] <- stringified
    } else if (which == "function") {
      ## I was thinking to just use do.call() here, but I want to use na.rm=TRUE
      ## and I am not certain how to get that set up correctly with do.call.
      stringified <- 0
      if (chosen_func == "mean") {
        stringified <- mean(input_df[[column]], na.rm = TRUE)
      } else if (chosen_func == "median") {
        stringified <- median(input_df[[column]], na.rm = TRUE)
      } else if (chosen_func == "max") {
        stringified <- max(input_df[[column]], na.rm = TRUE)
      } else if (chosen_func == "min") {
        stringified <- min(input_df[[column]], na.rm = TRUE)
      }
      output_entries[row] <- stringified
    } else {
      ## Assume a number was provided for the desired row
      output_entries[row] <- input_df[which, column]
    }
  } ## End for loop
  return(output_entries)
}

#' Produce plots of metadata factor(s) of interest.
#'
#' @param exp Input dataset.
#' @param column Currently a single, but soon multiple column(s) of metadata.
#' @param second_column Or perhaps put other columns here.
#' @param norm_column Normalize the data?
#' @param type Assume a vioiln unless otherwise specified.
#' @param scale Rescale the data?
#' @return ggplot and maybe some form of useful table.
#' @export
plot_metadata_factors <- function(exp, column = "hisatsinglemapped", second_column = NULL,
                                  norm_column = NULL, type = NULL, scale = "base10") {
  condition <- plotted <- NULL ## R CMD check
  design <- as.data.frame(colData(exp))
  color_choices <- get_colors_by_condition(exp)
  meta_plot <- NULL
  if (is.null(type)) {
    meta_plot <- ggplot(design, aes(x = .data[["condition"]], y = .data[[column]])) +
      ggplot2::geom_violin(aes(fill = factor(.data[["condition"]])), scale = "width") +
      ggplot2::geom_boxplot(aes(fill = factor(.data[["condition"]])), width = 0.2, size = 0.2, outlier.size = 1.5) +
      ggplot2::geom_jitter(height = 0, width = 0.1) +
      ggplot2::scale_fill_manual(values = as.character(color_choices), guide = "none") +
      ggplot2::theme_bw(base_size = base_size) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1.0))
    if (scale == "log2") {
      meta_plot <- meta_plot +
        ggplot2::coord_trans(y = "log2")
    }
  } else if (type == "ggstats") {
    if (scale == "log2") {
      design[["plotted"]] <- log2(design[[column]])
    } else {
      design[["plotted"]] <- design[[column]]
    }
    meta_plot <- ggstatsplot::ggbetweenstats(data = design, x = condition, y = plotted)
  }
  return(meta_plot)
}

#' Plot metadata factors as a sankey diagram.
#'
#' This provides two implementations of a sankey plot, one interactive
#' and one using ggplot2.
#'
#' @param design Metadata from which to extract the categories/numbers.
#' @param factors Factors/columns in the metadata to count and plot.
#' @param fill Use either the current or next node for coloring the transitions.
#' @param font_size Chosen font size, perhaps no longer needed?
#' @param node_width Make nodes more or less rectangular with this.
#' @param color_choices Either a named vector of states and colors, or NULL
#'  (in which case it will use viridis.)
#' @param drill_down When true, this will end in the product of the
#'  factor levels number of final states. (e.g. if there are 2 sexes,
#'  3 visits, and 4 genotypes, there will be 2, 6, 24 states going
#' from right to left).  If FALSE, there will be 2,3,4 states going
#' from right to left.
#' @return List containing a couple of plots, one interactive, one gg.
#' @export
plot_meta_sankey <- function(design, factors = c("condition", "batch"), fill = "node",
                             font_size = 18, node_width = 30,
                             color_choices = NULL,
                             drill_down = TRUE) {
  ## Recasting to a dataframe in case I pass a SE with a datatable/tibble/etc
  design <- as.data.frame(design)
  found <- factors %in% colnames(design)
  if (sum(found) < length(factors)) {
    missing <- factors[!found]
    message("These columns are not in the metadata: ", toString(missing))
    factors <- factors[found]
  }

  permutations <- c()
  states <- list()
  for (f in factors) {
    design[[f]] <- as.character(design[[f]])
    design[[f]] <- gsub(x = design[[f]], pattern = " ", replacement = "_")
    design[[f]] <- as.factor(design[[f]])
    state_levels <- levels(as.factor(design[[f]]))
    states[[f]] <- state_levels
    new_permutations <- tidyr::expand_grid(!!!states) |> purrr::pmap_chr(paste)
    permutations <- c(permutations, new_permutations)
  }

  start_node <- c("all", "0")
  names(start_node) <- c("name", "node")
  my_nodes <- data.frame(name = permutations)
  my_nodes[["node"]] <- rownames(my_nodes)
  my_nodes <- rbind(start_node, my_nodes)

  my_links <- data.frame(row.names = 0)
  my_links[["source"]] <- 0
  my_links[["target"]] <- 0
  my_links[["value"]] <- nrow(design)

  result <- list("all" = nrow(design))
  for (p in seq_len(length(permutations))) {
    element <- permutations[p]
    pieces <- strsplit(x = element, split = " ")[[1]]

    working_meta <- design
    for (cat_num in seq_len(length(pieces))) {
      category <- pieces[cat_num]
      factor <- factors[cat_num]
      idx <- working_meta[[factor]] == category
      working_meta <- working_meta[idx, ]
    }
    if (nrow(working_meta) > 0) {
      result[[element]] <- nrow(working_meta)

      target_node_idx <- my_nodes[["name"]] == element
      target_node <- my_nodes[target_node_idx, "node"]
      if (length(pieces) > 1) {
        source_node_pieces <- pieces[1:length(pieces) - 1]
        source_node_name <- stringi::stri_paste(source_node_pieces, collapse = " ")
        source_node_idx <- my_nodes[["name"]] == source_node_name
        source_node <- my_nodes[source_node_idx, "node"]
      } else {
        source_node <- "0"
      }
      link <- c(source_node, target_node, nrow(working_meta))
      my_links <- rbind(my_links, link)
    }
  }

  my_links <- my_links[-1, ]
  my_links[["value"]] <- as.numeric(my_links[["value"]])
  my_links[["source"]] <- as.numeric(my_links[["source"]])
  my_links[["target"]] <- as.numeric(my_links[["target"]])

  links_to_nodes <- merge(my_links, my_nodes, by.x = "target",
                          by.y = "node")
  links_to_nodes <- links_to_nodes[, c("name", "value")]

  sub_design <- design[, c(factors)]
  if (isTRUE(drill_down)) {
    for (col in seq_len(ncol(sub_design))) {
      if (col == 1) {
        next
      }
      col_name <- colnames(sub_design)[col]
      previous_col <- colnames(sub_design)[col - 1]
      new_values <- paste0(sub_design[[previous_col]], " ",
                           sub_design[[col_name]])
      sub_design[[col_name]] <- new_values
    }
  }

  plot_df <- sub_design %>%
    ggsankey::make_long(factors)
  plot_df[["name"]] <- plot_df[["node"]]
  plot_df <- merge(plot_df, links_to_nodes, by.x = "node", by.y = "name")
  plot_df[["name"]] <- paste0(plot_df[["name"]], ":",
                              plot_df[["value"]])

  color_fact <- NULL
  if (!is.null(color_choices)) {
    color_levels <- levels(as.factor(plot_df[["node"]]))
    all_colors <- unlist(color_choices)
    names(all_colors) <- gsub(x = names(all_colors),
                              pattern = ".*\\.", replacement = "")
    color_fact_idx <- names(all_colors) %in% color_levels
    color_fact <- all_colors[color_fact_idx]
    if (isTRUE(drill_down)) {
      ## Set the plot name and color names:
      plot_df[["name"]] <- gsub(x = plot_df[["node"]],
        pattern = "^.* (\\w+$)",
        replacement = "\\1")
      plot_df[["name"]] <- paste0(plot_df[["name"]], ":",
                                  plot_df[["value"]])
      color_level_suffixes <- gsub(x = color_levels,
        pattern = "^.* (\\w+$)",
        replacement = "\\1")
      color_fact_idx <- names(all_colors) %in% color_level_suffixes
      color_suffix_fact <- all_colors[color_fact_idx]
      new_color_fact <- rep("#000000", length(color_level_suffixes))
      names(new_color_fact) <- color_levels
      for (col in seq_len(length(new_color_fact))) {
        color_name <- color_level_suffixes[col]
        color_suffix <- as.character(color_suffix_fact[color_name])
        new_color_fact[col] <- color_suffix
      }
      color_fact <- new_color_fact
    }
  }
  retlist <- list(
    "design" = design,
    "factors" = factors,
    "observed_nodes" = unique(plot_df[["node"]]))

  if (fill == "node") {
    ggplt <- ggplot(plot_df, aes(x = x, next_x = next_x, node = node,
                                 next_node = next_node, fill = factor(node), label = name))
  } else if (fill == "next") {
    message("Filling to next node?")
    ggplt <- ggplot(plot_df, aes(x = x, next_x = next_x, node = node,
                                 next_node = next_node, fill = factor(next_node), label = name))
  }
  ggplt <- ggplt +
    ggsankey::geom_sankey(flow.alpha = 0.6,
                          node.color = "gray30") +
    ggsankey::geom_sankey_label(color = 1, fill = "white", size = 3.5)

  ## I want to figure out how to set up my own colors...
  if (is.null(color_choices)) {
    ggplt <- ggplt +
      ggplot2::scale_fill_viridis_d() +
      ggsankey::theme_sankey(base_size = 18) +
      ggplot2::labs(x = NULL) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5))
  } else {
    ggplt <- ggplt +
      ggplot2::scale_fill_manual(values = color_fact) +
      ggsankey::theme_sankey(base_size = font_size) +
      ggplot2::labs(x = NULL) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5))
  }
  retlist[["ggplot"]] <- ggplt
  class(retlist) <- "meta_sankey"
  return(retlist)
}
setGeneric("plot_meta_sankey")

#' Print a metadata sankey plot.
#'
#' @param x List containing the table of connected nodes and a ggplot2
#'  sankey.
#' @param ... Other args to match the generic.
#' @export
print.meta_sankey <- function(x, ...) {
  summary_string <- glue("A sankey plot describing the metadata of {nrow(x[['design']])} samples,
including {length(x[['observed_nodes']])} out of {length(x[['permutations']])} nodes \\
and traversing metadata factors:
{toString(x[['factors']])}.")
  message(summary_string)
  plot(x[["ggplot"]])
  return(invisible(x))
}

#' Feed an expt to a sankey plot.
#' @param design Metadata from which to extract the categories/numbers.
#' @param factors Factors/columns in the metadata to count and plot.
#' @param fill Use either the current or next node for coloring the transitions.
#' @param font_size Chosen font size, perhaps no longer needed?
#' @param node_width Make nodes more or less rectangular with this.
#' @param color_choices Either a named vector of states and colors, or NULL
#'  (in which case it will use viridis.)
#' @param drill_down When true, this will end in the product of the
#'  factor levels number of final states. (e.g. if there are 2 sexes,
#'  3 visits, and 4 genotypes, there will be 2, 6, 24 states going
#' from right to left).  If FALSE, there will be 2,3,4 states going
#' from right to left.
#' @export
setMethod(
  "plot_meta_sankey", signature = signature(design = "expt"),
  definition = function(design, factors = c("condition", "batch"),
                        fill = "node", font_size = 18, node_width = 30,
                        color_choices = NULL, drill_down = TRUE) {
    design <- pData(design)
    plot_meta_sankey(design, factors = factors, fill = fill,
                     font_size = font_size, node_width = node_width,
                     color_choices = color_choices, drill_down = drill_down)
  })

#' Using plot_meta_sankey on a SE is basically the same as an expt.
#'
#' @inheritParams plot_meta_sankey
setMethod(
  "plot_meta_sankey", signature = signature(design = "SummarizedExperiment"),
  definition = function(design, factors = c("condition", "batch"),
                        fill = "node", font_size = 18, node_width = 30,
                        color_choices = NULL, drill_down = TRUE) {
    design <- colData(design)
    plot_meta_sankey(design, factors = factors, fill = fill, font_size = font_size,
                     node_width = node_width, color_choices = color_choices,
                     drill_down = drill_down)
  })

## I might want to just delete this, I think the ggplot version is better.
plot_meta_interactive_sankey <- function(design, factors = c("condition", "batch"), fill = "node",
                                         font_size = 18, node_width = 30,
                                         color_choices = NULL, html = NULL,
                                         drill_down = TRUE) {

  found <- factors %in% colnames(design)
  if (sum(found) < length(factors)) {
    missing <- factors[!found]
    message("These columns are not in the metadata: ", toString(missing))
    factors <- factors[found]
  }

  permutations <- c()
  states <- list()
  for (f in factors) {
    state_levels <- levels(as.factor(design[[f]]))
    states[[f]] <- state_levels
    new_permutations <- tidyr::expand_grid(!!!states) |> purrr::pmap_chr(paste)
    permutations <- c(permutations, new_permutations)
  }

  start_node <- c("all", "0")
  names(start_node) <- c("name", "node")
  my_nodes <- data.frame(name = permutations)
  my_nodes[["node"]] <- rownames(my_nodes)
  my_nodes <- rbind(start_node, my_nodes)

  my_links <- data.frame(row.names = 0)
  my_links[["source"]] <- 0
  my_links[["target"]] <- 0
  my_links[["value"]] <- nrow(design)

  result <- list("all" = nrow(design))
  for (p in seq_len(length(permutations))) {
    element <- permutations[p]
    pieces <- strsplit(x = element, split = " ")[[1]]

    working_meta <- design
    for (cat_num in seq_len(length(pieces))) {
      category <- pieces[cat_num]
      factor <- factors[cat_num]
      idx <- working_meta[[factor]] == category
      working_meta <- working_meta[idx, ]
    }
    if (nrow(working_meta) > 0) {
      result[[element]] <- nrow(working_meta)

      target_node_idx <- my_nodes[["name"]] == element
      target_node <- my_nodes[target_node_idx, "node"]
      if (length(pieces) > 1) {
        source_node_pieces <- pieces[1:length(pieces) - 1]
        source_node_name <- stringi::stri_paste(source_node_pieces, collapse = " ")
        source_node_idx <- my_nodes[["name"]] == source_node_name
        source_node <- my_nodes[source_node_idx, "node"]
      } else {
        source_node <- "0"
      }
      link <- c(source_node, target_node, nrow(working_meta))
      my_links <- rbind(my_links, link)
    }
  }

  my_links <- my_links[-1, ]
  my_links[["value"]] <- as.numeric(my_links[["value"]])
  my_links[["source"]] <- as.numeric(my_links[["source"]])
  my_links[["target"]] <- as.numeric(my_links[["target"]])

  links_to_nodes <- merge(my_links, my_nodes, by.x = "target",
                          by.y = "node")
  links_to_nodes <- links_to_nodes[, c("name", "value")]


  plt <- NULL
  if (!is.null(html)) {
    plt = networkD3::sankeyNetwork(Links = my_links, Nodes = my_nodes,
                                   Source = "source", Target = "target",
                                   Value = "value", fontSize = font_size, nodeWidth = node_width)
  }
  retlist <- list(
    "permutations" = permutations,
    "plot" = plt)
  if (!is.null(html)) {
    retlist[["html"]] <- htmlwidgets::saveWidget(plt, file = html, selfcontained = TRUE)
  }
  return(retlist)
}

#' Given a table of meta data, read it in for use by create_expt().
#'
#' Reads an experimental design in a few different formats in preparation for
#' creating an expt.
#'
#' @param file Csv/xls file to read.
#' @param sep Used by read.csv, the separator
#' @param header Used by read.csv, is there a header?
#' @param sheet Used for excel/etc, which sheet to read?
#' @param comment Skip rows starting with this (in the first cell of the row if not a text file).
#' @param sanitize Sanitize the results?
#' @param sanitize_underscore including underscores!?
#' @param ... Arguments for arglist, used by sep, header and similar
#'  read_csv/read.table parameters.
#' @return Df of metadata.
#' @seealso [openxlsx] [readODS]
#' @export
read_metadata <- function(file, sep = ",", header = TRUE, sheet = 1, comment = "#", sanitize = TRUE, sanitize_underscore = FALSE,
                          ...) {
  arglist <- list(...)

  extension <- tools::file_ext(file)
  if (extension == "csv") {
    definitions <- read.csv(file = file, comment.char = comment,
                            sep = sep, header = header)
  } else if (extension == "tsv") {
    definitions <- try(readr::read_tsv(file, ...))
  } else if (extension == "xlsx") {
    ## xls = loadWorkbook(file, create = FALSE)
    ## tmp_definitions = readWorksheet(xls, 1)
    definitions <- try(openxlsx::read.xlsx(xlsxFile = file, sheet = sheet,
                                           detectDates = TRUE))

    if (class(definitions)[1] == "try-error") {
      definitions <- try(openxlsx::read.xlsx(xlsxFile = file, sheet = sheet,
                                             detectDates = FALSE))
      if (class(definitions)[1] == "try-error") {
        stop("Unable to read the metadata file: ", file)
      }
    }
  } else if (extension == "xls") {
    ## This is not correct, but it is a start
    definitions <- readxl::read_excel(path = file, sheet = sheet)
  } else if (extension == "ods") {
    definitions <- readODS::read_ods(path = file, sheet = sheet)
  } else {
    definitions <- read.table(file = file, sep = sep,
                              header = header)
  }

  if (!is.null(comment)) {
    first_column <- as.data.frame(definitions)[[1]]
    commented <- grepl(x = definitions[[1]], pattern = "^#")
    ## keep the un-commented lines.
    definitions <- definitions[! commented, ]
  }

  if (isTRUE(sanitize)) {
    colnames(definitions) <- gsub(x = colnames(definitions), pattern = "\\s+",
                                  replacement = "_")
    if (isTRUE(sanitize_underscore)) {
      colnames(definitions) <- tolower(gsub(pattern = "[[:punct:]]",
                                            replacement = "",
                                            x = colnames(definitions)))
    } else {
      colnames(definitions) <- tolower(gsub(pattern = "[^_[:^punct:]]",
                                            replacement = "",
                                            x = colnames(definitions), perl = TRUE))
    }
    ## I recently received a sample sheet with a blank sample ID column name...
    empty_idx <- colnames(definitions) == ""
    colnames(definitions)[empty_idx] <- "empty"
    colnames(definitions) <- make.names(colnames(definitions), unique = TRUE)
  }
  return(definitions)
}

#' Make an archive using a column from the metadata.
#'
#' I am hoping this will be useful for either backing up count tables or making
#' containerized versions of analyses.
#'
#' @param meta dataframe of the good stuff.
#' @param column Column containing filenames to archive.
#' @param output Output prefix for the tarball's name.
#' @param compression Actually, this might be a mistake, I think utils::tar takes 'gzip', not 'gz'?
#' @export
tar_meta_column <- function(meta, column = "hisatcounttable", output = NULL, compression = "xz") {
  start_list <- meta[[column]]
  file_list <- Sys.glob(start_list)
  include_list <- c()
  for (f in seq_along(file_list)) {
    file <- file_list[f]
    if (nchar(file) == 0) {
      mesg("There is no file matching: ", start_list[f], ".")
      next
    }
    if (is.na(file)) {
      mesg("The input file is NA for: ", start_list[f], ".")
      next
    }
    include_list <- c(include_list, file_list[f])
  }

  if (is.null(output)) {
    if (is.null(compression)) {
      output <- glue("{column}.tar")
    } else {
      output <- glue("{column}.tar.{compression}")
    }
  }
  mesg("Creating a tar archive of ", toString(include_list))
  if (is.null(compression)) {
    tarred <- utils::tar(output, files = include_list)
  } else {
    tarred <- utils::tar(output, files = include_list, compression = compression)
  }
  retlist <- list(
    "output" = output,
    "files_included" = include_list,
    "files_requested" = start_list)
  class(retlist) <- "meta_tarball"
  return(retlist)
}
setGeneric("tar_meta_column")

#' Make a tarball using a metadata column given an xlsx file.
#'
#' @param meta dataframe of the good stuff.
#' @param column Column containing filenames to archive.
#' @param output Output prefix for the tarball's name.
#' @param compression Actually, this might be a mistake, I think utils::tar takes 'gzip', not 'gz'?
#' @export
setMethod(
  "tar_meta_column", signature = signature(meta = "character"),
  definition = function(meta, column = "hisatcounttable",
                        output = NULL, compression = "xz") {
    meta <- read_metadata(meta)
    ## Someone snuck some spaces into my most recent sample sheet file column!
    meta <- sanitize_metadata(meta, columns = column, punct = FALSE, lower = FALSE, spaces = TRUE)
    result <- tar_meta_column(meta, column = column,
                              output = output, compression = compression)
    return(result)
  })

#' Generate an assembly annotation specification for use by gather_preprocessing_metadata()
#'
#' This is the default set of files/information that will be sought.  It is a bit much.
#' Each name of the returned list is one column in the final metadata.  The values within that
#' name are the relevant parameters for the associated dispatcher.
#'
#' The assembly pipeline I wrote for which this was written does the following:
#' 1.  Trimomatic (the assemblies I was doing were miseq phage).
#' 2.  Fastqc the trimmed reads.
#' 3.  Racer to correct sequencer-based errors.
#' 4.  Perform an initial classification with kraken vs. the standard database. (thus if there is contamination we can pick it up)
#' 5.  Use kraken to make a hypotehtical host for the phage and filter it.
#' 6.  Classify the remaining sequence with kraken vs a viral database.
#' 7.  Generate an initial assembly via unicycler.
#' 8.  Depth-filter said assembly.
#' 9.  Use Blast to search the ICTV for likely taxonomy.
#' 10. Count ORFs to define the +/- strands.
#' 11. Use Phageterm to define the DTRs and/or reorient the genome.
#' 12. Perform a taxonomy search on the assembled genome via phastaf (thus we can see if it is segmented or multiple genomes).
#' 13. Calculate coverage on a per-nucleotide basis.
#' 14. Search for likely terminases, and reorient the genome if phageterm (#11) failed.
#' 15. Create an initial annotation genbank file via prokka.
#' 16. Supplement the prokka ORFs via a trained prodigal run.
#' 17. Supplement them again via a promiscuous run of glimmer.
#' 18. Use phanotate as the arbiter of 'correct' phage ORFs. (e.g. the ORFs from #15-17 will only be used if they agree with and/or do not interfere with these).
#' 19. Merge the results from #15-18 into a single set of ORFs/genbank.
#' 20. Calculate the assembly kmer content via jellyfish.
#' 21. Look for t(m)RNAs via aragorn.
#' 22. Look for tRNAs via tRNAscan.
#' 23. Perform the set of blast/etc searches defined by trinotate.
#' 24. Look for MDR genes via abricate.
#' 25. Perform the set of blast/etc searches defined by interproscan.
#' 26. Cross reference the genome against the extant restriction enzyme catalog.
#' 27. Calculate the codon adaptation index of each ORF against the putative host from #5.
#' 28. Search for phage promoters.
#' 29. Search for Rho termination signals.
#' 30. Attempt to classify the phage's likelihood to be lysogenic/lytic via bacphlip.
#' 31. Search for strong RNA secondary structures via RNAfold.
#' 32. Merge the annotations collected from #21-29 into a larger genbank file.
#' 33. Repeat #32, but this time with feeling. (#32 adds comments with confidence intervals, this strips those out).
#' 34. Make an initial visualization of the assembly via cgview.
#' 35. Collect all the most likely useful stuff from above into a single archive.
#' 36. Clean up the mess.
#' @export
make_assembly_spec <- function() {
  specification <- list(
    ## First task performed is pretty much always trimming
    "assembly_fasta_nt" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/unicycler_assembly.fasta"),
    "assembly_genbank_annotated" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*mergeannot/{meta[['sampleid']]}.gbk"),
    "assembly_genbank_stripped" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*mergeannot/{meta[['sampleid']]}_stripped.gbk"),
    "assembly_cds_amino_acids" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*merge_cds_predictions/{meta[['sampleid']]}.faa"),
    "assembly_cds_nucleotides" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*merge_cds_predictions/{meta[['sampleid']]}.ffn"),
    "assembly_gff" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*merge_cds_predictions/{meta[['sampleid']]}.gff"),
    "assembly_tsv" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*merge_cds_predictions/{meta[['sampleid']]}.tsv"),
    "assembly_xls" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*mergeannot/{meta[['sampleid']]}.xlsx"),
    "input_r1" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/scripts/*trim_*.sh"),
    "input_r2" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/scripts/*trim_*.sh"),
    "trimomatic_input" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*trimomatic/*-trimomatic.stderr"),
    "trimomatic_output" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*trimomatic/*-trimomatic.stderr"),
    "trimomatic_ratio" = list(
      "column" = "trimomatic_percent"),
    ## Second task is likely error correction
    ##"racer_changed" = list(
    ##    "file" = "preprocessing/{meta[['sampleid']]}/outputs/*racer/racer.out"),
    "host_filter_species" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/host_species.txt"),
    ## After those, things can get pretty arbitrary...
    "hisat_genome_single_concordant" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*/hisat*_*.stderr"),
    "hisat_genome_multi_concordant" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*/hisat*_*.stderr"),
    "hisat_genome_single_all" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*/hisat*_*.stderr"),
    "hisat_genome_multi_all" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*/hisat*_*.stderr"),
    "hisat_count_table" = list(
      ## Note, I changed the output slightly for hisat in recent (202404) cyoa invocations
      ## So, for the moment I will make this less stringent to pick up both file names.
      ## Ergo the extra glob.
      ##"file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_{type}*.count.xz"),
      "file_new" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_*{type}*_{subtype}_{tag}_fcounts.csv.xz",
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_*{type}*_{subtype}_{tag}.count.xz"),
    "hisat_alignment" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_*{type}*.bam"),
    "deeptools_coverage" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_*{type}*.bw"),
    "jellyfish_count_table" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*jellyfish_*/*_matrix.csv.xz"),
    "jellyfish_observed" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*jellyfish_*/*_matrix.csv.xz"),
    "kraken_viral_classified" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*kraken_viral*/kraken.stderr"),
    "kraken_viral_unclassified" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*kraken_viral*/kraken.stderr"),
    "kraken_first_viral_species" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*kraken_viral*/kraken_report.txt"),
    "kraken_first_viral_species_reads" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*kraken_viral*/kraken_report.txt"),
    "kraken_standard_classified" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*kraken_standard*/kraken.stderr"),
    "kraken_standard_unclassified" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*kraken_standard*/kraken.stderr"),
    "kraken_first_standard_species" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*kraken_standard*/kraken_report.txt"),
    "kraken_first_standard_species_reads" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*kraken_standard*/kraken_report.txt"),
    "possible_host_species" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*filter_kraken_host/*.log"),
    "pernt_mean_coverage" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/??assembly_coverage_*/base_coverage.tsv"),
    "pernt_median_coverage" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/??assembly_coverage_*/base_coverage.tsv"),
    "pernt_min_coverage" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/??assembly_coverage_*/base_coverage.tsv"),
    "pernt_max_coverage" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/??assembly_coverage_*/base_coverage.tsv"),
    "salmon_stranded" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*salmon_*/salmon_*.stderr"),
    "salmon_mapped" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*salmon_*/salmon_*.stderr"),
    "shovill_contigs" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*shovill_*/shovill.log"),
    "shovill_length" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*shovill_*/shovill.log"),
    "shovill_estlength" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*shovill_*/shovill.log"),
    "shovill_minlength" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*shovill_*/shovill.log"),
    "unicycler_lengths" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*unicycler/*final_assembly.fasta"),
    "unicycler_relative_coverage" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*unicycler/*final_assembly.fasta"),
    "filtered_relative_coverage" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/??filter_depth/final_assembly.fasta"),
    "phastaf_num_hits" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*phastaf_*/phage.bed"),
    "phageterm_dtr_length" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*phageterm_*/phageterm_final_dtr.fasta"),
    "prodigal_positive_strand" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*prodigal_*/predicted_cds.gff"),
    "prodigal_negative_strand" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*prodigal_*/predicted_cds.gff"),
    "glimmer_positive_strand" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*glimmer/glimmer3.predict"),
    "glimmer_negative_strand" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*glimmer/glimmer3.predict"),
    "phanotate_positive_strand" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*phanotate/*_phanotate.tsv.xz"),
    "phanotate_negative_strand" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*phanotate/*_phanotate.tsv.xz"),
    "final_gc_content" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*prokka/{meta[['sampleid']]}.fna"),
    "interpro_signalp_hits" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
    "interpro_phobius_hits" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
    "interpro_pfam_hits" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
    "interpro_tmhmm_hits" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
    "interpro_cdd_hits" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
    "interpro_smart_hits" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
    "interpro_gene3d_hits" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
    "interpro_superfamily_hits" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*interproscan_*/{meta[['sampleid']]}.faa.tsv"),
    "tRNA_hits" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*prokka_*/{meta[['sampleid']]}.log"),
    "aragorn_tRNAs" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*aragorn/aragorn.txt"),
    "ictv_taxonomy" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*classify_*/*_filtered.tsv"),
    "ictv_accession" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*classify_*/*_filtered.tsv"),
    ##"ictv_family" = list(
    ##    "file" = "{basedir}/{meta[['sampleid']]}/outputs/*classify_*/*_filtered.tsv"),
    "ictv_genus" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*classify_*/*_filtered.tsv"),
    "notes" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/notes.txt")
  ) ## Finish the list on a separate line because this is already absurd enough.
  return(specification)
}

#' Generate a RNASeq specification for use by gather_preprocessing_metadata()
#'
#' This currently assumes the set of tools used by one doing RNASeq to be trimomatic,
#' fastqc, hisat2, and htseq.
#'
#' @param umi Include entries for umi-barcoded samples?
#' @export
make_rnaseq_spec <- function(umi = FALSE) {
  specification <- list(
    ## First task performed is pretty much always trimming
    "trimomatic_input" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*trimomatic/*-trimomatic.stderr"),
    "trimomatic_output" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*trimomatic/*-trimomatic.stderr"),
    "trimomatic_ratio" = list(
      "column" = "trimomatic_percent"),
    "fastqc_pct_gc" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*fastqc/*_fastqc/fastqc_data.txt"),
    "fastqc_most_overrepresented" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*fastqc/*_fastqc/fastqc_data.txt"),
    "kraken_matrix" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*kraken_*/kraken_report_matrix.tsv"),
    "hisat_rrna_input_reads" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*rRNA*.stderr"),
    "hisat_rrna_single_concordant" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*rRNA*.stderr"),
    "hisat_rrna_multi_concordant" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*rRNA*.stderr"),
    "hisat_rrna_pct_vs_trimmed" = list(
      "column" = "hisat_rrna_pct_vs_trimmed"),
    "hisat_rrna_percent_log" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*{type}*.stderr"),
    "hisat_genome_input_reads" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*{type}*.stderr"),
    "hisat_genome_single_concordant" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*{type}*.stderr"),
    "hisat_genome_multi_concordant" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*{type}*.stderr"),
    "hisat_genome_single_all" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*{type}*.stderr"),
    "hisat_genome_multi_all" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*{type}*.stderr"),
    "hisat_unmapped" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*{type}*.stderr"),
    "hisat_genome_percent_log" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*{type}*.stderr"),
    "hisat_genome_pct_vs_trimmed" = list(
      "column" = "hisat_genome_pct_vs_trimmed"),
    "hisat_observed_genes" = list(
      "file_new" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_{type}*_{subtype}_{tag}_fcounts.csv.xz",
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_{type}*_{subtype}_{tag}.count.xz"),
    "hisat_sum_genes" = list(
      "file_new" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_{type}*_{subtype}_{tag}_fcounts.csv.xz",
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_{type}*_{subtype}_{tag}.count.xz"),
    "hisat_sum_genes_pct" = list(
      "column" = "hisat_sum_genes_pct"),
    "hisat_observed_mean_exprs" = list(
      "file_new" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_{type}*_{subtype}_{tag}_fcounts.csv.xz",
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_{type}*_{subtype}_{tag}.count.xz"),
    "hisat_observed_median_exprs" = list(
      "file_new" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_{type}*_{subtype}_{tag}_fcounts.csv.xz",
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_{type}*_{subtype}_{tag}.count.xz"),
    "hisat_alignment" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_*{type}*.bam"),
    "deeptools_coverage" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_*{type}*.bw"),
    "salmon_stranded" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*salmon_{species}*/salmon_*.stderr"),
    "salmon_mapped" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*salmon_{species}*/salmon_*.stderr"),
    "salmon_percent" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*salmon_{species}*/salmon_*.stderr"),
    "salmon_observed_genes" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*salmon_{species}*/quant.sf"),
    ## The end should be the various output filenames.
    "input_r1" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/scripts/*trim_*.sh"),
    "input_r2" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/scripts/*trim_*.sh"),
    ## 202410: Given recent changes to FC/FH/FD handling, this will require some changes...
    ## Make cyoa create input_r1= output_r1= lines to make this easier.
    ##    "trimmed_r1" = list(
    ##      "file" = "{basedir}/{meta[['sampleid']]}/scripts/*trim_*.sh"),
    ##    "trimmed_r2" = list(
    ##      "file" = "{basedir}/{meta[['sampleid']]}/scripts/*trim_*.sh"),
    "hisat_count_table" = list(
      ## Note, I changed the output slightly for hisat in recent (202404) cyoa invocations
      ## So, for the moment I will make this less stringent to pick up both file names.
      ## Ergo the extra glob.
      ##"file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_{type}*.count.xz"),
      "file_new" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_*{type}*_{subtype}_{tag}_fcounts.csv.xz",
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_*{type}*_{subtype}_{tag}.count.xz"),
    "salmon_count_table" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*salmon_{species}*/quant.sf"),
    "bbmap_coverage_stats" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*bam2coverage*/coverage.tsv.xz"),
    "bbmap_coverage_per_nt" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*bam2coverage*/base_coverage.tsv.xz")
  )

  if (isTRUE(umi)) {
    specification[["umi_extract_r1_output"]] = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*umi_tools/r1_extracted.fastq.gz")
    specification[["umi_extract_r2_output"]] = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*umi_tools/r2_extracted.fastq.gz")
    specification[["umi_extract_matches"]] = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*umi_tools/extract.log")
    specification[["umi_dedup_output_bam"]] = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*umi_dedup/umi_tools_deduplicated.bam")
    specification[["umi_dedup_output_count"]] = list(
      "file_new" = "{basedir}/{meta[['sampleid']]}/outputs/*umi_dedup/umi_tools_deduplicated*.csv.xz",
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*umi_dedup/umi_tools_deduplicated*.count.xz")
    specification[["umi_dedup_chimeric"]] = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*umi_dedup/umi_dedup.stderr")
    specification[["umi_dedup_num_reads_in"]] = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*umi_dedup/umi_dedup.stderr")
    specification[["umi_dedup_num_reads_out"]] = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*umi_dedup/umi_dedup.stderr")
    specification[["umi_dedup_pct_reads"]] = list(
      "column" = "umi_dedup_pct_reads")
    specification[["umi_dedup_deduplicated_positions"]] = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*umi_dedup/umi_dedup.stderr")
    specification[["umi_dedup_mean_umi_per_pos"]] = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*umi_dedup/umi_dedup.stderr")
    specification[["umi_dedup_max_umi_per_pos"]] = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*umi_dedup/umi_dedup.stderr")
  }

  return(specification)
}

make_rnaseq_multibioproject <- function() {
  spec <- list(
    ## First task performed is pretty much always trimming
    ## "input_r1" = list(
    ##    "file" = "{basedir}/{meta[['sampleid']]}/scripts/*trim_*.sh"),
    ## "input_r2" = list(
    ##    "file" = "{basedir}/{meta[['sampleid']]}/scripts/*trim_*.sh"),
    "trimomatic_input" = list(
      "file" = "{basedir}/*/{meta[['sampleid']]}/outputs/*trimomatic/*-trimomatic.stderr"),
    "trimomatic_output" = list(
      "file" = "{basedir}/*/{meta[['sampleid']]}/outputs/*trimomatic/*-trimomatic.stderr"),
    "trimomatic_ratio" = list(
      "column" = "trimomatic_percent"),
    "salmon_stranded" = list(
      "file" = "{basedir}/*/{meta[['sampleid']]}/outputs/*salmon_{species}*/salmon_*.stderr"),
    "salmon_count_table" = list(
      "file" = "{basedir}/*/{meta[['sampleid']]}/outputs/*salmon_{species}*/quant.sf"),
    "salmon_mapped" = list(
      "file" = "{basedir}/*/{meta[['sampleid']]}/outputs/*salmon_{species}*/salmon_*.stderr"),
    "salmon_percent" = list(
      "file" = "{basedir}/*/{meta[['sampleid']]}/outputs/*salmon_{species}*/salmon_*.stderr")
  )
  return(spec)
}

#' Generate a DNASeq specification for use by gather_preprocessing_metadata()
#'
#' This currently assumes the set of tools used by one doing RNASeq to be trimomatic,
#' fastqc, hisat2, htseq, freebayes, and my variant post-processor.
#' @export
make_dnaseq_spec <- function() {
  specification <- list(
    "trimomatic_input" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*trimomatic/*-trimomatic.stderr"),
    "trimomatic_output" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*trimomatic/*-trimomatic.stderr"),
    "trimomatic_ratio" = list(
      "column" = "trimomatic_percent"),
    "fastqc_pct_gc" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*fastqc/*_fastqc/fastqc_data.txt"),
    "fastqc_most_overrepresented" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*fastqc/*_fastqc/fastqc_data.txt"),
    "hisat_genome_single_concordant" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*genome*.stderr"),
    "hisat_genome_multi_concordant" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*genome*.stderr"),
    "hisat_genome_single_all" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*genome*.stderr"),
    "hisat_genome_multi_all" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*genome*.stderr"),
    "hisat_genome_percent_log" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/hisat*_*{type}*.stderr"),
    "hisat_genome_pct_vs_trimmed" = list(
      "column" = "hisat_genome_percent_log"),
    "hisat_alignment" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_*{type}*.bam"),
    "deeptools_coverage" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_*{type}*.bw"),
    "gatk_unpaired" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/deduplication_stats.txt"),
    "gatk_paired" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/deduplication_stats.txt"),
    "gatk_supplementary" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/deduplication_stats.txt"),
    "gatk_unmapped" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/deduplication_stats.txt"),
    "gatk_unpaired_duplicates" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/deduplication_stats.txt"),
    "gatk_paired_duplicates" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/deduplication_stats.txt"),
    "gatk_paired_opt_duplicates" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/deduplication_stats.txt"),
    "gatk_duplicate_pct" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/deduplication_stats.txt"),
    "gatk_libsize" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/deduplication_stats.txt"),
    "freebayes_observed" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/all_tags*"),
    ## Put the various output files at the end.
    "input_r1" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/scripts/*trim_*.sh"),
    "input_r2" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/scripts/*trim_*.sh"),
    "freebayes_observed_file" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/all_tags*"),
    "hisat_count_table" = list(
      ## Note, I changed the output slightly for hisat in recent (202404) cyoa invocations
      ## So, for the moment I will make this less stringent to pick up both file names.
      ## Ergo the extra glob.
      ##"file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_{type}*.count.xz"),
      "file_new" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_*{type}*_{subtype}_{tag}_fcounts.csv.xz",
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*hisat*_{species}/{species}_*{type}*_{subtype}_{tag}.count.xz"),
    "deduplication_stats" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/deduplication_stats.txt"),
    "freebayes_variants_by_gene" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/variants_by_gene*"),
    "freebayes_variants_table" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/all_tags*"),
    "freebayes_modified_genome" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/{species}*.fasta"),
    "freebayes_bcf_file" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/{species}*.bcf"),
    "freebayes_penetrance_file" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*freebayes_{species}/variants_penetrance*"),
    "bedtools_coverage_file" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*bedtools_coverage_{species}/*.bed"),
    "bbmap_coverage_stats" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*bam2coverage*/coverage.tsv.xz"),
    "bbmap_coverage_per_nt" = list(
      "file" = "{basedir}/{meta[['sampleid']]}/outputs/*bam2coverage*/base_coverage.tsv.xz")
  )
  return(specification)
}

#' Steal transcript IDs from the first count table.
#'
#' @param meta Input metadata containing the salmon count table names.
#' @param annotations Extant set of gene annotations, likely from biomart.
#' @param meta_column metadata column with the filenames.
#' @param annot_gene_column Column of annotations with the gene IDs.
#' @param annot_tx_column Column of annotations with the transcript IDs.
#' @param keep_unique Drop the potential duplicate GIDs?
#' @return List containing modified annotations for the genes, transcripts,
#'  and the map between them.
#' @export
steal_salmon_tx_ids <- function(meta, annotations, meta_column = "salmon_count_table",
                                annot_gene_column = "ensembl_gene_id",
                                annot_tx_column = "ensembl_transcript_id",
                                keep_unique = TRUE) {
  ## FIXME: Use dispatch for this
  if ("preprocessing_metadata" %in% class(meta)) {
    meta <- meta[["new_meta"]]
  } else if ("character" %in% class(meta)) {
    meta <- extract_metadata(meta)
  }
  salmon_files <- meta[[meta_column]]
  first_file <- salmon_files[1]
  stolen_versions <- readr::read_tsv(first_file, skip = 1,
                                     col_names = c("Name", "Length", "EffectiveLength", "TPM", "NumReads"), col_types = c("cdddd"))
  stolen_versions[[annot_tx_column]] <- gsub(pattern = "^(.+)\\.\\d+$",
                                             replacement = "\\1",
                                             x = stolen_versions[["Name"]])
  stolen_versions <- as.data.frame(stolen_versions[, c(annot_tx_column, "Name")])
  colnames(stolen_versions) <- c(annot_tx_column, "salmon_tx_version")
  rownames(stolen_versions) <- stolen_versions[["salmon_tx_version"]]
  ## Add the new salmon_tx_version to the annotations.
  merged_annotations <- merge(annotations, stolen_versions, by = annot_tx_column, all.y = TRUE)

  ## Grab the tx_gene_map and set its column names.
  tx_gene_map <- merged_annotations[, c("salmon_tx_version", annot_gene_column)]
  colnames(tx_gene_map) <- c("transcript", "gene")
  if (isTRUE(keep_unique)) {
    kept <- !duplicated(tx_gene_map[["transcript"]])
    tx_gene_map <- tx_gene_map[kept, ]
  }
  na_genes <- is.na(tx_gene_map[["gene"]])
  if (sum(na_genes) > 0) {
    warning(sum(na_genes), " genes failed to match transcripts in the data, setting them to the txid.")
    tx_gene_map[na_genes, "gene"] <- tx_gene_map[na_genes, "transcript"]
  }
  rownames(tx_gene_map) <- make.names(tx_gene_map[["transcript"]], unique = TRUE)

  ## Set the rownames of the tx_annotations to the salmon ids
  tx_annotations <- merged_annotations
  rownames(tx_annotations) <- make.names(tx_annotations[["salmon_tx_version"]], unique = TRUE)

  ## Finally, set the rownames of the gene annotations.
  gene_annotations <- merged_annotations
  if (isTRUE(keep_unique)) {
    kept <- !duplicated(gene_annotations[[annot_gene_column]])
    gene_annotations <- gene_annotations[kept, ]
  }
  rownames(gene_annotations) <- make.names(gene_annotations[[annot_gene_column]], unique = TRUE)

  retlist <- list(
    "tx_annotations" = tx_annotations,
    "gene_annotations" = gene_annotations,
    "tx_gene_map" = tx_gene_map)
  return(retlist)
}

## EOF
