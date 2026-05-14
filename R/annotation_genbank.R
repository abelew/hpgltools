## annotation_genbank.R: A group of functions to extract annotations from
## GenBank. Probably the second most likely annotation source is NCBI, so this
## file seeks to simplify extracting annotations from genbank flat files and/or
## the NCBI web interface.

#' @include 01_hpgltools.R
NULL

#' Given a genbank accession, make a txDb object along with sequences, etc.
#'
#' Let us admit it, sometimes biomart is a pain.  It also does not have easily
#' accessible data for microbes.  Genbank does!
#'
#' This switched to using restez, which does not provide the helpful
#' data structures that genbankr did.  I am not sure if I can still
#' create orgdb/txdb entries therefore unless I yank that code out of
#' an old version of genbankr, or if I make use of txdbmaker etc.
#'
#' @param accession Accession to download and import.
#' @param db Entrez Database to query
#' @param file Use a file instead of downloading the accession?
#' @param sequence Download the sequence with the annotations?
#' @param type Extract entries of this type.
#' @param savetxdb Attempt saving a txdb object?
#' @param restez_db Location of local genbank used by restez.
#' @return List containing a txDb, sequences, and some other stuff which I
#'  haven't yet finalized.
#' @seealso [Biostrings] [GenomicFeatures] [genbankr::import()] [genbankr::readGenBank()]
#' @example inst/examples/annotation_genbank.R
#' @export
load_genbank_annotations <- function(accession = "AE009949", database = "nuccore",
                                     return_type = "gbwithparts", save_gb = TRUE,
                                     restez_db = "/sw/local/genbank/current") {
  restez_found <- try(restez::restez_path_set(filepath = restez_db), silent = TRUE)
  if ("try-error" %in% class(restez_found)) {
    message("The restez path: ", restez_db, " failed to initialized.")
    all_ids <- NULL
  } else {
    message("Collecting IDs from the local restez database at: ", restez_db, ".")
    all_ids <- restez::list_db_ids(n = NULL)
  }
  res <- list()
  if (accession %in% all_ids) {
    message("Collecting accession ", accession, " from the local database.")
    res <- restez::gb_record_get(accession)
  } else {
    message("Fetching the accession: ", accession, " from the ncbi ",
            database, " database.")
    res <- restez::entrez_fetch(db = database, id = accession, rettype = return_type)
    message("Downloaded ", stringr::str_length(res), " bytes.")
    if (isTRUE(save_gb)) {
      file_connection <- file(paste0(accession, ".gb"))
      written <- writeLines(res, file_connection)
      closed <- close(file_connection)
    }
  }
  gb_features <- restez::gb_extract(record = res, what = "features")
  message("Extracted ", length(gb_features), " from the genbank data.")
  gb_sequence <- restez::gb_extract(record = res, what = "sequence")
  message("Extracted ", stringr::str_length(gb_sequence),
          " characters from the gb sequence.")
  gb_definition <- restez::gb_extract(record = res, what = "definition")
  gb_version <- restez::gb_extract(record = res, what = "version")
  gb_biostring <- Biostrings::DNAStringSet(gb_sequence)

  feature_list <- list()
  feature_hits <- list()
  for (i in seq_along(gb_features)) {
    element <- gb_features[[i]]
    this_type <- element[["type"]]
    if (is.null(feature_list[[this_type]])) {
      feature_list[[this_type]] <- data.frame()
      feature_hits[[this_type]] <- 0
    }
    row <- as.data.frame(element)
    feature_hits[[this_type]] <- feature_hits[[this_type]] + 1
    if (feature_hits[[this_type]] == 1) {
      this_df <- row
      feature_list[[this_type]] <- this_df
    } else {
      feature_list[[this_type]] <- as.data.frame(
        data.table::rbindlist(list(feature_list[[this_type]], row), fill = TRUE))
    }
  }

  gr_gene_df <- feature_list[["gene"]]
  gr_gene_df[["strand"]] <- "+"
  neg_features <- grepl(x = gr_gene_df[["location"]], pattern = "complement")
  gr_gene_df[neg_features, "strand"] <- "-"
  gr_gene_df[["start"]] <- gsub(x = gr_gene_df[["location"]], pattern = "^complement\\(", replacement = "")
  gr_gene_df[["start"]] <- gsub(x = gr_gene_df[["start"]], pattern = "\\)$", replacement = "")
  gr_gene_df[["end"]] <- as.numeric(gsub(x = gr_gene_df[["start"]], pattern = "^([[:digit:]]+)\\.\\.([[:digit:]]+)$", replacement = "\\2"))
  gr_gene_df[["start"]] <- as.numeric(gsub(x = gr_gene_df[["start"]], pattern = "^([[:digit:]]+)?\\.\\.([[:digit:]]+)$", replacement = "\\1"))

  gr_cds_df <- feature_list[["CDS"]]
  gr_cds_df[["strand"]] <- "+"
  neg_features <- grepl(x = gr_cds_df[["location"]], pattern = "complement")
  gr_cds_df[neg_features, "strand"] <- "-"
  gr_cds_df[["start"]] <- gsub(x = gr_cds_df[["location"]], pattern = "^complement\\(", replacement = "")
  gr_cds_df[["start"]] <- gsub(x = gr_cds_df[["start"]], pattern = "\\)$", replacement = "")
  gr_cds_df[["end"]] <- as.numeric(gsub(x = gr_cds_df[["start"]], pattern = "^([[:digit:]]+)\\.\\.([[:digit:]]+)$", replacement = "\\2"))
  gr_cds_df[["start"]] <- as.numeric(gsub(x = gr_cds_df[["start"]], pattern = "^([[:digit:]]+)?\\.\\.([[:digit:]]+)$", replacement = "\\1"))

  gr_gene_df[["chr"]] <- gb_definition
  gr_cds_df[["chr"]] <- gb_definition

  gene_gr <- GenomicRanges::makeGRangesFromDataFrame(
    gr_gene_df, keep.extra.columns = TRUE)
  cds_gr <- GenomicRanges::makeGRangesFromDataFrame(
    gr_cds_df, keep.extra.columns = TRUE)

  ## Hmm it looks like restez mis-parses some genbank features.  That is super-annoying.
  ## It also does not keep the chromosome ID with the features, which
  ## will make creating a GRange difficult.
  ## I will therefore set this aside for now and return a data
  ## structure containing the outputs from restez and my somewhat
  ## mangled list of dataframes.

  retlist <- list(
    "definition" = gb_definition,
    "version" = gb_version,
    "sequence_vector" = gb_sequence,
    "sequence" = gb_biostring,
    "features" = gb_features,
    "feature_list" = feature_list,
    "gene_granges" = gene_gr,
    "cds_granges" = cds_gr)
  class(retlist) <- c("hpgltools::load_genbank_annotations", "list")
  return(retlist)
}

#' A genbank accession downloader scurrilously stolen from ape.
#'
#' This takes and downloads genbank accessions.
#'
#' Tested in test_40ann_biomartgenbank.R
#' In this function I stole the same functionality from the ape package and set
#' a few defaults so that it hopefully fails less often.
#'
#' @param accessions An accession -- actually a set of them.
#' @param write Write the files?  Otherwise return a list of the strings
#' @return A list containing the number of files downloaded and the character
#'   strings acquired.
#' @seealso [ape]
#' @examples
#'  written <- download_gbk(accessions = "AE009949")
#'  written$written_file
#' @author The ape authors with some modifications by atb.
#' @export
download_gbk <- function(accessions = "AE009949", write = TRUE) {
  num_acc <- length(accessions)
  nrequest <- num_acc %/% 400 + as.logical(num_acc %% 400)
  downloaded <- character(0)
  num_downloaded <- 0
  strings <- list()
  written_file <- NULL
  for (i in seq_len(nrequest)) {
    a <- (i - 1) * 400 + 1
    b <- 400 * i
    if (i == nrequest) {
      b <- num_acc
    }
    accession <- accessions[i]

    url <- paste0(
        "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=",
        paste(accessions[a:b], collapse = ","), "&rettype=gb&retmode=text&report=gbwithparts")

    dl_file <- glue("{accession}.gb")
    data <- try(download.file(url = url, destfile = dl_file, quiet = TRUE))
    scanned <- NULL
    if (class(data) != "try-error") {
      scanned <- try(scan(file = dl_file, what = "", sep = "\n", quiet = TRUE))
    }
    if (class(scanned) != "try-error") {
      downloaded <- c(downloaded, scanned)
      num_downloaded <- num_downloaded + 1
    }
    strings[[accession]] <- downloaded
    written_file <- glue("{accessions[a]}.gb")
    if (isTRUE(write)) {
      file_connection <- file(written_file)
      writeLines(downloaded, file_connection)
      close(file_connection)
    }
  } ## End of for loop
  retlist <- list(
      "written_file" = written_file,
      "num_success" = num_downloaded,
      "strings" = strings)
  return(retlist)
}

## EOF
