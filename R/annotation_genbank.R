## annotation_genbank.r: A group of functions to extract annotations from
## GenBank. Probably the second most likely annotation source is NCBI, so this
## file seeks to simplify extracting annotations from genbank flat files and/or
## the NCBI web interface.

#' Given a genbank accession, make a txDb object along with sequences, etc.
#'
#' Let us admit it, sometimes biomart is a pain.  It also does not have easily
#' accessible data for microbes.  Genbank does!
#'
#' Tested in test_40ann_biomartgenbank.R and test_70expt_spyogenes.R
#' This primarily sets some defaults for the genbankr service in order to
#' facilitate downloading genomes from genbank and dumping them into a local
#' txdb instance.
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
load_genbank_annotations <- function(accession = "AE009949", db = "nucleotide",
                                     file = NULL, sequence = TRUE, type = "CDS",
                                     savetxdb = FALSE, restez_db = "/sw/local/genbank/current") {
  restez::restez_path_set(filepath = restez_db)
  all_ids <- restez::list_db_ids(n = NULL)
  res <- list()
  if (accession %in% all_ids) {
    res <- restez::gb_record_get(accession)
  } else {
    res <- restez::entrez_fetch(db = db, id = accession, rettype = db)
  }
  data <- restez::gb_extract(record = res, what = "features")
  ret <- data.frame()
  hits <- 0
  for (i in seq_along(data)) {
    element <- data[[i]]
    if (element[["type"]] == type) {
      hits <- hits + 1
      row <- as.data.frame(element)
      if (hits == 1) {
        ret <- row
      } else {
        ret <- as.data.frame(data.table::rbindlist(list(ret, row), fill = TRUE))
      }
    } else {
      next
    }
  }
  if (!is.null(ret[["location"]]) && is.null(ret[["start"]])) {
    ret[["start"]] <- gsub(x = ret[["location"]], pattern = "^([[:digit:]]+).*$",
                           replacement = "\\1")
    ret[["end"]] <- gsub(x = ret[["location"]], pattern = "^.*\\.\\.([[:digit:]]+)$",
                         replacement = "\\1")
  }
  return(ret)
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
