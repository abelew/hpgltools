#' Count n-mers in a given data set using Biostrings
#'
#' This just calls PDict() and vcountPDict() on a sequence database given a
#' pattern and number of mismatches.  This may be used by divide_seq()
#' normalization.
#'
#' @param genome Sequence database, genome in this case.
#' @param pattern Count off this string.
#' @param mismatch How many mismatches are acceptable?
#' @return Set of counts by sequence.
#' @export
count_nmer <- function(genome, pattern = "ATG", mismatch = 0) {
  seq_obj <- Biostrings::getSeq(genome)
  dict <- Biostrings::PDict(pattern, max.mismatch = mismatch)
  result <- as.data.frame(Biostrings::vcountPDict(dict, seq_obj))
  rownames(result) <- pattern
  colnames(result) <- names(seq_obj)
  result <- t(result)
  return(result)
}
setGeneric("count_nmer")

#' Given an eupathdb species lacking UTR boundaries, extract an arbitrary region
#' before/after each gene.
#'
#' This is a very domain-specific function.
#'
#' @param species_name Species name for which to query the eupathdb.
#' @param entry EuPathDB metadatum entry.
#' @param webservice If specified, makes the query faster, I always used
#'  tritrypdb.org.
#' @param padding Number of nucleotides to gather.
#' @param ... Extra arguments for the various EuPathDB functions.
#' @return Set of padding UTR sequences/coordinates.
gather_eupath_utrs_padding <- function(species_name = "Leishmania major", entry = NULL,
                                       webservice = "tritrypdb", padding = 200, ...) {
  if (!is.null(entry)) {
    pkg_names <- EuPathDB::get_eupath_pkgnames(entry)
  } else {
    entry <- EuPathDB::get_eupath_entry(species_name, webservice = webservice)
    pkg_names <- EuPathDB::get_eupath_pkgnames(entry)
  }
  bsgenome_name <- pkg_names[["bsgenome"]]
  orgdb_name <- pkg_names[["orgdb"]]
  if (!isTRUE(pkg_names[["bsgenome_installed"]])) {
    genome_installedp <- EuPathDB::make_eupath_bsgenome(species = species_name, entry = entry,
                                                        ...)
  }
  if (!isTRUE(pkg_names[["orgdb_installed"]])) {
    orgdb_installedp <- EuPathDB::make_eupath_orgdb(entry = entry,
                                                    ...)
  }

  ##lib_result <- sm(library(orgdb_name, character.only = TRUE))
  lib_result <- sm(requireNamespace(orgdb_name))
  att_result <- sm(try(attachNamespace(orgdb_name), silent = TRUE))
  ##lib_result <- sm(library(bsgenome_name, character.only = TRUE))
  lib_result <- sm(requireNamespace(bsgenome_name))
  att_result <- sm(try(attachNamespace(bsgenome_name), silent = TRUE))
  orgdb <- get0(orgdb_name)
  bsgenome <- get0(bsgenome_name)
  wanted_fields <- c("annot_gene_location_text", "annot_gene_name",
                     "annot_gene_product", "annot_gene_type")
  annot <- load_orgdb_annotations(orgdb, keytype = "gid", fields = wanted_fields)
  annot <- EuPathDB::extract_gene_locations(annot[["genes"]])
  annot_complete <- complete.cases(annot)
  annot_df <- annot[annot_complete, ]
  annot_df[["chr_length"]] <- 0

  ## Gather the chromosome lengths to make sure we don't pass them.
  genome_info <- BiocGenerics::as.data.frame(BSgenome::seqinfo(bsgenome))
  chr_names <- rownames(genome_info)
  for (l in seq_along(chr_names)) {
    name <- chr_names[l]
    len <- genome_info[l, "seqlengths"]
    chr_idx <- annot_df[["chromosome"]] == name
    annot_df[chr_idx, "chr_length"] <- len
  }

  result <- gather_utrs_padding(bsgenome, annot_df, padding = padding,
                                ...)
  return(result)
}

#' Take a BSgenome and data frame of chr/start/end/strand, provide 5' and 3'
#' padded sequence.
#'
#' For some species, we do not have a fully realized set of UTR boundaries, so
#' it can be useful to query some arbitrary and consistent amount of sequence
#' before/after every CDS sequence.  This function can provide that information.
#' Note, I decided to use tibble for this so that if one accidently prints too
#' much it will not freak out.
#'
#' @param bsgenome BSgenome object containing the genome of interest.
#' @param annot_df Annotation data frame containing all the entries of
#'  interest, this is generally extracted using a function in the
#'  load_something_annotations() family (load_orgdb_annotations() being the most
#'  likely).
#' @param gid Specific GID(s) to query.
#' @param name_column Give each gene a name using this column.
#' @param chr_column Column name of the chromosome names.
#' @param start_column Column name of the start information.
#' @param end_column Ibid, end column.
#' @param strand_column Ibid, strand.
#' @param type_column Subset the annotation data using this column, if not null.
#' @param gene_type Subset the annotation data using the type_column with this
#'  type.
#' @param padding Return this number of nucleotides for each gene.
#' @param ... Arguments passed to child functions (I think none currently).
#' @return Dataframe of UTR, CDS, and UTR+CDS sequences.
#' @export
gather_utrs_padding <- function(bsgenome, annot_df, gid = NULL, name_column = "gid",
                                chr_column = "chromosome", start_column = "start",
                                end_column = "end", strand_column = "strand",
                                type_column = "annot_gene_type",
                                gene_type = "protein coding", padding = c(80, 160), ...) {
  arglist <- list(...)
  retlist <- list()
  if (!is.null(type_column)) {
    ## Pull the entries which are of the desired type.
    all_gene_idx <- annot_df[[type_column]] == gene_type
    annot_df <- annot_df[all_gene_idx, ]
  }
  if (!is.null(gid)) {
    annot_df <- annot_df[gid, ]
  }
  na_idx <- is.na(annot_df[[strand_column]])
  if (sum(na_idx) > 0) {
    warning("Dropping ", sum(na_idx), " rows.")
    annot_df <- annot_df[!na_idx, ]
  }
  plusone_idx <- annot_df[[strand_column]] == 1
  if (sum(one_idx) > 0) {
    annot_df[plusone_idx, strand_column] <- "+"
  }
  minusone_idx <- annot_df[[strand_column]] == -1
  if (sum(minusone_idx) > 0) {
    annot_df[minusone_idx, strand_column] <- "-"
  }
  annot_df[[start_column]] <- as.numeric(annot_df[[start_column]])
  annot_df[[end_column]] <- as.numeric(annot_df[[end_column]])

  if (length(padding) == 2) {
    fivep_padding <- padding[1]
    threep_padding <- padding[2]
    padding <- max(fivep_padding, threep_padding)
    ### (low is start - fivep) 5'--->3'
    annot_df[plusone_idx, "low_boundary"] <- annot_df[plusone_idx, start_column] - fivep_padding
    ## 3'<---5' (high is end + fivep)
    annot_df[minusone_idx, "high_boundary"] <- annot_df[minusone_idx, end_column] + fivep_padding
    ## 5'--->3' (high is end + threep)
    annot_df[plusone_idx, "high_boundary"] <- annot_df[plusone_idx, end_column] + threep_padding
    ## (low is start - threep) 3'<---5'
    annot_df[minusone_idx, "low_boundary"] <- annot_df[minusone_idx, start_column] - threep_padding
  } else if (length(padding) == 1) {
    fivep_paddig <- padding
    threep_padding <- padding
    ## Here is the thing I absolutely cannot remember:
    ## start is _always_ a lower number than end.
    annot_df[["low_boundary"]] <- annot_df[[start_column]] - padding
    annot_df[["high_boundary"]] <- annot_df[[end_column]] + padding
  } else {
    stop("I only understand padding with 1 or 2 elements.")
  }

  ## Check for NAs again
  na_idx <- is.na(annot_df[["low_boundary"]])
  if (sum(na_idx) > 0) {
    warning("Dropping ", sum(na_idx), " rows.")
    annot_df <- annot_df[!na_idx, ]
  }

  ## Add some logic to make sure we do not have sequences which
  ## expand beyond the chromosome boundaries.
  sub_zero_idx <- annot_df[["low_boundary"]] < 1
  if (sum(sub_zero_idx) > 0) {
    message("Found ", sum(sub_zero_idx), " genes which are less than ",
            padding, " nt. from the beginning of the chromosome.")
    annot_df[sub_zero_idx, "low_boundary"] <- 1
  }
  past_end_idx <- annot_df[["high_boundary"]] > annot_df[["chr_length"]]
  if (sum(past_end_idx) > 0) {
    message("Found ", sum(past_end_idx), " genes which are less than ",
            padding, " nt. from the end of the chromosome.")
    annot_df[past_end_idx, "high_boundary"] <- annot_df[past_end_idx, "chr_length"]
  }

  plus_idx <- annot_df[[strand_column]] == "+"
  minus_idx <- annot_df[[strand_column]] == "-"
  plus_entries <- annot_df[plus_idx, ]
  minus_entries <- annot_df[minus_idx, ]

  pluses_fivep <- GenomicRanges::GRanges(
                                     seqnames = S4Vectors::Rle(plus_entries[[chr_column]]),
                                     ranges = IRanges::IRanges(
                                                           start = plus_entries[["low_boundary"]],
                                                           end = plus_entries[[start_column]]),
                                     strand = S4Vectors::Rle(plus_entries[[strand_column]]),
                                     name = S4Vectors::Rle(plus_entries[[name_column]]))
  range_names <- as.character(GenomicRanges::seqnames(pluses_fivep))
  genome_names <- as.character(names(bsgenome))
  keep_idx <- range_names %in% genome_names
  pluses_fivep <- pluses_fivep[keep_idx]
  retlist[["fiveprime_plus_granges"]] <- pluses_fivep

  pluses_threep <- GenomicRanges::GRanges(
                                      seqnames = S4Vectors::Rle(plus_entries[[chr_column]]),
                                      ranges = IRanges::IRanges(
                                                            start = plus_entries[[end_column]],
                                                            end = plus_entries[["high_boundary"]]),
                                      strand = S4Vectors::Rle(plus_entries[[strand_column]]),
                                      name = S4Vectors::Rle(plus_entries[[name_column]]))
  range_names <- as.character(GenomicRanges::seqnames(pluses_threep))
  genome_names <- as.character(names(bsgenome))
  keep_idx <- range_names %in% genome_names
  pluses_threep <- pluses_threep[keep_idx]
  retlist[["threeprime_plus_granges"]] <- pluses_threep

  pluses_cds <- GenomicRanges::GRanges(
                                   seqnames = S4Vectors::Rle(plus_entries[[chr_column]]),
                                   ranges = IRanges::IRanges(
                                                         start = plus_entries[[start_column]],
                                                         end = plus_entries[[end_column]]),
                                   strand = S4Vectors::Rle(plus_entries[[strand_column]]),
                                   name = S4Vectors::Rle(plus_entries[[name_column]]))
  range_names <- as.character(GenomicRanges::seqnames(pluses_cds))
  genome_names <- as.character(names(bsgenome))
  keep_idx <- range_names %in% genome_names
  pluses_cds <- pluses_cds[keep_idx]
  retlist[["cds_plus_granges"]] <- pluses_cds

  pluses_all <- GenomicRanges::GRanges(
                                   seqnames = S4Vectors::Rle(plus_entries[[chr_column]]),
                                   ranges = IRanges::IRanges(
                                                         start = plus_entries[["low_boundary"]],
                                                         end = plus_entries[["high_boundary"]]),
                                   strand = S4Vectors::Rle(plus_entries[[strand_column]]),
                                   name = S4Vectors::Rle(plus_entries[[name_column]]))
  range_names <- as.character(GenomicRanges::seqnames(pluses_all))
  genome_names <- as.character(names(bsgenome))
  keep_idx <- range_names %in% genome_names
  pluses_all <- pluses_all[keep_idx]

  plus_fivep_seqstrings <- BSgenome::getSeq(bsgenome, pluses_fivep)
  names(plus_fivep_seqstrings) <- pluses_fivep$name
  plus_threep_seqstrings <- BSgenome::getSeq(bsgenome, pluses_threep)
  names(plus_threep_seqstrings) <- pluses_threep$name
  plus_cds_seqstrings <- BSgenome::getSeq(bsgenome, pluses_cds)
  names(plus_cds_seqstrings) <- pluses_cds$name
  plus_all_seqstrings <- BSgenome::getSeq(bsgenome, pluses_all)
  names(plus_all_seqstrings) <- pluses_all$name

  minuses_fivep <- GenomicRanges::GRanges(
                                      seqnames = S4Vectors::Rle(minus_entries[[chr_column]]),
                                      ranges = IRanges::IRanges(
                                                            start = minus_entries[[end_column]],
                                                            end = minus_entries[["high_boundary"]]),
                                      strand = S4Vectors::Rle(minus_entries[[strand_column]]),
                                      name = S4Vectors::Rle(minus_entries[[name_column]]))
  range_names <- as.character(GenomicRanges::seqnames(minuses_fivep))
  genome_names <- as.character(names(bsgenome))
  keep_idx <- range_names %in% genome_names
  minuses_fivep <- minuses_fivep[keep_idx]
  retlist[["fiveprime_minus_granges"]] <- minuses_fivep

  minuses_threep <- GenomicRanges::GRanges(
                                       seqnames = S4Vectors::Rle(minus_entries[[chr_column]]),
                                       ranges = IRanges::IRanges(
                                                             start = minus_entries[["low_boundary"]],
                                                             end = minus_entries[[start_column]]),
                                       strand = S4Vectors::Rle(minus_entries[[strand_column]]),
                                       name = S4Vectors::Rle(minus_entries[[name_column]]))
  range_names <- as.character(GenomicRanges::seqnames(minuses_threep))
  keep_idx <- range_names %in% genome_names
  minuses_threep <- minuses_threep[keep_idx]
  retlist[["threeprime_minus_granges"]] <- minuses_threep

  minuses_cds <- GenomicRanges::GRanges(
                                    seqnames = S4Vectors::Rle(minus_entries[[chr_column]]),
                                    ranges = IRanges::IRanges(
                                                          start = minus_entries[[start_column]],
                                                          end = minus_entries[[end_column]]),
                                    strand = S4Vectors::Rle(minus_entries[[strand_column]]),
                                    name = S4Vectors::Rle(minus_entries[[name_column]]))
  range_names <- as.character(GenomicRanges::seqnames(minuses_cds))
  keep_idx <- range_names %in% genome_names
  minuses_cds <- minuses_cds[keep_idx]
  retlist[["cds_minus_granges"]] <- minuses_cds

  minuses_all <- GenomicRanges::GRanges(
                                    seqnames = S4Vectors::Rle(minus_entries[[chr_column]]),
                                    ranges = IRanges::IRanges(
                                                          start = minus_entries[["low_boundary"]],
                                                          end = minus_entries[["high_boundary"]]),
                                    strand = S4Vectors::Rle(minus_entries[[strand_column]]),
                                    name = S4Vectors::Rle(minus_entries[[name_column]]))
  range_names <- as.character(GenomicRanges::seqnames(minuses_all))
  keep_idx <- range_names %in% genome_names
  minuses_all <- minuses_all[keep_idx]

  minus_fivep_seqstrings <- BSgenome::getSeq(bsgenome, minuses_fivep)
  names(minus_fivep_seqstrings) <- minuses_fivep$name
  minus_threep_seqstrings <- BSgenome::getSeq(bsgenome, minuses_threep)
  names(minus_threep_seqstrings) <- minuses_threep$name
  minus_cds_seqstrings <- BSgenome::getSeq(bsgenome, minuses_cds)
  names(minus_cds_seqstrings) <- minuses_cds$name
  minus_all_seqstrings <- BSgenome::getSeq(bsgenome, minuses_all)
  names(minus_all_seqstrings) <- minuses_all$name

  ## These provide data frames of the sequence lexically before/after every gene.
  plus_fivep_tbl <- tibble::as_tibble(as.data.frame(plus_fivep_seqstrings))
  plus_fivep_tbl[["name"]] <- names(plus_fivep_seqstrings)
  colnames(plus_fivep_tbl) <- c("sequence", "ID")
  plus_fivep_tbl <- plus_fivep_tbl[, c("ID", "sequence")]
  retlist[["fiveprime_plus_table"]] <- plus_fivep_tbl

  plus_threep_tbl <- tibble::as_tibble(as.data.frame(plus_threep_seqstrings))
  plus_threep_tbl[["name"]] <- names(plus_threep_seqstrings)
  colnames(plus_threep_tbl) <- c("sequence", "ID")
  plus_threep_tbl <- plus_threep_tbl[, c("ID", "sequence")]
  retlist[["threeprime_plus_table"]] <- plus_threep_tbl

  minus_fivep_tbl <- tibble::as_tibble(as.data.frame(minus_fivep_seqstrings))
  minus_fivep_tbl[["name"]] <- names(minus_fivep_seqstrings)
  colnames(minus_fivep_tbl) <- c("sequence", "ID")
  minus_fivep_tbl <- minus_fivep_tbl[, c("ID", "sequence")]
  retlist[["fiveprime_minus_table"]] <- minus_fivep_tbl

  minus_threep_tbl <- tibble::as_tibble(as.data.frame(minus_threep_seqstrings))
  minus_threep_tbl[["name"]] <- names(minus_threep_seqstrings)
  colnames(minus_threep_tbl) <- c("sequence", "ID")
  minus_threep_tbl <- minus_threep_tbl[, c("ID", "sequence")]
  retlist[["threeprime_minus_table"]] <- minus_threep_tbl

  plus_cds_tbl <- tibble::as_tibble(as.data.frame(plus_cds_seqstrings))
  plus_cds_tbl[["name"]] <- names(plus_cds_seqstrings)
  colnames(plus_cds_tbl) <- c("sequence", "ID")
  plus_cds_tbl <- plus_cds_tbl[, c("ID", "sequence")]
  retlist[["cds_plus_table"]] <- plus_cds_tbl

  minus_cds_tbl <- tibble::as_tibble(as.data.frame(minus_cds_seqstrings))
  minus_cds_tbl[["name"]] <- names(minus_cds_seqstrings)
  colnames(minus_cds_tbl) <- c("sequence", "ID")
  minus_cds_tbl <- minus_cds_tbl[, c("ID", "sequence")]
  retlist[["cds_minus_table"]] <- minus_cds_tbl

  return(retlist)
}

#' Get UTR sequences using information provided by TxDb and fiveUTRsByTranscript
#'
#' For species like Mus musculus, load_orgdb_annotations(Mus.musculus)
#' should return a list including the requisite GRanges for the 5'/3' UTRs.
#'
#' @param bsgenome A BSGenome instance containing the encoded genome.
#' @param fivep_utr Locations of the 5' UTRs.
#' @param threep_utr Locations of the 3' UTRs.
#' @param start_column What column in the annotation data contains the starts?
#' @param strand_column What column in the annotation data contains the sequence strands?
#' @param end_column Column in the data with the end locations.
#' @param chr_column Column in the df with the chromosome names.
#' @param name_column Finally, where are the gene names?
#' @param ... Parameters passed to child functions.
#' @return UTRs!
gather_utrs_txdb <- function(bsgenome, fivep_utr = NULL, threep_utr = NULL,
                             start_column = "start", end_column = "end",
                             strand_column = "strand", chr_column = "seqnames",
                             name_column = "group_name", ...) {
  fivep_df <- as.data.frame(fivep_utr)
  threep_df <- as.data.frame(threep_utr)


  ## These are named 'smaller' and 'larger' because this is extracting the ranges
  ## which are numerically smaller/larger than the CDS entry, not 5'/3' because
  ## this does not yet take into account the strand information.
  fivep_ranges <- GenomicRanges::GRanges(
                                     seqnames = S4Vectors::Rle(fivep_df[[chr_column]]),
                                     ranges = IRanges::IRanges(
                                                           start = fivep_df[[start_column]],
                                                           end = fivep_df[[end_column]] + 2),
                                     strand = S4Vectors::Rle(fivep_df[[strand_column]]),
                                     name = S4Vectors::Rle(fivep_df[[name_column]]))
  threep_ranges <- GenomicRanges::GRanges(
                                      seqnames = S4Vectors::Rle(threep_df[[chr_column]]),
                                      ranges = IRanges::IRanges(
                                                            start = threep_df[[start_column]],
                                                            end = threep_df[[end_column]] + 2),
                                      strand = S4Vectors::Rle(threep_df[[strand_column]]),
                                      name = S4Vectors::Rle(threep_df[[name_column]]))
  fivep_seqstrings <- BSgenome::getSeq(bsgenome, fivep_ranges)
  threep_seqstrings <- BSgenome::getSeq(bsgenome, threep_ranges)
  names(fivep_seqstrings) <- fivep_ranges$name
  names(threep_seqstrings) <- threep_ranges$name

  ## These provide data frames of the sequence lexically before/after every gene.
  fivep_seqdt <- data.table::as.data.table(fivep_seqstrings)
  fivep_seqdt[["name"]] <- as.character(fivep_ranges$name)
  colnames(fivep_seqdt) <- c("sequence", "name")
  threep_seqdt <- data.table::as.data.table(threep_seqstrings)
  threep_seqdt[["name"]] <- as.character(threep_ranges$name)
  colnames(threep_seqdf) <- c("sequence", "name")
  ## These provide the strand and location information for the same.
  fivep_infodt <- data.table::as.data.table(fivep_ranges)
  threep_infodt <- data.table::as.data.table(threep_ranges)
  ## Bring together the location and sequence information.
  all_fivep <- cbind(fivep_infodt, fivep_seqdt)
  all_threep <- cbind(threep_infodt, threep_seqdt)

  retlist <- list(
      "five_prime" = all_fivep,
      "three_prime" = all_threep)
  return(retlist)
}

#' Gather some simple sequence attributes.
#'
#' This extends the logic of the pattern searching in pattern_count_genome() to
#' search on some other attributes.
#'
#' @param fasta Genome encoded as a fasta file.
#' @param gff Optional gff of annotations (if not provided it will just ask the
#'  whole genome).
#' @param type Column of the gff file to use.
#' @param key What type of entry of the gff file to key from?
#' @return List of data frames containing gc/at/gt/ac contents.
#' @seealso [Biostrings] [Rsamtools]
#' @examples
#'  pa_data <- get_paeruginosa_data()
#'  pa_fasta <- pa_data[["fasta"]]
#'  pa_gff <- pa_data[["gff"]]
#'  pa_attribs <- sequence_attributes(pa_fasta, gff = pa_gff)
#'  head(pa_attribs)
#' @export
sequence_attributes <- function(fasta, gff = NULL, type = "gene", key = NULL) {
  rawseq <- Rsamtools::FaFile(fasta)
  if (is.null(key)) {
    key <- c("ID", "locus_tag")
  }
  if (is.null(gff)) {
    entry_sequences <- Biostrings::getSeq(rawseq)
  } else {
    ## entries <- rtracklayer::import.gff3(gff, asRangedData = FALSE)
    entries <- rtracklayer::import.gff(gff)
    ## keep only the ones of the type of interest (gene).
    type_entries <- entries[entries$type == type, ]
    ## Set some hopefully sensible names.
    names(type_entries) <- rownames(type_entries)
    ## Get the sequence from the genome for them.
    entry_sequences <- Biostrings::getSeq(rawseq, type_entries)
    tmp_entries <- as.data.frame(type_entries)
    ## Give them some sensible names
    for (k in key) {
      if (!is.null(tmp_entries[[k]])) {
        names(entry_sequences) <- tmp_entries[[k]]
      }
    }
  }
  attribs <- data.frame(
      "gc" = Biostrings::letterFrequency(entry_sequences, "CG", as.prob = TRUE),
      "at" = Biostrings::letterFrequency(entry_sequences, "AT", as.prob = TRUE),
      "gt" = Biostrings::letterFrequency(entry_sequences, "GT", as.prob = TRUE),
      "ac" = Biostrings::letterFrequency(entry_sequences, "AC", as.prob = TRUE),
      stringsAsFactors = FALSE)
  rownames(attribs) <- names(entry_sequences)
  colnames(attribs) <- c("gc", "at", "gt", "ac")
  return(attribs)
}

make_bsgenome_from_fasta <- function(fasta, pkgname, title, organism, common_name, provider, build_dir = "build", genome_prefix = NULL, author = "Ashton Trey Belew <abelew@gmail.com>", URL = "https://github.com/abelew/hpgltools", installp = TRUE) {
  if (!file.exists(build_dir)) {
    created <- dir.create(build_dir, recursive = TRUE)
  }
  input <- Biostrings::readDNAStringSet(fasta)
  output_list <- list()
  sequence_names <- "c("
  message(" Writing chromosome files, this is slow for fragmented scaffolds.")
  show_progress <- interactive() && is.null(getOption("knitr.in.progress"))
  #if (isTRUE(show_progress)) {
  #  bar <- utils::txtProgressBar(style = 3)
  #}
  genome_prefix <- NULL
  for (index in seq_len(length(input))) {
   # if (isTRUE(show_progress)) {
   #   pct_done <- index / length(input)
   #   setTxtProgressBar(bar, pct_done)
   # }
    chr <- names(input)[index]
    chr_name <- strsplit(chr, split = " ")[[1]][1]
    if (is.null(genome_prefix)) {
      genome_prefix <- gsub(x = chr_name, pattern = "(.*)\\.(.*)", replacement = "\\1")
    }
    ## chr_file <- file.path(bsgenome_dir, paste0(chr_name, ".fa"))
    chr_file <- file.path(build_dir, glue::glue("{chr_name}.fa"))
    output <- Biostrings::writeXStringSet(input[index], chr_file, append = FALSE,
                                          compress = FALSE, format = "fasta")
    output_list[[chr_name]] <- chr_file
    sequence_names <- paste0(sequence_names, '"', chr_name, '", ')
  }
  #if (isTRUE(show_progress)) {
  #  close(bar)
  #}
  message("Finished writing ", length(input), " contigs.")
  sequence_names <- gsub(pattern = ", $", replacement = ")", x = sequence_names)

  ## Now start creating the DESCRIPTION file
  desc_file <- file.path(build_dir, "DESCRIPTION")
  descript <- desc::description$new("!new")
  descript$set(Package = pkgname)
  descript$set(Title = title)
  descript$set(Author = author)
  version_string <- format(Sys.time(), "%Y.%m")
  descript$set(Version = version_string)
  descript$set(Maintainer = author)
  descript$set(Description = glue::glue("A BSgenome for {title}."))
  descript$set(License = "Artistic-2.0")
  descript$set(URL = "https://eupathdb.org")
  descript$set(seqs_srcdir = build_dir)
  descript$set(seqnames = sequence_names)
  descript$set(organism = organism)
  descript$set(common_name = common_name)
  descript$set(provider = provider)
  descript$set(release_date = format(Sys.time(), "%Y%m%d"))
  descript$set(BSgenomeObjname = common_name)
  ## descript$set(provider_version = glue::glue("{fasta_hostname} {db_version}"))
  ## descript$set(release_name = as.character(db_version))
  ## descript$set(release_name = glue::glue("{fasta_hostname} {db_version}"))
  ## descript$set(genome = genome_prefix)
  descript$set(genome = "hg38")  ## There is some weird BS with this field, so I am faking it.
  descript$set(organism_biocview = common_name)
  descript$del("LazyData")
  descript$del("Authors@R")
  descript$del("URL")
  descript$del("BugReports")
  descript$del("Encoding")
  description_file <- file.path(build_dir, "DESCRIPTION")
  descript$write(description_file)

  ## Generate the package, this puts it into the cwd.
  message(" Calling forgeBSgenomeDataPkg().")
  ## Otherwise I get error in cannot find uniqueLetters (this seems to be a new development)
  ## Invoking library(Biostrings") annoys R CMD check, but I am not sure there is a good
  ## way around that due to limitations of Biostrings, lets see.
  uniqueLetters <- Biostrings::uniqueLetters
  tt <- try(do.call("library", as.list("Biostrings")), silent = TRUE)
  pkg_builder <- BSgenome::forgeBSgenomeDataPkg(description_file)
  if ("try-error" %in% class(pkg_builder)) {
    message("forgeBSgenomeDataPkg failed with error: ")
    message("A likely reason is too many open files, which may be changed in /etc/sysctl.conf")
    print(pkg_builder)
    pkg_builder <- NULL
    return(NULL)
  }

  built <- NULL
  workedp <- ! "try-error" %in% class(pkg_builder)
  if (isTRUE(workedp)) {
    built <- try(devtools::build(pkgname, quiet = TRUE))
    workedp <- ! "try-error" %in% class(built)
  }
  if (isTRUE(installp) && isTRUE(workedp)) {
    installed <- try(devtools::install_local(built))
  }
  retlist <- list(
    "name" = pkgname,
    "archive" = built,
    "installed" = installed,
    "contigs" = length(input))
  class(retlist) <- "built_bsgenome"
  return(retlist)
}

#' Extract CDS sequences from a genome and set of annotations.
#'
#' Given a BSGenome and some annotations, write out the CDS entries.
#'
#' @param genome BSGenome containing the raw sequence.
#' @param annot Annotation dataframe.
#' @param ids Set of annotations to write, write them all if null.
#' @param output Fasta file to write.
#' @param strand_column Column name with the strand information.
#' @param chr_column Column name with the chromosomes.
#' @param start_column Column with the start positions.
#' @param end_column Column with the end positions.
#' @param name_column Names of the CDS
#' @param name_prefix Prefix to add to the entries.
#' @export
write_cds_entries <- function(genome, annot, ids = NULL, output = "all_cds.fasta",
                              strand_column = "strand", chr_column = "chromosome",
                              start_column = "start", end_column = "end",
                              name_column = "rownames", name_prefix = "lpanamensis_mcol") {
  if (class(genome)[1] == "character") {
    genome <- Rsamtools::FaFile(genome)
  }
  seq_obj <- Biostrings::getSeq(genome)

  na_idx <- is.na(annot[[strand_column]])
  annot[na_idx, strand_column] <- "undef"
  annot[[strand_column]] <- as.character(annot[[strand_column]])
  na_idx <- is.na(annot[[start_column]])
  annot[na_idx, start_column] <- 0
  annot[[start_column]] <- as.numeric(annot[[start_column]])
  na_idx <- is.na(annot[[end_column]])
  annot[na_idx, end_column] <- 1
  annot[[end_column]] <- as.numeric(annot[[end_column]])
  na_idx <- is.na(annot[[chr_column]])
  annot[na_idx, chr_column] <- "undef"
  annot[[chr_column]] <- as.character(annot[[chr_column]])

  wanted_entries_idx <- rownames(annot) %in% ids
  message("Found ", sum(wanted_entries_idx), " in the annotation data.")
  wanted <- annot[wanted_entries_idx, ]

  if (name_column == "rownames") {
    wanted[["ID"]] <- rownames(wanted)
    name_column <- "ID"
  }
  if (!is.null(name_prefix)) {
    wanted[[name_column]] <- paste0(name_prefix, "_", wanted[[name_column]])
  }

  wanted_ranges <- GenomicRanges::GRanges(
                                      seqnames = S4Vectors::Rle(wanted[[chr_column]]),
                                      ranges = IRanges::IRanges(start = wanted[[start_column]],
                                                                end = wanted[[end_column]]),
                                      strand = S4Vectors::Rle(wanted[[strand_column]]),
                                      name = S4Vectors::Rle(wanted[[name_column]]))

  wanted_seqstrings <- Biostrings::getSeq(genome, wanted_ranges)
  names(wanted_seqstrings) <- wanted[[name_column]]

  written <- Biostrings::writeXStringSet(wanted_seqstrings, output, append = FALSE,
                                         compress = FALSE, format = "fasta")
}

## EOF
