## annotation_orgdb.r: Extract annotation from orgdb packages.  OrgDB is the
## standard SQLite based method for packaging annotation data in R.  These
## functions seek to make extracting information of interest from them easier.

#' Load organism annotation data from an orgdb sqlite package.
#'
#' Creates a dataframe gene and transcript information for a given set of gene
#' ids using the AnnotationDbi interface.
#'
#' This defaults to a few fields which I have found most useful, but the brave
#' or pathological can pass it 'all'.
#'
#' @param orgdb OrganismDb instance.
#' @param gene_ids Search for a specific set of genes?
#' @param include_go Ask the Dbi for gene ontology information?
#' @param keytype mmm the key type used?
#' @param strand_column There are a few fields I want to gather by default:
#'  start, end, strand, chromosome, type, and name; but these do not
#'  necessarily have consistent names, use this column for the chromosome
#'  strand.
#' @param start_column Use this column for the gene start.
#' @param end_column Use this column for the gene end.
#' @param chromosome_column Use this column to identify the chromosome.
#' @param type_column Use this column to identify the gene type.
#' @param name_column Use this column to identify the gene name.
#' @param fields Columns included in the output.
#' @param sum_exon_widths Perform a sum of the exons in the data set?
#' @return Table of geneids, chromosomes, descriptions, strands, types, and lengths.
#' @seealso [AnnotationDbi] [AnnotationDbi::select()] [GenomicFeatures]
#' @examples
#'  hs_orgdb_annot <- load_orgdb_annotations()
#'  summary(hs_orgdb_annot$genes)
#' @export
load_orgdb_annotations <- function(orgdb = NULL, gene_ids = NULL, include_go = FALSE,
                                   keytype = "ensembl", strand_column = "cdsstrand",
                                   start_column = "cdsstart", end_column = "cdsend",
                                   chromosome_column = "cdschrom",
                                   type_column = "gene_type", name_column = "cdsname",
                                   fields = NULL, sum_exon_widths = FALSE) {
  if (is.null(orgdb)) {
    message("Assuming Homo.sapiens.")
    org_loaded <- do.call("library", args = list("package" = "Homo.sapiens", "character.only" = TRUE))
    orgdb <- get0("Homo.sapiens")
  } else if (class(orgdb) == "character") {
    org_loaded <- do.call("library", list("package" = orgdb, "character.only" = TRUE))
    orgdb <- get0(orgdb)
  }
  keytype <- toupper(keytype)
  strand_column <- toupper(strand_column)
  start_column <- toupper(start_column)
  end_column <- toupper(end_column)
  chromosome_column <- toupper(chromosome_column)
  type_column <- toupper(type_column)
  name_column <- toupper(name_column)
  ## Caveat: if fields was NULL, now it is character(0)
  fields <- toupper(fields)
  all_fields <- AnnotationDbi::columns(orgdb)
  chosen_fields <- c()

  if (! name_column %in% all_fields) {
    a_name <- grepl(pattern = "NAME", x = all_fields)
    new_name_column <- all_fields[a_name][1]
    message("Unable to find ", name_column, ", setting it to ", new_name_column, ".")
    name_column <- new_name_column
  }
  if (! type_column %in% all_fields) {
    message("Unable to find ", type_column, " in the db, removing it.")
    type_column <- NULL
  }
  if (! chromosome_column %in% all_fields) {
    message("Unable to find ", chromosome_column, " in the db, removing it.")
    chromosome_column <- NULL
  }
  if (! strand_column %in% all_fields) {
    message("Unable to find ", strand_column, " in the db, removing it.")
    strand_column <- NULL
  }
  if (! start_column %in% all_fields) {
    message("Unable to find ", start_column, " in the db, removing it.")
    start_column <- NULL
  }
  if (! end_column %in% all_fields) {
    message("Unable to find ", end_column, " in the db, removing it.")
    end_column <- NULL
  }

  if (length(fields) == 0) {
    chosen_fields <- c(name_column, type_column, chromosome_column, strand_column,
                       start_column, end_column)
  } else if (length(fields) == 1 && grepl(x = fields, pattern = "\\^")) {
    field_idx <- grepl(x = all_fields, pattern = toupper(fields))
    fields <- all_fields[field_idx]
    chosen_fields <- c(name_column, type_column, chromosome_column, strand_column,
                       start_column, end_column, fields)
  } else {
    chosen_fields <- c(name_column, type_column, chromosome_column, strand_column,
                       start_column, end_column, fields)
  }

  if ("ALL" %in% chosen_fields) {
    message("Selecting the following fields, this might be too many: \n",
            toString(all_fields))
    chosen_fields <- all_fields
  } else {
    if (sum(chosen_fields %in% all_fields) != length(chosen_fields)) {
      missing_idx <- ! chosen_fields %in% all_fields
      missing_fields <- chosen_fields[missing_idx]
      found_fields <- chosen_fields %in% all_fields
      chosen_fields <- chosen_fields[found_fields]
      message("Some requested columns are not available: ", toString(missing_fields), ".")
      message("The following are available: ", toString(all_fields))
    }
  }

  ## Gene IDs
  if (is.null(gene_ids)) {
    gene_ids <- try(AnnotationDbi::keys(orgdb, keytype = keytype))
    if (class(gene_ids) == "try-error") {
      if (grepl(x = gene_ids[[1]], pattern = "Invalid keytype")) {
        valid_keytypes <- AnnotationDbi::keytypes(orgdb)
        stop("Try using valid keytypes: ", toString(valid_keytypes))
      } else {
        stop("There was an error getting the gene ids.")
      }
    } else {
      message("Extracted all gene ids.")
    }
  }
  ## Note querying by "GENEID" will exclude noncoding RNAs
  message("Attempting to select: ", toString(chosen_fields))
  gene_info <- try(AnnotationDbi::select(
                                      x = orgdb,
                                      keys = gene_ids,
                                      keytype = keytype,
                                      columns = chosen_fields))
  if (class(gene_info) == "try-error") {
    message("Select statement failed, this is commonly because there is no join",
            " between the transcript table and others.")
    message("Thus it says some stupid crap about 'please add gtc to the interpolator'",
            " which I think references select-method.R in GenomicFeatures.")
    message("So, try replacing columns with stuff like 'tx*' with 'cds*'?")
    stop()
  }

  ## Compute total transcript lengths (for all exons)
  ## https://www.biostars.org/p/83901/
  gene_exons <- try(GenomicFeatures::exonsBy(orgdb, by = "gene"), silent = TRUE)
  if (class(gene_exons) == "try-error") {
    gene_exons <- NULL
  }
  transcripts <- try(GenomicFeatures::transcripts(orgdb), silent = TRUE)
  if (class(transcripts) == "try-error") {
    transcripts <- NULL
  }
  fivep_utr <- try(GenomicFeatures::fiveUTRsByTranscript(orgdb, use.names = TRUE), silent = TRUE)
  if (class(fivep_utr) == "try-error") {
    fivep_utr <- NULL
  }
  threep_utr <- try(GenomicFeatures::threeUTRsByTranscript(orgdb, use.names = TRUE), silent = TRUE)
  if (class(threep_utr) == "try-error") {
    threep_utr <- NULL
  }
  colnames(gene_info) <- tolower(colnames(gene_info))
  if (isTRUE(sum_exon_widths)) {
    message("Summing exon lengths, this takes a while.")
    lengths <- lapply(gene_exons, function(x) {
      sum(BiocGenerics::width(GenomicRanges::reduce(x)))
    })
    message("Adding exon lengths to the gene_exons.")
    lengths <- as.data.frame(unlist(lengths), stringsAsFactors = FALSE)
    colnames(lengths) <- "transcript_length"
    gene_info <- merge(gene_info, lengths, by.x = keytype, by.y = "row.names")
  }
  rownames(gene_info) <- make.names(gene_info[[1]], unique = TRUE)

  retlist <- list(
      "genes" = gene_info,
      "gene_exons" = gene_exons,
      "transcripts" = transcripts,
      "fivep_utr" = fivep_utr,
      "threep_utr" = threep_utr)
  class(retlist) <- "orgdb_annotations"
  return(retlist)
}

#' Avoid annoyingly large print results of orgdb data.
#'
#' @param x Result from load_orgdb_annotations().
#' @param ... pass along args
#' @export
print.orgdb_annotations <- function(x, ...) {
  result_string <- glue("A set of orgdb annotations including: {nrow(x[['genes']])} gene annotations.")
  return(result_string)
}

#' Retrieve GO terms associated with a set of genes.
#'
#' AnnotationDbi provides a reasonably complete set of GO mappings between gene
#' ID and ontologies.  This will extract that table for a given set of gene
#' IDs.
#'
#' This is a nice way to extract GO data primarily because the Orgdb data sets
#' are extremely fast and flexible, thus by changing the keytype argument, one
#' may use a lot of different ID types and still score some useful ontology data.
#'
#' @param orgdb OrganismDb instance.
#' @param gene_ids Identifiers of the genes to retrieve annotations.
#' @param keytype The mysterious keytype returns yet again to haunt my dreams.
#' @param columns The set of columns to request.
#' @param guess_columns Instead of a set of specific columns, use grep to find anything with 'go'
#' @return Data frame of gene IDs, go terms, and names.
#' @seealso [AnnotationDbi] [GO.db]
#' @example inst/examples/annotation_orgdb.R
#' @author I think Keith provided the initial implementation of this, but atb
#'  messed with it pretty extensively.
#' @export
load_orgdb_go <- function(orgdb = NULL, gene_ids = NULL, keytype = "ensembl",
                          columns = c("go", "goall", "goid"), guess_columns = FALSE) {
                          ## columns = "go", rbind = TRUE) {
  if (is.null(orgdb)) {
    message("Assuming Homo.sapiens.")
    org_loaded <- do.call("library", list("package" = "Homo.sapiens", "character.only" = TRUE))
    orgdb <- get0("Homo.sapiens")
  } else if ("character" %in% class(orgdb)) {
    org_loaded <- do.call("library", list("package" = orgdb, "character.only" = TRUE))
    orgdb <- get0(orgdb)
  }
  tt <- sm(requireNamespace("GO.db"))
  keytype <- toupper(keytype)
  columns <- toupper(columns)
  if (isTRUE(guess_columns)) {
    available <- AnnotationDbi::keytypes(orgdb)
    column_idx <- grepl(x = available, pattern = "GO")
    columns <- available[column_idx]
  }
  if (is.null(gene_ids)) {
    gene_ids <- try(AnnotationDbi::keys(orgdb, keytype = keytype), silent = TRUE)
    if (class(gene_ids) == "try-error") {
      avail_types <- AnnotationDbi::keytypes(orgdb)
      if ("GID" %in% avail_types) {
        message("The chosen keytype was not available.  Using 'GID'.")
        keytype <- "GID"
        gene_ids <- AnnotationDbi::keys(orgdb, keytype = keytype)
      } else {
        keytype <- avail_types[[1]]
        message("Neither the chosen keytype, nor 'GID' was available.
The available keytypes are: ", toString(avail_types), "choosing ", keytype, ".")
        gene_ids <- AnnotationDbi::keys(orgdb, keytype = keytype)
      }
    }
  }
  if (class(orgdb)[[1]] == "OrganismDb") {
    message("This is an organismdbi, that should be ok.")
  } else if (class(orgdb)[[1]] == "OrgDb" | class(orgdb)[[1]] == "orgdb") {
    message("This is an orgdb, good.")
  } else {
    stop("This requires either an organismdbi or orgdb instance, not ", class(orgdb)[[1]])
  }
  available_columns <- AnnotationDbi::columns(orgdb)
  chosen_columns <- c()
  for (col in columns) {
    if (col %in% available_columns) {
      chosen_columns <- c(chosen_columns, col)
    }
  }
  if (is.null(chosen_columns)) {
    stop("Did not find any of: ", toString(columns),
         " in the set of available columns: ", toString(available_columns))
  }

  go_terms <- data.frame()
  passed_columns <- 0
  observed_colnames <- c()
  for (column in chosen_columns) {
    mesg("Getting all rows from ", column, ".")
    column_terms <- try(AnnotationDbi::select(x = orgdb,
                                              keys = gene_ids,
                                              keytype = keytype,
                                              columns = column))
    if (class(column_terms) == "try-error") {
      if (grep(pattern = "Invalid keytype", x = go_terms[[1]])) {
        message("Here are the possible keytypes:")
        message(toString(AnnotationDbi::keytypes(orgdb)))
      }
    } else {
      passed_columns <- passed_columns + 1
    }
    if (passed_columns == 1) {
      observed_colnames <- colnames(column_terms)
      go_terms <- column_terms
    } else {
      mesg("Adding ", column, " entries.")
      colnames(column_terms) <- observed_colnames
      go_terms <- rbind.data.frame(go_terms, column_terms)
      dup_idx <- duplicated(go_terms)
      mesg("Deduplicating, removing ", sum(dup_idx), " rows.")
      go_terms <- go_terms[!dup_idx, ]
    }
  }
  if (passed_columns == 0) {
    stop("None of the go columns provided information.")
  }

  if ("GO" %in% chosen_columns) {
    go_terms <- go_terms[!is.na(go_terms[["GO"]]), ]
    go_term_names <- sm(AnnotationDbi::select(x = GO.db::GO.db,
                                              keys = unique(go_terms[["GO"]]),
                                              columns = c("TERM", "GOID", "ONTOLOGY")))
    go_terms <- merge(go_terms, go_term_names, by.x = "GO", by.y = "GOID")
  }
  return(go_terms)
}

load_txdb_annotations <- function(txdb = NULL, gene_ids = NULL, types = c("tx", "exon", "cds")) {
  if (is.null(txdb)) {
    message("Assuming Homo.sapiens.")
    tx_loaded <- do.call("library", args = list("package" = "TxDb.Hsapiens.UCSC.hg19.knownGene", "character.only" = TRUE))
    txdb <- get0("TxDb.Hsapiens.UCSC.hg19.knownGene")
  } else if (class(txdb) == "character") {
    tx_loaded <- do.call("library", list("package" = txdb, "character.only" = TRUE))
    txdb <- get0(txdb)
  }
  possible_keytypes <- keytypes(txdb)
  possible_columns <- columns(txdb)
  retlist <- list()
  gene_ids <- try(AnnotationDbi::keys(txdb, keytype = "GENEID"))
  retlist[["genes"]] <- gene_ids
  for (type in types) {
    type <- toupper(type)
    retlist[[type]] <- data.frame()
    keytype <- paste0(type, "ID")
    type_ids <- AnnotationDbi::keys(txdb, keytype = keytype)
    wanted_idx <- grepl(pattern = glue("^{type}"), x = possible_columns)
    wanted_columns <- possible_columns[wanted_idx]
    type_info <- try(AnnotationDbi::select(
      x = txdb,
      keys = gene_ids,
      keytype = "GENEID",
      columns = wanted_columns))
    retlist[[type]] <- type_info
  }
  class(retlist) <- "hpgltools::load_txdb_annotations"
  return(retlist)
}

#' Given 2 species names from the eupathdb, make orthology tables betwixt them.
#'
#' The eupathdb provides such a tremendous wealth of information.  For me
#' though, it is difficult sometimes to boil it down into just the bits of
#' comparison I want for 1 species or between 2 species.  A singularly common
#' question I am asked is: "What are the most similar genes between species x
#' and y among these two arbitrary parasites?"  There are lots of ways to poke
#' at this question: run BLAST/fasta36, use biomart, query the ortholog tables
#' from the eupathdb, etc.  However, in all these cases, it is not trivial to
#' ask the next question:  What about: a:b and b:a?
#' This function attempts to address that for the case of two eupath species
#' from the same domain. (tritrypdb/fungidb/etc.)  It does however assume that
#' the sqlite package has been installed locally, if not it suggests you run the
#' make_organismdbi function in order to do that.
#'
#' One other important caveat: this function assumes queries in the format
#' 'table_column' where in this particular instance, the table is further
#' assumed to be the ortholog table.
#'
#' @param db Species name (subset) from one eupath database.
#' @param master Primary keytype to use for indexing the various tables.
#' @param query_species A list of exact species names to search for.  If uncertain
#'  about them, add print_speciesnames=TRUE and be ready for a big blob of
#'  text.  If left null, then it will pull all species.
#' @param id_column What column in the database provides the set of ortholog IDs?
#' @param org_column What column provides the species name?
#' @param group_column Ortholog group column name.
#' @param name_column Name of the gene for this group.
#' @param count_column Name of the column with the count of species represented.
#' @param print_speciesnames Dump the species names for diagnostics?
#' @param webservice Which eupathdb project to query?
#' @return A big table of orthoMCL families, the columns are:
#'  \enumerate{
#'   \item  GID: The gene ID
#'   \item  ORTHOLOG_ID: The gene ID of the associated ortholog.
#'   \item  ORTHOLOG_SPECIES: The species of the associated ortholog.
#'   \item  ORTHOLOG_URL: The OrthoMCL group ID's URL.
#'   \item  ORTHOLOG_COUNT: The number of all genes from all species represented in
#'   this group.
#'   \item  ORTHOLOG_GROUP: The family ID
#'   \item  QUERIES_IN_GROUP: How many of the query species are represented in this
#'   group?
#'   \item  GROUP_REPRESENTATION: ORTHOLOG_COUNT / the number of possible species.
#'  }
#' @author atb
#' @export
extract_eupath_orthologs <- function(db, master = "GID", query_species = NULL,
                                     id_column = "ORTHOLOGS_GID",
                                     org_column = "ORTHOLOGS_ORGANISM",
                                     group_column = "ANNOT_GENE_ORTHOMCL_NAME",
                                     name_column = "ORTHOLOGS_PRODUCT",
                                     count_column = "ORTHOLOGS_COUNT",
                                     print_speciesnames = FALSE,
                                     webservice = "eupathdb") {

  pkg <- NULL
  if (class(db)[1] == "OrgDb") {
    pkg <- db
  } else {
    stop("I only understand orgdbs or the name of a species.")
  }

  columns <- c(id_column, group_column, org_column, name_column, count_column)
  columns <- toupper(columns)
  gene_set <- AnnotationDbi::keys(pkg, keytype = master)
  column_set <- AnnotationDbi::columns(pkg)
  column_intersect <- columns %in% column_set
  if (sum(column_intersect) == length(columns)) {
    message("Found all the required columns!")
  } else {
    missing_idx <- ! columns %in% column_set
    missing <- columns[missing_idx]
    message("Some columns were missing: ", toString(missing))
    message("Removing them, which may end badly.")
    columns <- columns[column_intersect]
  }
  all_orthos <- AnnotationDbi::select(x = pkg, keytype = master,
                                      keys = gene_set, columns = columns)
  all_orthos[[org_column]] <- as.factor(all_orthos[[org_column]])
  num_possible <- 1
  species_names <- levels(all_orthos[[org_column]])
  if (is.null(query_species)) {
    query_species <- species_names
  } else if (! query_species %in% species_names) {
    warning("Did not find the desired species in the set of all species.")
    query_species <- species_names
  }
  num_possible <- length(species_names)
  message("There are ", num_possible, " possible species in this group.")

  if (isTRUE(print_speciesnames)) {
    print(toString(species_names))
    return(invisible())
  }

  ## Now pull out the species of interest
  found_species <- 0
  for (sp in query_species) {
    if (sp %in% all_orthos[[org_column]]) {
      message("Found species: ", sp)
    } else {
      message("Did not find species: ", sp)
    }
  }
  kept_orthos_idx <- all_orthos[[org_column]] %in% query_species
  kept_orthos <- all_orthos[kept_orthos_idx, ]
  ## The following is not possible if we used the orthologslite table.
  ## In fact, the orthologslite table is (I am realizing) quite a disappointment.
  ## I might remove that query and just force the much slower orthologs table as it
  ## provides much more useful information.
  if (is.null(all_orthos[["ORTHOLOGS_COUNT"]])) {
    kept_orthos_dt <- data.table::as.data.table(kept_orthos)
  } else {
    colnames(kept_orthos) <- c(master, "ORTHOLOGS_ID", "ORTHOLOGS_GROUP", "ORTHOLOGS_SPECIES",
                               "ORTHOLOGS_NAME", "ORTHOLOGS_COUNT")
    kept_orthos[["ORTHOLOGS_COUNT"]] <- as.integer(kept_orthos[["ORTHOLOGS_COUNT"]])
    GID <- NULL
    kept_orthos_dt <- data.table::as.data.table(kept_orthos) %>%
      dplyr::group_by(GID) %>%
      dplyr::add_count(GID)
    colnames(kept_orthos_dt) <- c(master, "ORTHOLOGS_ID", "ORTHOLOGS_GROUP",
                                  "ORTHOLOGS_SPECIES", "ORTHOLOGS_NAME", "ORTHOLOGS_COUNT",
                                  "QUERIES_IN_GROUP")
    kept_orthos_dt[["ORTHOLOGS_REPRESENTATION"]] <- kept_orthos_dt[["ORTHOLOGS_COUNT"]] / num_possible
    num_queries <- length(query_species)
  }
  return(kept_orthos_dt)
}

#' Map AnnotationDbi keys from one column to another.
#'
#' Given a couple of keytypes, this provides a quick mapping across them.  I
#' might have an alternate version of this hiding in the gsva code, which
#' requires ENTREZIDs.  In the mean time, this creates a dataframe of the mapped
#' columns for a given set of gene ids using the in a sqlite instance.
#'
#' @param orgdb OrganismDb instance.
#' @param gene_ids Gene identifiers for retrieving annotations.
#' @param mapto Key to map the IDs against.
#' @param keytype Choose a keytype, this will yell if it doesn't like your choice.
#' @return a table of gene information
#' @seealso [AnnotationDbi]
#' @example inst/examples/annotation_orgdb.R
#' @author Keith Hughitt with changes by atb.
#' @export
map_orgdb_ids <- function(orgdb, gene_ids = NULL, mapto = "ensembl",
                          keytype = "geneid") {
  if (is.null(orgdb)) {
    message("Assuming Homo.sapiens.")
    org_loaded <- do.call("library", list("package" = "Homo.sapiens", "character.only" = TRUE))
    orgdb <- get0("Homo.sapiens")
  } else if ("character" %in% class(orgdb)) {
    org_loaded <- do.call("library", list("package" = orgdb, "character.only" = TRUE))
    orgdb <- get0(orgdb)
  }
  mapto <- toupper(mapto)
  keytype <- toupper(keytype)
  avail_keytypes <- AnnotationDbi::keytypes(orgdb)
  found_keys <- sum(mapto %in% avail_keytypes)
  if (found_keys < length(mapto)) {
    warning("The chosen keytype ", mapto, " is not in this orgdb.")
    warning("Try some of the following instead: ", toString(avail_keytypes), ".")
    warning("Going to pull all the availble keytypes, which is probably not what you want.")
    mapto <- avail_keytypes
  }

  test_masterkey <- sum(keytype %in% avail_keytypes)
  if (test_masterkey != 1) {
    warning("The chosen master key ", keytype, " is not in this orgdb.")
    warning("Try some of the following instead: ", toString(avail_keytypes), ".")
    warning("I am going to choose one arbitrarily, which is probably not what you want.")
    if ("ENTREZID" %in% avail_keytypes) {
      keytype <- "ENTREZID"
      message("Using entrezid as the master key.")
    } else if ("ENSEMBLID" %in% avail_keytypes) {
      keytype <- "ENSEMBLID"
      message("Using ensemblid as the master key.")
    } else
      stop("Could not think of a usable master key.")
  }

  ## If no gene ids were chosen, grab them all.
  if (is.null(gene_ids)) {
    gene_ids <- AnnotationDbi::keys(orgdb, keytype = keytype)
  }
  ## Gene info
  ## Note querying by "GENEID" will exclude noncoding RNAs
  gene_info <- AnnotationDbi::select(x = orgdb, keytype = keytype,
                                     keys = gene_ids, columns = mapto)
  colnames(gene_info) <- tolower(colnames(gene_info))
  return(gene_info)
}

#' Iterate over keytypes looking for matches against a set of IDs.
#'
#' Sometimes, one does not know what the correct keytype is for a given set of
#' IDs.  This will hopefully find them.
#'
#' @param ids Set of gene IDs to seek.
#' @param orgdb Orgdb instance to iterate through.
#' @param verbose talky talk
#' @return Likely keytype which provides the desired IDs.
#' @seealso [org.Dm.eg.db]
#' @example inst/examples/annotation_orgdb.R
#' @export
guess_orgdb_keytype <- function(ids, orgdb = NULL, verbose = FALSE) {
  if (is.null(orgdb)) {
    message("Assuming Homo.sapiens.")
    lib <- do.call(what = "library", args = list("package" = "Homo.sapiens", "character.only" = TRUE))
    orgdb <- get0("Homo.sapiens")
  } else if ("character" %in% class(orgdb)) {
    lib <- do.call(what = "library", args = list("package" = orgdb, "character.only" = TRUE))
    orgdb <- get0(orgdb)
  }
  found_ids <- 0
  current_type <- NULL
  possible_keytypes <- AnnotationDbi::keytypes(orgdb)
  for (k in seq_along(possible_keytypes)) {
    type <- possible_keytypes[k]
    possible_keys <- try(AnnotationDbi::keys(x = orgdb, keytype = type), silent = TRUE)
    if (class(possible_keys)[1] == "try-error") {
      possible_keys <- ""
    }
    this_type_found <- sum(ids %in% possible_keys)
    if (isTRUE(verbose)) {
      message("Keytype: ", type, " has ", this_type_found, " keys.")
    }
    if (this_type_found == length(ids)) {
      return(type)
    } else if (this_type_found > found_ids) {
      current_type <- type
      found_ids <- this_type_found
    }
  }
  if (found_ids == 0) {
    message("Did not find your IDs using any keytype in the orgdb.")
  } else {
    message("The best choice was ", current_type, " which has ", found_ids,
            " out of ", length(ids), " ids.")
  }
  return(current_type)
}

#' Guess the orgdb from a genusspecies.
#'
#' Given a name like 'mmusculus', guess the orgdb package name.
#' @param species Input species
#' @param genus and genus.
#' @examples
#'  guess <- map_species_orgdb("hsapiens")
map_species_orgdb <- function(species, genus = NULL) {
  if (is.null(genus)) {
    ## Then assume things like 'hsapiens'.
    chars <- strsplit(species, "")[[1]]
    first <- toupper(chars[1])
    second <- tolower(chars[2])
  } else {
    first_chars <- strsplit(genus, "")[[1]]
    second_chars <- strsplit(species, "")[[1]]
    first <- toupper(first_chars[1])
    second <- tolower(second_chars[1])
  }
  guess <- glue("org.{first}{second}.eg.db")
  return(guess)
}

#' Get an orgdb from an AnnotationHub taxonID.
#'
#' Ideally, annotationhub will one day provide a one-stop shopping source for a
#' tremendous wealth of curated annotation databases, sort of like a
#' non-obnoxious biomart.  But for the moment, this function is more
#' fragile than I would like.
#'
#' @param ahid TaxonID from AnnotationHub
#' @param title Title for the annotation hub instance
#' @param species Species to download
#' @param type Datatype to download
#' @return An Orgdb instance
#' @seealso [AnnotationHub] [S4Vectors]
#' @examples
#' \dontrun{
#'  org <- mytaxIdToOrgDb(species = "Leishmania", type = "TxDb")
#' }
#' @export
orgdb_from_ah <- function(ahid = NULL, title = NULL, species = NULL, type = "OrgDb") {
  ## Other available types:
  tt <- sm(loadNamespace("AnnotationHub"))
  ah <- sm(AnnotationHub::AnnotationHub())
  message("Available types: \n", toString(levels(as.factor(ah$rdataclass))))

  if (!is.null(type)) {
    ah <- AnnotationHub::query(x = ah, pattern = type)
  }
  if (is.null(title) & is.null(species) & is.null(ahid)) {
    message("Going to attempt to find a human database.  I hope this is what you want!")
    hits <- grepl(pattern = "Hs\\.eg\\.db", x = ah$title)
    ahid <- names(ah)[hits]
  } else if (!is.null(species)) {
    ## Then we got a species
    possible <- ah$species
    titles <- ah$title
    hits_idx <- grepl(pattern = species, x = possible)
    first_true <- which.max(hits_idx)
    first_true_name <- titles[first_true]
    hits <- names(ah)[hits_idx]
    message("The possible hits are: \n",
            toString(hits), "\nchoosing: ", hits[1],
            "\nwhich is ", first_true_name)
    ahid <- hits[1]
  } else if (!is.null(title)) {
    ## We got a title
    possible <- ah$title
    hits_idx <- grepl(pattern = title, x = possible)
    first_true <- which.max(hits_idx)
    first_true_name <- possible[first_true]
    hits <- names(ah)[hits_idx]
    message("The possible hits are: \n",
            toString(hits), "\nchoosing: ", hits[1],
            "\nwhich is ", first_true_name)
    ahid <- hits[1]
  }

  ah_names <- names(ah)
  ah_titles <- ah$title
  hit_idx <- ah_names == ahid
  hit_num <- which.max(hit_idx)
  hit_title <- ah_titles[hit_num]
  message("Chose ", ahid, " which is ", hit_title, ".")
  res <- ah[[ahid]]
  return(res)
}

## EOF
