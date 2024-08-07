#' Use AnnotationDbi to translate geneIDs from type x to type y.
#'
#' This is intended to convert all the IDs in a geneSet from one ID type to
#' another and giving back the geneSet with the new IDs.
#' FIXME: This should use convert_ids() to simplify itself
#'
#' One caveat: this will collapse redundant IDs via unique().
#'
#' @param gsc geneSetCollection with IDs of a type one wishes to change.
#' @param orgdb Annotation object containing the various IDs.
#' @param from_type Name of the ID which your gsc is using.  This can probably
#'  be automagically detected...
#' @param to_type Name of the ID you wish to use.
#' @return Fresh gene set collection replete with new names.
#' @seealso [AnnotationDbi] [guess_orgdb_keytypes()] [convert_ids()] [GSEABase]
#' @export
convert_gsc_ids <- function(gsc, orgdb = "org.Hs.eg.db", from_type = NULL, to_type = "ENTREZID") {
  message("Converting the rownames() of the expressionset to ", to_type, ".")
  ##tt <- sm(library(orgdb, character.only = TRUE))
  lib_result <- sm(requireNamespace(orgdb))
  att_result <- sm(try(attachNamespace(orgdb), silent = TRUE))
  orgdb <- get0(orgdb)
  gsc_lst <- as.list(gsc)
  new_gsc <- list()
  if (isTRUE(verbose)) {
    bar <- utils::txtProgressBar(style = 3)
  }
  for (g in seq_len(gsc)) {
    if (isTRUE(verbose)) {
      pct_done <- g / length(gsc_lst)
      setTxtProgressBar(bar, pct_done)
    }
    gs <- gsc[[g]]
    old_ids <- GSEABase::geneIds(gs)
    if (is.null(from_type)) {
      from_type <- guess_orgdb_keytype(old_ids, orgdb)
    }
    new_ids <- convert_ids(old_ids, from = from_type, to = to_type, orgdb = orgdb)
    GSEABase::geneIds(gs) <- unique(new_ids)
    ## gsc_lst[[g]] <- gs
    new_gsc[[g]] <- gs
  }
  if (isTRUE(verbose)) {
    close(bar)
  }
  gsc <- GSEABase::GeneSetCollection(new_gsc)
  return(gsc)
}

#' Get an Identifier function from a shorthand name.
#'
#' I am hoping to write one for EuPathDB and some other source, thus the switch.
#'
#' @param type String name for the identifier in question.
get_identifier <- function(type) {
  id_function <- GSEABase::NullIdentifier
  switchret <- switch(
      type,
      "entrez" = {
        id_function <- GSEABase::EntrezIdentifier
      },
      "symbol" = {
        id_function <- GSEABase::SymbolIdentifier
      },
      "ensembl" = {
        id_function <- GSEABase::ENSEMBLIdentifier
      },
      "uniprot" = {
        id_function <- GSEABase::UniprotIdentifier
      },
      "refseq" = {
        id_function <- GSEABase::RefseqIdentifier
      },
      "name" = {
        id_function <- GSEABase::GenenameIdentifier
      },
      {
        id_function <- GSEABase::NullIdentifier
      })
  return(id_function)
}

#' Create a metadata dataframe of msigdb data, this hopefully will be usable to
#' fill the fData slot of a gsva returned expressionset.
#'
#' @param msig_db Filename containing the MSigDB metadata.
#' @param wanted_meta Choose metadata columns of interest.
#' @return list containing 2 data frames: all metadata from broad, and the set
#'  matching the sig_data GeneSets.
#' @seealso [xml2] [rvest]
#' @export
get_msigdb_metadata <- function(msig_db = "msigdb_v6.2.xml",
                                wanted_meta = c("ORGANISM", "DESCRIPTION_BRIEF", "AUTHORS", "PMID")) {
  list_result <- list()
  all_data <- parse_msigdb(msig_db)

  ## all_data is my dataframe of xml annotations, lets extract the wanted columns before
  ## adding it to the gsva result expressionset.
  sub_data <- all_data
  if (is.null(wanted_meta)) {
    message("Not subsetting the msigdb metadata, the wanted_meta argument was 'all'.")
  } else if (wanted_meta[1] == "all") {
    message("Not subsetting the msigdb metadata, the wanted_meta argument was NULL.")
  } else if (class(wanted_meta)[1] == "character") {
    sub_data <- sub_data[, wanted_meta]
  }

  return(sub_data)
}

#' Parse either xml or sqlite data from MSigDB.  I think I will likely remove the
#' xml version as I think the msigdb xml files are poorly formatted.
#'
#' @param filename Input file
parse_msigdb <- function(filename) {
  ex <- tools::file_ext(filename)
  if (ex == "xml") {
    parse_msigdb_xml(filename)
  } else if (ex == "db") {
    parse_msigdb_sqlite(filename)
  } else {
    stop("I do not know this file type.")
  }
}

#' Extract metadata from the msigdb sqlite data.
#'
#' @param filename db file.
#'
#' @importFrom RSQLite dbConnect dbFetch dbClearResult dbSendQuery
parse_msigdb_sqlite <- function(filename) {
  db <- RSQLite::dbConnect(RSQLite::SQLite(), dbname = filename)
  author_query <- dbSendQuery(conn = db, "SELECT * FROM author")
  author_df <- dbFetch(author_query)
  cleared <- dbClearResult(author_query)
  collection_query <- dbSendQuery(conn = db, "SELECT * FROM collection")
  collection_df <- dbFetch(collection_query)
  cleared <- dbClearResult(collection_query)
  pub_query <- dbSendQuery(conn = db, "SELECT * FROM publication")
  pub_df <- dbFetch(pub_query)
  cleared <- dbClearResult(pub_query)
  gs_query <- dbSendQuery(conn = db, "SELECT * FROM gene_set")
  gs_df <- dbFetch(gs_query)
  cleared <- dbClearResult(gs_query)
  detail_query <- dbSendQuery(conn = db, "SELECT * FROM gene_set_details")
  detail_df <- dbFetch(detail_query)
  cleared <- dbClearResult(detail_query)
  closed <- RSQLite::dbDisconnect(db)
  combined <- merge(gs_df, detail_df, by.x = "id", by.y = "gene_set_id")
  combined <- merge(combined, pub_df, by = "id", all.x = TRUE)
  rownames(combined) <- combined[["standard_name"]]
  combined[["standard_name"]] <- NULL
  return(combined)
}

## Unfortunately, as of ~ 2023 the xml provided by the MsigDB appears
## to by syntactically invalid XML.  As a result, the following dies
## immediately with "Unescaped '<' is not allowed".

#' Extract fun experimental metadata from a MSigDB xml file.
#'
#' @param filename input file.
parse_msigdb_xml <- function(filename) {
  msig_result <- xml2::read_xml(x = filename)

  ##db_data <- rvest::xml_nodes(x = msig_result, xpath = "//MSIGDB")
  db_data <- rvest::html_elements(x = msig_result, xpath = "//MSIGDB")
  db_name <- rvest::html_attr(x = db_data, name = "NAME")
  db_ver <- rvest::html_attr(x = db_data, name = "VERSION")
  db_date <- rvest::html_attr(x = db_data, name = "BUILD_DATE")

  ##genesets <- rvest::xml_nodes(x = msig_result, "GENESET")
  genesets <- rvest::html_elements(x = msig_result, "GENESET")
  row_names <- rvest::html_attr(x = genesets, name = "STANDARD_NAME")
  column_names <- names(rvest::html_attrs(x = genesets[[1]]))
  all_data <- data.frame(row.names = row_names)
  for (i in seq(from = 2, to = length(column_names))) {
    c_name <- column_names[i]
    c_data <- rvest::html_attr(x = genesets, name = c_name)
    all_data[[c_name]] <- c_data
  }

}

#' Load signatures from either a gmt file, xml file, or directly from the
#' GSVAdata data set in R.
#'
#' There are a bunch of places from which to acquire signature data.  This
#' function attempts to provide a single place to load them.  The easiest way to
#' get up to date signatures is to download them from msigdb and set the
#' signatures parameter to the downloaded filename.
#'
#' @param signatures Either the filename downloaded or the variable's name as
#'  found in the environment created by data_pkg.
#' @param data_pkg Used when signatures is not a filename to load a data
#'  package, presumably GSVAdata.
#' @param signature_category Probably not needed unless you download a signature
#'  file containing lots of different categories.
#' @param id_type Specify the ID type in the data.
#' @return signature dataset which may be used by gsva()
#' @seealso [GSEABase]
#' @export
load_gmt_signatures <- function(signatures = "c2BroadSets", data_pkg = "GSVAdata",
                                signature_category = "c2", id_type = "entrez") {
  sig_data <- NULL

  id_function <- get_identifier(id_type)

  if (class(signatures)[1] == "character" && grepl(pattern = "\\.gmt$", x = signatures)) {
    sig_data <- GSEABase::getGmt(signatures,
                                 collectionType = GSEABase::BroadCollection(category = signature_category),
                                 geneIdType = id_function())
  } else if (class(signatures)[1] == "character" && grepl(pattern = "\\.xml$", x = signatures)) {
    gsc <- GSEABase::getBroadSets(signatures)
    types <- sapply(gsc, function(elt) GSEABase::bcCategory(GSEABase::collectionType(elt)))
    sig_data <- gsc[types == signature_category]
  } else if (class(signatures)[1] == "character") {
    lib_result <- sm(requireNamespace(data_pkg))
    att_result <- sm(try(attachNamespace(data_pkg), silent = TRUE))
    lst <- list("list" = signatures, "package" = data_pkg)
    test <- do.call("data", as.list(signatures, lst))
    sig_data <- get0(signatures)
  } else if (class(signatures)[1] != "GeneSetCollection") {
    stop("The data must be a GeneSetCollection.")
  } else {
    sig_data <- signatures
  }
  mesg("Loaded ", length(sig_data), " gene sets.")
  return(sig_data)
}

#' Create a gene set collection from a set of arbitrary IDs.
#'
#' This function attempts to simplify the creation of a gsva compatible
#' GeneSet.  Some important caveats when working with gsva, notably the gene IDs
#' we use are not usually compatible with the gene IDs used by gsva, thus the
#' primary logic in this function is intended to bridge these IDs.
#'
#' @param first_ids The required IDs for a single set.
#' @param second_ids Potentially null optionally used for a second, presumably
#'  contrasting set.
#' @param annotation_name Orgdb annotation, used to translate IDs to the required type.
#' @param researcher_name Prefix of the name for the generated set(s).
#' @param study_name Second element in the name of the generated set(s).
#' @param category_name Third element in the name of the generated set(s).
#' @param phenotype_name Optional phenotype data for the generated set(s).
#' @param identifier_type ID type to use in the gene set.
#' @param organism Set the organism for the gsc object.
#' @param pair_names The suffix of the generated set(s).
#' @param current_id What type of ID is the data currently using?
#' @param required_id What type of ID should the use?
#' @param min_gmt_genes Minimum number of genes in the set for consideration.
#' @return Small list comprised of the created gene set collection(s).
#' @seealso [GSEABase]
#' @importFrom GSEABase GeneSet GeneColorSet
#' @export
make_gsc_from_ids <- function(first_ids, second_ids = NULL, annotation_name = "org.Hs.eg.db",
                              researcher_name = "elsayed", study_name = "macrophage",
                              category_name = "infection", phenotype_name = NULL,
                              identifier_type = "entrez", organism = NULL,
                              pair_names = "up", current_id = "ENSEMBL",
                              required_id = "ENTREZID", min_gmt_genes = 10) {
  first <- NULL
  second <- NULL
  do_first <- TRUE
  do_second <- TRUE
  if (is.null(current_id) | is.null(required_id)) {
    first <- first_ids
    current_id <- ""
    second <- second_ids
    required_id <- ""
  }
  if (current_id == required_id) {
    first <- first_ids
    second <- second_ids
  } else {
    message("Converting the rownames() of the expressionset to ENTREZID.")
    ## tt <- sm(try(do.call("library", as.list(orgdb)), silent = TRUE))
    lib_result <- sm(requireNamespace(annotation_name))

    att_result <- sm(try(attachNamespace(annotation_name), silent = TRUE))
    first_ids <- sm(AnnotationDbi::select(x = get0(annotation_name),
                                          keys = first_ids,
                                          keytype = current_id,
                                          columns = c(required_id)))

    first_idx <- complete.cases(first_ids)
    if (sum(first_idx) == 0) {
      do_first <- FALSE
    } else {
      if(!all(first_idx)) {
        message(sum(first_idx) == FALSE,
                " ENSEMBL ID's didn't have a matching ENTEREZ ID in this database from first list of IDs given. Dropping them now.")
      }
      first_ids <- first_ids[first_idx, ]
      first <- first_ids[[required_id]]
    }

    if (!is.null(second_ids)) {
      second_ids <- sm(AnnotationDbi::select(x = get0(annotation_name),
                                             keys = second_ids,
                                             keytype = current_id,
                                             columns = c(required_id)))
      second_idx <- complete.cases(second_ids)
      if (sum(second_idx) == 0) {
        do_second <- FALSE
      } else {
        if(!all(second_idx)) {
          message(sum(second_idx == FALSE), " ENSEMBL ID's didn't have a matching ENTEREZ ID in this database from second list of IDs given. Dropping them now.")
        }
        second_ids <- second_ids[second_idx, ]
        second <- second_ids[[required_id]]
      }
    }
  }

  all_colored <- NULL
  fst_gsc <- NULL
  sec_gsc <- NULL
  fst <- NULL
  sec <- NULL
  identifier <- NULL
  included_first <- unique(first)
  if (!is.null(second)) {
    included_second <- unique(second)
    if (length(included_second) >= min_gmt_genes) {
      do_second <- TRUE
      sec <- data.frame(row.names = included_second)
    }
  }

  if (length(included_first) >= min_gmt_genes) {
    fst <- data.frame(row.names = unique(first))
    if (is.null(phenotype_name)) {
      phenotype_name <- "unknown"
    }
    fst[["direction"]] <- pair_names[1]
    fst[["phenotype"]] <- phenotype_name
    study_name <- gsub(x = study_name, pattern = "_", replacement = "")
    category_name <- gsub(x = category_name, pattern = "_", replacement = "")
    set_prefix <- glue("{researcher_name}_{study_name}_{category_name}")
    fst_name <- toupper(glue("{set_prefix}_{pair_names[1]}"))
    identifier <- get_identifier(identifier_type)
    first_args <- list("type" = identifier(),
                       "setName" = fst_name,
                       "geneIds" = as.character(rownames(fst)))
    if (!is.null(organism)) {
      first_args[["organism"]] <- organism
    }
    required <- requireNamespace("GSEABase")
    loaded <- loadNamespace("GSEABase")
    loaded <- try(attachNamespace("GSEABase"), silent = TRUE)
    fst_gsc <- do.call("GeneSet", first_args)
  } else {
    message("There are: ", length(included_first), " genes in the first set.")
    message("According to the min_gmt_genes parameter, this is insufficient.")
    fst_gsc <- NULL
  }

  if (isTRUE(do_second)) {
    if (is.null(phenotype_name)) {
      phenotype_name <- "unknown"
    }
    if (is.na(pair_names[2]) & pair_names[1] == "up") {
      message("Setting the second pair_names to 'down'.")
      ## Then we can assume it is down.
      pair_names <- c("up", "down")
    }
    sec_name <- toupper(glue("{set_prefix}_{pair_names[2]}"))
    sec[["direction"]] <- pair_names[2]
    sec[["phenotype"]] <- phenotype_name
    both <- rbind(fst, sec)
    color_name <- toupper(glue("{set_prefix}_{phenotype_name}"))
    second_args <- list("type" = identifier(),
                     "setName" = sec_name,
                     "geneIds" = as.character(rownames(sec)))
    colored_args <- list(
        "type" = identifier(),
        "setName" = color_name,
        "geneIds" = rownames(both),
        "phenotype" = phenotype_name,
        "geneColor" = as.factor(both[["direction"]]),
        "phenotypeColor" = as.factor(both[["phenotype"]]))
    if (!is.null(organism)) {
      second_args[["organism"]] <- organism
      colored_args[["organism"]] <- organism
    }
    sec_gsc <- base::do.call("GeneSet", second_args)
    all_colored = base::do.call("GeneColorSet", colored_args)
  } else {
    ## End testing if we should do the second gsc
    message("There are: ", length(included_second), " genes in the second set.")
    message("According to the min_gmt_genes parameter, this is insufficient.")
    sec_gsc <- NULL
  }
  retlst <- list()
  retlst[[fst_name]] <- fst_gsc
  if (isTRUE(do_second)) {
    ## Then there is a colored/down
    retlst[[sec_name]] <- sec_gsc
    retlst[["colored"]] <- all_colored
  }
  return(retlst)
}

#' Given a pairwise result, make a gene set collection.
#'
#' If I want to play with gsva and friends, then I need GeneSetCollections!  To
#' that end, this function uses extract_significant_genes() in order to gather
#' sets of genes deemed 'significant'.  It then passes these sets to
#' make_gsc_from_ids().
#'
#' @param pairwise A pairwise result, or combined de result, or extracted genes.
#' @param according_to When getting significant genes, use this method.
#' @param annotation_name Name of the annotation database for the data.
#' @param pair_names Describe the contrasts of the GSC: up vs. down, high vs. low, etc.
#' @param category_name What category does the GSC describe?
#' @param phenotype_name When making color sets, use this phenotype name.
#' @param set_name A name for the created gene set.
#' @param color Make a colorSet?
#' @param current_id Usually we use ensembl IDs, but that does not _need_ to be the case.
#' @param required_id gsva uses entrezids by default.
#' @param ... Extra arguments for extract_significant_genes().
#' @return List containing 3 GSCs, one containing both the ups/downs called
#'  'colored', one of the ups, and one of the downs.
#' @seealso [combine_de_tables()] [extract_significant_genes()] [make_gsc_from_ids()]
#'  [GSEABase]
#' @export
make_gsc_from_pairwise <- function(pairwise, according_to = "deseq", annotation_name = "org.Hs.eg.db",
                                   pair_names = c("ups", "downs"), category_name = "infection",
                                   phenotype_name = "parasite", set_name = "elsayed_macrophage",
                                   color = TRUE, current_id = "ENSEMBL", required_id = "ENTREZID",
                                   ...) {

  if (class(pairwise)[1] == "data.frame") {
    ups <- pairwise
  } else if (class(pairwise)[1] == "all_pairwise") {
    message("Invoking combine_de_tables().")
    combined <- sm(combine_de_tables(pairwise, ...))
    message("Invoking extract_significant_genes().")
    updown <- sm(extract_significant_genes(combined,
                                           according_to = according_to,
                                           ...)[[according_to]])
    ups <- updown[["ups"]]
    downs <- updown[["downs"]]
  } else if (class(pairwise)[1] == "combined_de") {
    message("Invoking extract_significant_genes().")
    updown <- sm(extract_significant_genes(pairwise,
                                           according_to = according_to,
                                           ...)[[according_to]])
    ups <- updown[["ups"]]
    downs <- updown[["downs"]]
  } else if (class(pairwise)[1] == "sig_genes") {
    updown <- pairwise[[according_to]]
    ups <- updown[["ups"]]
    downs <- updown[["downs"]]
  } else if (class(pairwise)[1] == "character") {
    message("Invoking make_gsc_from_ids().")
    ret <- make_gsc_from_ids(pairwise, annotation_name = annotation_name,
                             pair_names = pair_names, category_name = category_name,
                             phenotype_name = phenotype_name, current_id = current_id,
                             required_id = required_id, ...)
    return(ret)
  } else {
    stop("I do not understand this type of input.")
  }

  ## The rownames() of the expressionset must be in ENTREZIDs for gsva to work.
  ## tt <- sm(library(annotation_name, character.only = TRUE))
  lib_result <- sm(requireNamespace(annotation_name))
  att_result <- sm(try(attachNamespace(annotation_name), silent = TRUE))

  ## Check that the current_id and required_id are in the orgdb.
  pkg <- get0(annotation_name)
  available <- AnnotationDbi::columns(pkg)
  if (! current_id %in% available) {
    warning("The column: ", current_id, " is not in the set of available orgdb columns.")
    message("Here are the available columns:")
    print(available)
    message("Arbitrarily choosing 'GID'.")
    current_id <- "GID"
  }

  if (! required_id %in% available) {
    warning("The column: ", required_id, " is not in the set of available orgdb columns.")
    message("Here are the available columns:")
    print(available)
    message("Arbitrarily choosing 'GID'.")
    required_id <- "GID"
  }

  up_lst <- list()
  down_lst <- list()
  colored_lst <- list()
  for (c in seq_len(ups)) {
    name <- names(ups)[c]
    up <- ups[[name]]
    down <- downs[[name]]
    if (class(up)[1] == "character") {
      up_ids <- up
      down_ids <- down
    } else if (class(up)[1] == "data.frame") {
      up_ids <- rownames(up)
      down_ids <- rownames(down)
    }
    if (current_id == required_id) {
      up[[required_id]] <- rownames(up)
      down[[required_id]] <- rownames(down)
    } else {
      mesg("Converting the rownames() of the expressionset to ", required_id, ".")
      up_ids <- sm(AnnotationDbi::select(x = get0(annotation_name),
                                         keys = up_ids,
                                         keytype = current_id,
                                         columns = c(required_id)))
      up_idx <- complete.cases(up_ids)
      up_ids <- up_ids[up_idx, ]
      up <- merge(up, up_ids, by.x = "row.names", by.y = current_id)
      if (!is.null(down_ids)) {
        down_ids <- sm(AnnotationDbi::select(x = get0(annotation_name),
                                             keys = down_ids,
                                             keytype = current_id,
                                             columns = c(required_id)))
        down_idx <- complete.cases(down_ids)
        down_ids <- down_ids[down_idx, ]
        down <- merge(down, down_ids, by.x = "row.names", by.y = current_id)
      } else {
        down <- NULL
      }
    }
    up[["direction"]] <- pair_names[1]
    up[["phenotype"]] <- phenotype_name
    down[["direction"]] <- pair_names[2]
    down[["phenotype"]] <- phenotype_name
    both <- rbind(up, down)
    shared_ids <- up[[required_id]] %in% down[[required_id]]
    if (sum(shared_ids) > 0) {
      warning("There are ", sum(shared_ids), " shared IDs in the up and down lists.")
      shared_id <- up[shared_ids, ][[required_id]]
      shared_rows <- both[[required_id]] == shared_id
      both <- both[!shared_rows, ]
    }

    dup_elements <- duplicated(up[[required_id]])
    if (sum(dup_elements) > 0) {
      warning("There are ", sum(dup_elements),
              " non-unique elements in the IDs of the up data.")
      up <- up[!dup_elements, ]
    }
    dup_elements <- duplicated(down[[required_id]])
    if (sum(dup_elements) > 0) {
      warning("There are ", sum(dup_elements),
              " non-unique elements in the IDs of the down data.")
      down <- down[!dup_elements, ]
    }
    dup_elements <- duplicated(both[[required_id]])
    if (sum(dup_elements) > 0) {
      warning("There are ", sum(dup_elements),
              " non-unique elements in the IDs of the shared.")
      both <- both[!dup_elements, ]
    }

    ## Choose the Identifier for the colorsets.  For the moment just make it either entrez or null.
    ## Perhaps this logic should be added to get_identifier()?
    identifier <- GSEABase::NullIdentifier()
    if (grepl(x=tolower(required_id), pattern="entrez")) {
      identifier <- GSEABase::EntrezIdentifier()
    } else if (grepl(x=tolower(required_id), pattern="ensemb")) {
      identifier <- GSEABase::ENSEMBLIdentifier()
    }

    set_prefix <- glue("{set_name}_{category_name}")
    color_set_name <- toupper(glue("{set_prefix}_{phenotype_name}"))
    up_name <- toupper(glue("{set_prefix}_{pair_names[1]}"))

    colored_ids <- both[[required_id]]
    if (is.null(colored_ids)) {
      colored_gsc <- NULL
    } else {
      colored_gsc <- GSEABase::GeneColorSet(
                                   identifier,
                                   setName = color_set_name,
                                   geneIds = as.character(both[[required_id]]),
                                   phenotype = phenotype_name,
                                   geneColor = as.factor(both[["direction"]]),
                                   phenotypeColor = as.factor(both[["phenotype"]]))
      colored_lst[[name]] <- colored_gsc
    }

    up_ids <- up[[required_id]]
    message(up_ids)
    message("here?")
    if (is.null(up_ids)) {
      up_gsc <- NULL
    } else {
      up_gsc <- GSEABase::GeneSet(
                              identifier,
                              setName = up_name,
                              geneIds = as.character(up[[required_id]]))
      up_lst[[name]] <- up_gsc
    }

    down_gsc <- NULL
    down_ids <- down[[required_id]]
    if (!is.null(down_ids)) {
      down_lst[[name]] <- down_gsc
      if (!is.null(pair_names[2])) {
        down_name <- toupper(glue("{set_prefix}_{pair_names[2]}"))
        down_gsc <- GSEABase::GeneSet(
                                  identifier,
                                  setName = down_name,
                                  geneIds = as.character(down[[required_id]]))
        down_lst[[name]] <- down_gsc
        message("Created down_gsc")
      }
    }
  } ## End of the for loop.

  message("Got to retlist")
  retlst <- list(
      "colored" = colored_lst,
      "up" = up_lst,
      "down" = down_lst)
  return(retlst)
}

make_gsc_from_significant <- function(significant, according_to = "deseq", annotation_name = "org.Hs.eg.db",
                                      category_name = "infection",
                                      phenotype_name = "parasite", set_name = "elsayed_macrophage",
                                      color = TRUE, current_id = "ENSEMBL", required_id = "ENTREZID",
                                      min_gmt_genes = 10, ...) {
  ups <- list()
  downs <- list()
  if ("all" %in% according_to || length(according_to) > 1) {
    message("For now, limiting this to deseq.")
    according_to = "deseq"
  }
  ups <- significant[[according_to]][["ups"]]
  downs <- significant[[according_to]][["downs"]]

  set_prefix <- glue("{set_name}_{category_name}")
  color_set_name <- toupper(glue("{set_prefix}_{phenotype_name}"))
  up_name <- toupper(glue("{set_prefix}_up"))
  down_name <- toupper(glue("{set_prefix}_down"))

  ## Choose the Identifier for the colorsets.  For the moment just make it either entrez or null.
  identifier <- GSEABase::NullIdentifier()
  if (grepl(x=tolower(required_id), pattern="entrez")) {
    identifier <- GSEABase::EntrezIdentifier()
  } else if (grepl(x=tolower(required_id), pattern="ensemb")) {
    identifier <- GSEABase::ENSEMBLIdentifier()
  }

  contrasts <- names(ups)
  up_gscs <- list()
  down_gscs <- list()
  both_gscs <- list()
  up_ids <- list()
  down_ids <- list()
  both_ids <- list()
  for (contrast in contrasts) {
    up_ids[[contrast]] <- collect_gsc_ids_from_df(ups[[contrast]], annotation_name,
                                                  current_id, required_id)
    message("Number of up IDs in contrast ", contrast, ": ", length(up_ids[[contrast]]), ".")
    down_ids[[contrast]] <- collect_gsc_ids_from_df(downs[[contrast]], annotation_name,
                                                    current_id, required_id)
    message("Number of down IDs in contrast ", contrast, ": ", length(down_ids[[contrast]]), ".")
    both_factor_names <- c(up_ids[[contrast]], down_ids[[contrast]])
    both_factor_values <- c(rep("up", length(up_ids[[contrast]])),
                            rep("down", length(down_ids[[contrast]])))
    names(both_factor_values) <- both_factor_names
    both_ids[[contrast]] <- as.factor(both_factor_values)

    up_gscs[[contrast]] <- down_gscs[[contrast]] <- both_gscs[[contrast]] <- NULL
    if (!is.null(up_ids[[contrast]]) & length(up_ids[[contrast]]) >= min_gmt_genes) {
      up_name <- glue("{contrast}_up")
      up_gscs[[contrast]] <- GSEABase::GeneSet(identifier, setName = up_name,
                                               geneIds = as.character(up_ids[[contrast]]))
    }
    if (!is.null(down_ids[[contrast]]) & length(up_ids[[contrast]]) >= min_gmt_genes) {
      down_name <- glue("{contrast}_down")
      down_gscs[[contrast]] <- GSEABase::GeneSet(identifier, setName = down_name,
                                                 geneIds = as.character(down_ids[[contrast]]))
    }

    if (!is.null(both_ids[[contrast]])) {
      both_name <- glue("{contrast}_both")
      both_elements <- names(both_ids[[contrast]])
      dup_elements <- duplicated(both_elements)
      if (sum(dup_elements) > 0) {
        warning("There are ", sum(dup_elements),
                " non-unique elements in the IDs of the shared.")
        unique_names <- both_elements[!dup_elements]
        both_ids[[contrast]] <- both_ids[[contrast]][unique_names]
      }

      if (length(both_ids[[contrast]]) >= min_gmt_genes) {
        both_names <- names(both_ids[[contrast]])
        both_gscs[[contrast]] <- GSEABase::GeneColorSet(
                                               identifier,
                                               setName = both_name,
                                               geneIds = as.character(both_names),
                                               phenotype = phenotype_name,
                                               geneColor = both_ids[[contrast]])
      }
    }
  } ## End looking at every contrast.
  retlist <- list(
      "ups" = up_gscs,
      "downs" = down_gscs,
      "both" = both_gscs)
  return(retlist)
}

collect_gsc_ids_from_df <- function(df, orgdb, current_id, required_id) {
  ret_gsc <- NULL
  if (nrow(df) == 0) {
    return(NULL)
  }
  if (current_id == required_id) {
    return(rownames(df))
  }

  ## The rownames() of the expressionset must be in ENTREZIDs for gsva to work.
  ## tt <- sm(library(orgdb, character.only = TRUE))
  lib_result <- sm(requireNamespace(orgdb))
  att_result <- sm(try(attachNamespace(orgdb), silent = TRUE))

  ## Check that the current_id and required_id are in the orgdb.
  pkg <- get0(orgdb)
  available <- AnnotationDbi::columns(pkg)
  if (! current_id %in% available) {
    warning("The column: ", current_id, " is not in the set of available orgdb columns.")
    message("Here are the available columns:")
    print(available)
    message("Arbitrarily choosing 'GID'.")
    current_id <- "GID"
  }

  if (! required_id %in% available) {
    warning("The column: ", required_id, " is not in the set of available orgdb columns.")
    message("Here are the available columns:")
    print(available)
    message("Arbitrarily choosing 'GID'.")
    required_id <- "GID"
  }

  mesg("Converting the rownames() of the df to ", required_id, ".")
  req_ids <- sm(AnnotationDbi::select(x = get0(orgdb),
                                      keys = rownames(df),
                                      keytype = current_id,
                                      columns = c(required_id)))
  idx <- complete.cases(req_ids)
  req_ids <- req_ids[idx, required_id]
  return(req_ids)
}



#' Given a pairwise result, make a gene set collection.
#'
#' If I want to play with gsva and friends, then I need GeneSetCollections!
#' Much like make_gsc_from_significant(), this function extract the genes deemed
#' 'abundant' and generates gene sets accordingly.
#'
#' @param pairwise A pairwise result, or combined de result, or extracted genes.
#' @param according_to When getting significant genes, use this method.
#' @param annotation_name Name of the annotation dataset.
#' @param researcher_name Prefix of the name for the generated set(s).
#' @param study_name Second element in the name of the generated set(s).
#' @param category_name Third element in the name of the generated set(s).
#' @param phenotype_name Optional phenotype data for the generated set(s).
#' @param pair_names The suffix of the generated set(s).
#' @param current_id What type of ID is the data currently using?
#' @param required_id What type of ID should the use?
#' @param ... Extra arguments for extract_abundant_genes().
#' @return List containing 3 GSCs, one containing both the highs/lows called
#'  'colored', one of the highs, and one of the lows.
#' @seealso [extract_abundant_genes()] [make_gsc_from_ids()] [GSEABase]
#' @export
make_gsc_from_abundant <- function(pairwise, according_to = "deseq", annotation_name = "org.Hs.eg.db",
                                   researcher_name = "elsayed", study_name = "macrophage",
                                   category_name = "infection", phenotype_name = NULL,
                                   pair_names = "high", current_id = "ENSEMBL",
                                   required_id = "ENTREZID", ...) {
  highs <- list()
  lows <- list()
  if (class(pairwise)[1] == "data.frame") {
    highs <- pairwise
  } else if (class(pairwise)[1] == "all_pairwise") {
    message("Invoking combine_de_tables().")
    combined <- sm(combine_de_tables(pairwise, ...))
    message("Invoking extract_significant_genes().")
    highs <- sm(extract_abundant_genes(
        combined, according_to = according_to, ...)[[according_to]])
    lows <- sm(extract_abundant_genes(
        combined, according_to = according_to, least = TRUE, ...)[[according_to]])
  } else if (class(pairwise)[1] == "combined_de") {
    message("Invoking extract_significant_genes().")
    highs <- sm(extract_abundant_genes(
        pairwise, according_to = according_to, ...)[[according_to]])
    lows <- sm(extract_abundant_genes(
        pairwise, according_to = according_to, least = TRUE, ...)[[according_to]])
  } else if (class(pairwise)[1] == "character") {
    message("Invoking make_gsc_from_ids().")
    ret <- make_gsc_from_ids(pairwise, annotation_name = annotation_name,
                             pair_names = pair_names, category_name = category_name,
                             phenotype_name = phenotype_name,
                             current_id = current_id, required_id = required_id, ...)
    return(ret)
  } else {
    stop("I do not understand this type of input.")
  }

  ## The rownames() of the expressionset must be in ENTREZIDs for gsva to work.
  ## tt <- sm(library(annotation_name, character.only = TRUE))
  lib_result <- sm(requireNamespace(annotation_name))
  att_result <- sm(try(attachNamespace(annotation_name), silent = TRUE))
  high_lst <- list()
  low_lst <- list()
  colored_lst <- list()
  for (c in seq_len(high_lst[["abundances"]])) {
    name <- names(high_lst[["abundances"]])[c]
    high <- high_lst[["abundances"]][[name]]
    low <- low_lst[["abundances"]][[name]]
    if (class(high) == "character") {
      high_ids <- high
      low_ids <- low
    } else if (class(high) == "data.frame") {
      high_ids <- rownames(high)
      low_ids <- rownames(low)
    }
    if (current_id == required_id) {
      high[[required_id]] <- rownames(high)
      low[[required_id]] <- rownames(low)
    } else {
      message("Converting the rownames() of the expressionset to ENTREZID.")
      high_ids <- sm(AnnotationDbi::select(x = get0(annotation_name),
                                           keys = high_ids,
                                           keytype = current_id,
                                           columns = c(required_id)))
      high_idx <- complete.cases(high_ids)
      high_ids <- high_ids[high_idx, ]
      high <- merge(high, high_ids, by.x = "row.names", by.y = current_id)
      if (!is.null(low_ids)) {
        low_ids <- sm(AnnotationDbi::select(x = get0(annotation_name),
                                            keys = low_ids,
                                            keytype = current_id,
                                            columns = c(required_id)))
        low_idx <- complete.cases(low_ids)
        low_ids <- low_ids[low_idx, ]
        low <- merge(low, low_ids, by.x = "row.names", by.y = current_id)
      } else {
        low <- NULL
      }
    }
    high[["direction"]] <- pair_names[1]
    high[["phenotype"]] <- phenotype_name
    low[["direction"]] <- pair_names[2]
    low[["phenotype"]] <- phenotype_name

    both <- rbind(high, low)
    shared_ids <- high[[required_id]] %in% low[[required_id]]
    if (sum(shared_ids) > 0) {
      warning("There are ", sum(shared_ids), " shared IDs in the high and low lists.")
      shared_id <- high[shared_ids, ][[required_id]]
      shared_rows <- both[[required_id]] == shared_id
      both <- both[!shared_rows, ]
    }

    duplicated_elements <- duplicated(high[[required_id]])
    if (sum(duplicated_elements) > 0) {
      warning("There are ", sum(duplicated_elements),
              " non-unique elements in the IDs of the high data.")
      high <- high[!duplicated_elements, ]
    }
    duplicated_elements <- duplicated(low[[required_id]])
    if (sum(duplicated_elements) > 0) {
      warning("There are ", sum(duplicated_elements),
              " non-unique elements in the IDs of the low data.")
      low <- low[!duplicated_elements, ]
    }
    duplicated_elements <- duplicated(both[[required_id]])
    if (sum(duplicated_elements) > 0) {
      warning("There are ", sum(duplicated_elements),
              " non-unique elements in the IDs of the shared.")
      both <- both[!duplicated_elements, ]
    }

    set_prefix <- glue("{set_name}_{category_name}")
    color_set_name <- toupper(glue("{set_prefix}_{phenotype_name}"))
    high_name <- toupper(glue("{set_prefix}_{pair_names[1]}"))
    colored_gsc <- GSEABase::GeneColorSet(
                                 GSEABase::EntrezIdentifier(),
                                 setName = color_set_name,
                                 geneIds = as.character(both[[required_id]]),
                                 phenotype = phenotype_name,
                                 geneColor = as.factor(both[["direction"]]),
                                 phenotypeColor = as.factor(both[["phenotype"]]))
    colored_lst[[name]] <- colored_gsc
    high_gsc <- GSEABase::GeneSet(
                              GSEABase::EntrezIdentifier(),
                              setName = high_name,
                              geneIds = as.character(high[[required_id]]))
    high_lst[[name]] <- high_gsc
    low_gsc <- NULL
    low_lst[[name]] <- low_gsc
    if (!is.null(pair_names[2])) {
      low_name <- toupper(glue("{set_prefix}_{pair_names[2]}"))
      low_gsc <- GSEABase::GeneSet(
                               GSEABase::EntrezIdentifier(),
                               setName = low_name,
                               geneIds = as.character(low[[required_id]]))
      low_lst[[name]] <- low_gsc
    }
  } ## End of the for loop.

  retlst <- list(
      "colored" = colored_lst,
      "high" = high_lst,
      "low" = low_lst)
  return(retlst)
}

## EOF
