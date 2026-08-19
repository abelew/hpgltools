## annotation_sql.R: Load annotation information directly from SQL databases.

## Recently I have been having a hard time getting information from
## biomart.  I am apparently an old crank, and so sort of assume the
## fault lies in a combination of the scourge of LLMs and the effects
## of dear leader.
## I am sure there is a more rational explanation, but the result is
## the unfortunate fact that relatively recently it has become much
## more likely that my attempts to load a table of annotation
## information from biomart will fail.

## I therefore dusted off my SQL hat and started poking around on
## ftp.ensembl.org/pub/mysql/
## I figured out how one may populate a local database instance with
## the tables found there and put that code in CYOA in the file
## 'Downloaders.pm'

## The functions in this file are sisters to that, they provide some
## simple select statements so that I may collect arbitrary
## information from the biomart tables.  I do not yet fully understand
## their organization, so these will almost certainly be rudimentary
## until I finish expanding my understanding.

## The good news: a table that takes 3 minutes < x <= never-finished
## when querying useast.ensembl.org took < 1 second here.

#' Translate 'normal' column IDs to whatever the heck is going on in
#' the biomart schema.
#'
#' There are a lot of biomart tables (53,100), almost all of them are
#' full of alignment information used to infer gene homology
#' relationships.  The remaining tables are the fun stuff and fall
#' into a few categories: gene, transcript, translation.
#'
#' The tables also have a peculiar naming scheme which I have not yet
#' figured out; the end result is that if one wants the gene_id column
#' for mmusculus, one does not invoke
#' >'SELECT gene_id from mmusculus_gene'
#' but instead
#' 'SELECT gene_id_1020_key from mmusculus_gene_ensembl__gene__main'
#' I am sure there are good reasons for this, but I do not know them.
#' So, this function exists to take 'normal' names, invoke describe
#' table, and use some regex shenanigans to find the weirdo names.
#'
#' @param wanted List of desired tables and columns
#' @param species Species to query
#' @param dbh Existing database handle to the local (or not I guess)
#'  db.
#' @return A new list with the weirdo table and column names.
make_biomart_wanted <- function(wanted = NULL, species = "mmusculus", dbh = NULL) {
  if (is.null(wanted)) {
    wanted <- list(
      "gene__main" = c("gene_id", "stable_id", "version", "biotype", "display_label",
                       "description", "seq_region_start", "seq_region_end", "seq_region_strand"),
      "transcript__main" = c("gene_id", "stable_id_1066", "version", "biotype", "db_display_name",
                             "description", "transcript_count", "display_label",
                             "gene__main_stable_id_version", "transcript__main_stable_id_version")
     )
  }

  actual_wanted <- list()
  for (table in names(wanted)) {
    message("Getting columns from ", table, ".")
    actual_table <- glue("{species}_gene_ensembl__{table}")
    actual_wanted[[actual_table]] <- c()
    column_query <- DBI::dbSendQuery(dbh, glue("DESCRIBE {actual_table}"))
    description <- DBI::dbFetch(column_query)
    cleared <- DBI::dbClearResult(column_query)
    simple_columns <- wanted[[table]]
    actual_columns <- c()
    for (simple_column in simple_columns) {
      actual_column_idx <- grepl(x = description[["Field"]], pattern = simple_column)
      actual_column_name <- description[actual_column_idx, "Field"][1]
      actual_columns <- c(actual_columns, actual_column_name)
    }
    actual_wanted[[actual_table]] <- actual_columns
  }
  retlist <- list(
    "original" = wanted,
    "actual" = actual_wanted)
  return(retlist)
}

#' Load biomart data from a MariaDB copy.
#'
#' Given some tables and columns, select out all rows in the
#' appropriate biomart tables.  I use this to create species
#' annotation tables for RNASeq and similar tasks; as a result I
#' pretty much only ever want everything.  Once I need subsets I will
#' fill in the 'WHERE' clauses.
#'
#' @param species The tables are prefixed by species, so pick one.
#' @param overwrite This can write out a local rda file of the
#'  results, overwrite it?
#' @param do_save Save the output to a local rda file?
#' @param db_name Name of the local database.
#' @param dh_host Where to send queries.
#' @param db_user User at the database.
#' @param db_port This assumes MariaDB's 3306, I was thinking to use
#'  postgres but did not know if their schemas had mysql-specific
#'  stuff.  I may try messing with this in the future.
#' @param db_type Type of DBI connection to create.
#' @param db_password If you want to mess with the data, change the
#'  user and this.
#' @param wanted List of tables and columns of interest.
#' @param savefile Filename to which to save the results.
#' @return List with the same organization as 'wanted' containing
#'  dataframes with the columns specified in 'wanted'
load_biomart_sql_annotations <- function(species = "hsapiens", overwrite = FALSE, do_save = TRUE,
                                         db_name = "ensembl_mart", db_host = "localhost",
                                         db_user = "anonymous",
                                         db_port = 3306, db_type = "MariaDB",
                                         db_password = NULL, wanted = NULL, savefile = NULL) {
  ## TODO: I think the column names will be the weirdo column names
  ## from the biomart tables; they should be changed to those
  ## specified by 'wanted'.
  sql_annotations <- NULL
  if (file.exists(savefile) && isFALSE(overwrite)) {
    fresh <- new.env()
    message("The biomart annotations file already exists, loading from it.")
    ## load_string <- paste0("load('", savefile, "', envir = fresh)")
    load_string <- glue("load('{savefile}', envir = fresh)")
    eval(parse(text = load_string))
    biomart_annotations <- fresh[["sql_annotations"]]
    retlist <- list(
      "annotation" = sql_annotations,
      "host" = db_host,
      "database" = db_name,
      "species" = species)
    class(retlist) <- "hpgltools::load_biomart_sql_annotations"
    return(retlist)
  }

  dbh <- DBI::dbConnect(drv = RMariaDB::MariaDB(), username = db_user,
                        password = db_password, host = db_host,
                        dbname = db_name, port = db_port)
  if (is.null(wanted)) {
    wanted <- make_biomart_wanted(species = species, dbh = dbh)
  }

  original <- wanted[["original"]]
  actual <- wanted[["actual"]]

  sql_annotations <- list()
  for (table_idx in seq_along(names(wanted))) {
    wanted_name <- names(original)[table_idx]
    actual_name <- names(actual)[table_idx]
    wanted_columns <- original[[wanted_name]]
    actual_columns <- actual[[actual_name]]
    actual_string <- toString(actual_columns)
    query_result <- DBI::dbSendQuery(dbh, glue("SELECT {actual_string} from {actual_name}"))
    new_df <- DBI::dbFetch(query_result)
    colnames(new_df) <- wanted_columns
    sql_annotations[[wanted_name]] <- new_df
    cleared <- DBI::dbClearResult(query_result)
  }

  if (is.null(savefile)) {
    savefile <- glue("{species}_biomart_sql_annotations.rda")
  }

  if (isTRUE(do_save)) {
    message("Saving annotations to ", savefile, ".")
    save(list = ls(pattern = "sql_annotations"), file = savefile)
    message("Finished save().")
  }

  retlist <- list(
    "annotation" = sql_annotations,
    "host" = db_host,
    "database" = db_name)
  class(retlist) <- c("hpgltools::load_biomart_sql_annotations", "list")
  return(retlist)
}

#' I guess if I were smart I would do a join in SQL; but I am still
#' learning the structure of these tables.
merge_biomart_tables <- function(sql_annotations, primary_key = "gene_id") {
  final_df <- sql_annotations[["gene__main"]]
  additional_df <- sql_annotations[["transcript__main"]]
  final_df <- merge(final_df, additional_df, by = "gene_id", all.y = TRUE)
  return(final_df)
}
