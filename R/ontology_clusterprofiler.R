## ontology_clusterprofiler.r: Methods to simplify using clusterProfiler.  I
## think clusterprofiler is probably the most complete and comprehensive GSEA
## toolkit.  It is not necessarily as easy to use as I might prefer.  This seeks
## to fill in some corner case.

#' Run simple_clusterprofiler on every table from extract_significant_genes()
#'
#' Most of the options are taken from simple_cprofiler and passed directly down.
#'
#' @param sig Result from extract_significant_genes.
#' @param tables Result from combine_de_tables.
#' @param according_to Use this result type for the clusterprofiler searches.
#' @param together Concatenate the up/down genes into one set?
#' @param orgdb Name of the DBI containing the annotations for this organism.
#' @param orgdb_from Column from which to convert the gene IDs.
#' @param orgdb_to Column to which to convert, this must match the ontology IDS.
#' @param go_level Ignore categories above this level in the ontology tree.
#' @param pcutoff p-value significance cutoff.
#' @param qcutoff FDR adjusted significance cutoff.
#' @param fc_column Column containing the fold-change values.
#' @param second_fc_column A fallback column for FC values.
#' @param internal This changes the output data structure, but I forget how.
#' @param updown Seek out categories with increased enrichment,  or decreased.
#' @param permutations Run x permutations in clusterProfiler.
#' @param min_groupsize Ignore groups with less than x genes.
#' @param kegg_prefix Prefix of this kegg organism ID.
#' @param kegg_organism Full name of this organism when querying KEGG.
#' @param do_gsea Perform a full gene set enrichment analysis.
#' @param categories Plot this number of categories by default.
#' @param do_david Do a david over representation search? (the java
#'  DAVID interface is kind of broken, this should stay FALSE)
#' @param do_kegg Attempt a KEGG over representation analysis.
#' @param padj_type FDR correction method.
#' @param plot_type Choose specific plot(s).
#' @param do_reactome Attempt a reactome over representation analysis.
#' @param excel Output xlsx filename.
#' @param organism String name of the organism.
#' @param internal I dunno
#' @param max_groupsize Ignore groups which are too big.
#' @param padj_type Use this FDR
#' @param do_reactome what it says on the tin.
#' @param do_dose Attempt disease ontology search.
#' @param do_mesh Attempt MESH search.
#' @param do_msigdb Attempt mSigDB search.
#' @param mesh_category Use this category for MESH.
#' @param mesh_dbname Use this MESH sub-database.
#' @param msigdb_category Use this mSigDB sub-database.
#' @param msig_db Use this database file for the msigdb data.
#' @param ... Arguments to pass to simple_clusterprofiler().
#' @export
all_cprofiler <- function(sig, tables, according_to = "deseq", together = FALSE, todo = NULL,
                          orgdb = "org.Hs.eg.db", orgdb_from = NULL, orgdb_to = "ENTREZID",
                          go_level = 3, pcutoff = 0.05, qcutoff = 0.1, fc_column = "logFC",
                          second_fc_column = "deseq_logfc", internal = FALSE, updown = "up",
                          permutations = 1000, min_groupsize = 5, max_groupsize = 500,
                          kegg_prefix = NULL, kegg_organism = NULL, organism = "human",
                          categories = 12, padj_type = "BH",
                          mesh_category = "C", mesh_dbname = "gendoo",
                          msigdb_category = "C2", msig_db = NULL,
                          kegg_universe = NULL, reactome_organism = NULL,
                          mesh_db = NULL, ah_data = NULL, signature_data = NULL,
                          signature_df = NULL, de_table_namedf = NULL,
                          sig_genes_namedf = NULL,
                          excel = "excel/all_cp.xlsx", ...) {
  arglist <- list(...)
  ret <- list()
  input_up <- list()
  input_down <- list()
  source <- "significant"

  xlsx_dir <- dirname(excel)
  xlsx_base <- gsub(x = basename(excel), pattern = "\\.[[:alpha:]]{3,}$", replacement = "")
  ## Check if this came from extract_significant_genes or extract_abundant_genes.
  fc_col <- paste0(according_to, "_logfc")
  if (!is.null(sig[[according_to]][["ups"]])) {
    input_up <- sig[[according_to]][["ups"]]
    input_down <- sig[[according_to]][["downs"]]
  } else if (!is.null(sig[["abundances"]])) {
    source <- "abundance"
    input_up <- sig[["abundances"]][[according_to]][["high"]]
    input_down <- sig[["abundances"]][[according_to]][["low"]]
  } else {
    stop("I do not understand this input.")
  }

  sig_names <- names(input_up)
  for (i in seq_along(sig_names)) {
    name <- sig_names[i]
    table <- tables[["data"]][[name]]
    mesg("Starting ", name, ".")
    retname_up <- paste0(name, "_up")
    retname_down <- paste0(name, "_down")
    up <- input_up[[name]]
    down <- input_down[[name]]
    up_elements <- 0
    down_elements <- 0
    if (source == "abundance") {
      up <- names(up)
      down <- names(down)
      up_elements <- length(up)
      down_elements <- length(down)
    } else {
      up_elements <- nrow(up)
      down_elements <- nrow(down)
    }
    if (isTRUE(together)) {
      if (source == "abundance") {
        up <- c(up, down)
        up_elements <- up_elements + down_elements
        down <- c()
        down_elements <- 0
      } else {
        up <- rbind(up, down)
        up_elements <- nrow(up)
        down <- data.frame()
        down_elements <- 0
      }
    }
    if (up_elements > 0) {
      chosen_up_xlsx <- file.path(xlsx_dir, glue("{xlsx_base}_{retname_up}.xlsx"))
      up <- as.data.frame(up)
      table <- as.data.frame(table)
      args <- list(
        "sig_genes" = up, "de_table" = table, "orgdb" = orgdb, "todo" = todo,
        "orgdb_from" = orgdb_from, "orgdb_to" = orgdb_to, "go_level" = go_level,
        "pcutoff" = pcutoff, "qcutoff" = qcutoff, "fc_column" = fc_column,
        "second_fc_column" = second_fc_column, "internal" = internal, "updown" = updown,
        "permutations" = permutations, "min_groupsize" = min_groupsize, "kegg_prefix" = kegg_prefix,
        "kegg_organism" = kegg_organism, "categories" = categories, "padj_type" = padj_type,
        "excel" = chosen_up_xlsx, "organism" = organism,
        "max_groupsize" = max_groupsize, "mesh_category" = mesh_category,
        "mesh_dbname" = mesh_dbname, "msigdb_category" = msigdb_category, "msig_db" = msig_db,
        "kegg_universe" = kegg_universe, "reactome_organism" = reactome_organism,
        "mesh_db" = mesh_db, "ah_data" = ah_data, "signature_data" = signature_data,
        "signature_df" = signature_df, "de_table_namedf" = de_table_namedf,
        "sig_genes_namedf" = sig_genes_namedf)
      simple_cl <- try(do.call(what = "simple_clusterprofiler", args = args))
      if (! "try-error" %in% class(simple_cl)) {
        kegg_universe  <- simple_cl[["kegg_universe"]]
        reactome_organism <- simple_cl[["reactome_organism"]]
        mesh_db <- simple_cl[["mesh_db"]]
        ah_data <- simple_cl[["ah_data"]]
        signature_data <- simple_cl[["signature_data"]]
        signature_df <- simple_cl[["signature_df"]]
        de_table_namedf <- simple_cl[["de_table_namedf"]]
        sig_genes_namedf <- simple_cl[["sig_genes_namedf"]]
        orgdb_from <- simple_cl[["orgdb_from"]]
        orgdb_to <- simple_cl[["orgdb_to"]]
        todo <- simple_cl[["todo"]]
        ret[[retname_up]] <- simple_cl
      } else {
        ret[[retname_up]] <- NULL
      }
    } else {
      ret[[retname_up]] <- NULL
    }
    if (down_elements > 0) {
      slept <- Sys.sleep(10)
      chosen_down_xlsx <- file.path(xlsx_dir, glue("{xlsx_base}_{retname_down}.xlsx"))
      args <- list(
        "sig_genes" = down, "de_table" = table, "orgdb" = orgdb, "todo" = todo,
        "orgdb_from" = orgdb_from, "orgdb_to" = orgdb_to, "go_level" = go_level,
        "pcutoff" = pcutoff, "qcutoff" = qcutoff, "fc_column" = fc_column,
        "second_fc_column" = second_fc_column, "internal" = internal, "updown" = updown,
        "permutations" = permutations, "min_groupsize" = min_groupsize, "kegg_prefix" = kegg_prefix,
        "kegg_organism" = kegg_organism, "categories" = categories, "padj_type" = padj_type,
        "excel" = chosen_up_xlsx, "organism" = organism,
        "max_groupsize" = max_groupsize, "mesh_category" = mesh_category,
        "mesh_dbname" = mesh_dbname, "msigdb_category" = msigdb_category, "msig_db" = msig_db,
        "kegg_universe" = kegg_universe, "reactome_organism" = reactome_organism,
        "mesh_db" = mesh_db, "ah_data" = ah_data, "signature_data" = signature_data,
        "signature_df" = signature_df, "de_table_namedf" = de_table_namedf,
        "sig_genes_namedf" = sig_genes_namedf)
      simple_cl <- try(do.call(what = "simple_clusterprofiler", args = args))
      if (! "try-error" %in% simple_cl) {
        kegg_universe  <- simple_cl[["kegg_universe"]]
        reactome_organism <- simple_cl[["reactome_organism"]]
        mesh_db <- simple_cl[["mesh_db"]]
        ah_data <- simple_cl[["ah_data"]]
        signature_data <- simple_cl[["signature_data"]]
        signature_df <- simple_cl[["signature_df"]]
        de_table_namedf <- simple_cl[["de_table_namedf"]]
        sig_genes_namedf <- simple_cl[["sig_genes_namedf"]]
        orgdb_from <- simple_cl[["orgdb_from"]]
        orgdb_to <- simple_cl[["orgdb_to"]]
        todo <- simple_cl[["todo"]]
        ret[[retname_down]] <- simple_cl
      } else {
        ret[[retname_down]] <- NULL
      }
      #ret[[retname_down]] <- sm(simple_clusterprofiler(down, first_col = fc_col))
    } else {
      ret[[retname_down]] <- NULL
    }
  }
  class(ret) <- "all_cprofiler"
  return(ret)
}

#' Use a heuristic to guess the appropriate keytype when using clusterProfiler::enricher().
#'
#' @param org orgdb containing the potential keys
#' @param from starting keytype
#' @param sig_genes Input gene set
#' @param to new keytype
#' @param possible_keys Set of keys to test
#' @param universe Universe of possible genes.
guess_bitr_keytype <- function(org, from, sig_genes = NULL, to = "ENTREZ",
                               possible_keys = NULL, universe = NULL) {
  if (is.null(universe)) {
    universe <- AnnotationDbi::keys(org)
  }
  if (is.null(possible_keys)) {
    possible_keys <- AnnotationDbi::keytypes(org)
  }
  from <- NULL
  chosen_gene_df <- data.frame()
  chosen_sig_df <- data.frame()
  test_genes_df <- data.frame()
  test_sig_df <- data.frame()
  num_sig <- 0
  num_hits <- 0
  orgdb_sig_from <- from
  for (k in seq_along(possible_keys)) {
    key <- possible_keys[k]
    test_genes_df <- sm(try(clusterProfiler::bitr(universe, fromType = key,
                                                  toType = to, OrgDb = org), silent = TRUE))
    if (class(test_genes_df) == "try-error") {
      test_genes_df <- data.frame()
    }
    if (!is.null(sig_genes)) {
      test_sig_df <- sm(try(clusterProfiler::bitr(sig_genes, fromType = key,
                                                  toType = to, OrgDb = org), silent = TRUE))
      if (class(test_sig_df) == "try-error") {
        test_sig_df <- data.frame()
      }
    }
    test_num_hits <- nrow(test_genes_df)
    if (test_num_hits > num_hits) {
      from <- key
      num_hits <- test_num_hits
      chosen_gene_df <- test_genes_df
    }
    test_sig_hits <- nrow(test_sig_df)
    if (test_sig_hits > num_sig) {
      orgdb_sig_from <- key
      num_sig <- test_sig_hits
      chosen_sig_df <- test_sig_df
    }
  }
  mesg("Chose keytype: ", from, " for all genes because it had ", num_hits,
       " out of ", length(universe), " genes.")
  if (!is.null(sig_genes)) {
    mesg("Chose keytype: ", orgdb_sig_from, " for sig genes because it had ", num_sig,
         " out of ", length(sig_genes), " genes.")
  }
  retlist <- list(
    "gene_df" = chosen_gene_df,
    "sig_df" = chosen_sig_df,
    "universe_from" = from,
    "sig_from" = orgdb_sig_from)
  return(retlist)
}

#' I cannot be trusted to type 'cluster'
#'
#' @param ... Passed to simple_clusterprofiler()
#' @export
simple_cprofiler <- function(...) {
  simple_clusterprofiler(...)
}

#' Perform the array of analyses in the 2016-04 version of clusterProfiler
#'
#' The new version of clusterProfiler has a bunch of new toys.  However, it is
#' more stringent in terms of input in that it now explicitly expects to receive
#' annotation data in terms of a orgdb object.  This is mostly advantageous, but
#' will probably cause some changes in the other ontology functions in the near
#' future.  This function is an initial pass at making something similar to my
#' previous 'simple_clusterprofiler()' but using these new toys.
#'
#' @param sig_genes Dataframe of genes deemed 'significant.'
#' @param de_table Dataframe of all genes in the analysis, primarily for GSEA.
#' @param orgdb Name of the orgDb used for gathering annotation data.
#' @param orgdb_from Name of a key in the orgdb used to cross reference to entrez IDs.
#' @param orgdb_to List of keys to grab from the orgdb for cross referencing
#'  ontologies.
#' @param go_level How deep into the ontology tree should this dive for over
#'  expressed categories.
#' @param pcutoff P-value cutoff for 'significant' analyses.
#' @param qcutoff Q-value cutoff for 'significant' analyses.
#' @param fc_column When extracting vectors of all genes, what column should be used?
#' @param second_fc_column When extracting vectors of all genes, what column
#'  should be tried the second time around?
#' @param updown Include the less than expected ontologies?
#' @param permutations How many permutations for GSEA-ish analyses?
#' @param min_groupsize Minimum size of an ontology before it is included.
#' @param kegg_prefix Many KEGG ids need a prefix before they will cross reference.
#' @param kegg_organism Choose the 3 letter KEGG organism name here.
#' @param do_gsea Perform gsea searches?
#' @param categories How many categories should be plotted in bar/dot plots?
#' @param excel Print the results to an excel file?
#' @param do_david Attempt to use the DAVID database for a search?
#' @param do_kegg Perform kegg search?
#' @param david_id Which column to use for cross-referencing to DAVID?
#' @param david_user Default registered username to use.
#' @param organism String name of the organism.
#' @param internal I dunno
#' @param max_groupsize Ignore groups which are too big.
#' @param padj_type Use this FDR
#' @param do_reactome what it says on the tin.
#' @param do_dose Attempt disease ontology search.
#' @param do_mesh Attempt MESH search.
#' @param do_msigdb Attempt mSigDB search.
#' @param mesh_category Use this category for MESH.
#' @param mesh_dbname Use this MESH sub-database.
#' @param msigdb_category Use this mSigDB sub-database.
#' @param msig_db Use this database file for the msigdb data.
#' @return a list
#' @seealso [clusterProfiler] [AnnotationDbi] [KEGGREST]
#' @examples
#' \dontrun{
#'  holyasscrackers <- simple_clusterprofiler(gene_list, all_genes, "org.Dm.eg.db")
#' }
#' @export
simple_clusterprofiler <- function(sig_genes, de_table = NULL, orgdb = "org.Hs.eg.db",
                                   todo = NULL, orgdb_from = NULL, orgdb_to = "ENTREZID",
                                   go_level = 3, pcutoff = 0.05, organism = "human", qcutoff = 0.1,
                                   fc_column = "logFC", second_fc_column = "deseq_logfc",
                                   internal = FALSE, updown = "up", permutations = 1000,
                                   min_groupsize = 5, max_groupsize = 500, kegg_prefix = NULL,
                                   kegg_organism = NULL, categories = 12, excel = NULL,
                                   david_id = "ENTREZ_GENE_ID", padj_type = "BH",
                                   david_user = "abelew@umd.edu", mesh_category = "C",
                                   mesh_dbname = "gendoo", msigdb_category = "C2", msig_db = NULL,
                                   kegg_universe = NULL, reactome_organism = NULL,
                                   mesh_db = NULL, ah_data = NULL, signature_data = NULL,
                                   signature_df = NULL, de_table_namedf = NULL,
                                   sig_genes_namedf = NULL) {
  loaded <- sm(requireNamespace(package = "clusterProfiler", quietly = TRUE))
  loaded <- sm(requireNamespace(package = "DOSE", quietly = TRUE))
  org <- NULL
  ## Start off by figuring out what was given, an OrgDb or the name of one.
  ## S4ME!
  if (is.null(todo)) {
    todo <- list(
      "enrich_go" = TRUE,
      "gse_go" = TRUE,
      "enrich_kegg" = TRUE,
      "gse_kegg" = TRUE,
      "enrich_reactome" = TRUE,
      "gse_reactome" = TRUE,
      "enrich_david" = FALSE,
      "gse_david" = FALSE,
      "enrich_dose" = TRUE,
      "gse_dose" = TRUE,
      "enrich_mesh" = TRUE,
      "gse_mesh" = TRUE,
      "enrich_msigdb" = TRUE,
      "gse_msigdb" = TRUE)
  }
  if (class(orgdb)[[1]] == "OrgDb") {
    org <- orgdb
  } else if ("character" %in% class(orgdb)) {
    sm(requireNamespace(orgdb))
    org <- loadNamespace(orgdb) ## put the orgDb instance into an environment
    org <- org[[orgdb]] ## Then extract it
  } else {
    stop("Need either the name of an orgdb package or the orgdb itself.")
  }

  ## It is likely that we will need to query multiple different keys from the OrgDb
  ## So gather them now for later reference.
  mapper_keys <- AnnotationDbi::keytypes(org)
  ## If we must, we can extract the set of all genes from the orgdb
  ## However, this means we may not do a GSEA analysis
  universe_genes <- AnnotationDbi::keys(org)
  if (is.null(de_table)) {
    do_gsea <- FALSE
    all_genenames <- universe_genes
  } else {
    all_genenames <- rownames(de_table)
    shared_genes <- all_genenames %in% universe_genes
  }
  sig_genenames <- rownames(sig_genes)
  orgdb_to <- toupper(orgdb_to)

  test_genes_df <- NULL
  test_sig_df <- NULL
  num_sig <- 0
  num_hits <- 0
  orgdb_sig_from <- orgdb_from
  ## I think this logic should be split into a separate function and used elsewhere
  if (is.null(orgdb_from)) {
    mesg("Testing available OrgDb keytypes for the best mapping to entrez.")
    for (k in mapper_keys) {
      test_genes_df <- sm(try(clusterProfiler::bitr(
        all_genenames, fromType = k, toType = orgdb_to, OrgDb = org), silent = TRUE))
      test_sig_df <- sm(try(clusterProfiler::bitr(
        sig_genenames, fromType = k, toType = orgdb_to, OrgDb = org), silent = TRUE))
      if (class(test_genes_df) == "try-error") {
        test_genes_df <- data.frame()
      }
      if (class(test_sig_df) == "try-error") {
        test_sig_df <- data.frame()
      }
      test_num_hits <- nrow(test_genes_df)
      if (test_num_hits > num_hits) {
        orgdb_from <- k
        num_hits <- test_num_hits
        de_table_namedf <- test_genes_df
      }
      test_sig_hits <- nrow(test_sig_df)
      if (test_sig_hits > num_sig) {
        orgdb_sig_from <- k
        num_sig <- test_sig_hits
        sig_genes_namedf <- test_sig_df
      }
    }
    mesg("Chose keytype: ", orgdb_from, " for all genes because it had ", num_hits,
         " out of ", length(all_genenames), " genes.")
    mesg("Chose keytype: ", orgdb_sig_from, " for sig genes because it had ", num_sig,
         " out of ", length(sig_genenames), " genes.")
    ## I want to replace this blob with the following function.
    ##comparison <- guess_bitr_keytype(org, orgdb_from, sig_genes = sig_genenames, to = orgdb_to,
    ##                                 possible_keys = mapper_keys, universe = all_genenames)

  } else if (orgdb_from == orgdb_to) {
    de_table_namedf <- data.frame(from = all_genenames, to = all_genenames)
    sig_genes_namedf <- data.frame(from = sig_genenames, to = sig_genenames)
    colnames(de_table_namedf) <- c(orgdb_to, "same")
    colnames(sig_genes_namedf) <- c(orgdb_to, "same")
  } else { ## If we do have a column for the OrgDB
    de_table_namedf <- sm(try(clusterProfiler::bitr(
      all_genenames, fromType = orgdb_from, toType = orgdb_to, OrgDb = org), silent = TRUE))
    sig_genes_namedf <- sm(try(clusterProfiler::bitr(
      sig_genenames, fromType = orgdb_from, toType = orgdb_to, OrgDb = org), silent = TRUE))
  }

  if (is.null(sig_genes[[fc_column]]) && is.null(sig_genes[[second_fc_column]])) {
    stop("The fold change column provided no genes, try another column in the data set.")
  } else if (is.null(sig_genes[[fc_column]])) {
    fc_column <- second_fc_column
  }
  gsea_fc_column <- fc_column
  if (is.null(de_table[[gsea_fc_column]]) && is.null(de_table[[second_fc_column]])) {
    message("Unable to find the fold-change column in the de table, not doing gsea.")
    todo[["gse_go"]] <- FALSE
    todo[["gse_kegg"]] <- FALSE
    todo[["gse_reactome"]] <- FALSE
    todo[["gse_david"]] <- FALSE
    todo[["gse_dose"]] <- FALSE
    todo[["gse_mesh"]] <- FALSE
    todo[["gse_msigdb"]] <- FALSE
  } else if (is.null(de_table[[gsea_fc_column]])) {
    gsea_fc_column <- second_fc_column
  }
  ## Acquire the set of IDs against which all queries need to be made
  universe_to <- AnnotationDbi::keys(org, keytype = orgdb_to)
  ## And the set of similar IDs mapped against the significance table.
  all_gene_list <- de_table_namedf[[orgdb_to]]
  all_gene_drop <- !is.na(all_gene_list)
  all_gene_list <- all_gene_list[all_gene_drop]
  sig_gene_list <- sig_genes_namedf[[orgdb_to]]
  sig_gene_drop <- !is.na(sig_gene_list)
  sig_gene_list <- sig_gene_list[sig_gene_drop]
  if (is.null(sig_gene_list)) {
    stop("No genes were found between the significant genes and the universe.")
  }
  ## Now we have a universe of geneIDs and significant IDs
  ## Let us perform some analyses...
  mesg("Calculating GO groups.")
  ggo_mf <- ggo_bp <- ggo_cc <- NULL ## GO groups by ontology
  ego_all_mf <- ego_all_bp <- ego_all_cc <- NULL ## GO enrichment data structures
  ego_all_mf_df <- ego_all_bp_df <- ego_all_cc_df <- NULL ## GO enrichment full dfs
  ego_sig_mf_df <- ego_sig_bp_df <- ego_sig_cc_df <- NULL ## pvalue cutoff dfs
  if (isTRUE(todo[["enrich_go"]])) {
    ggo_mf <- sm(clusterProfiler::groupGO(
      gene = sig_gene_list, OrgDb = org, keyType = orgdb_to, ont = "MF", level = go_level))
    ggo_bp <- sm(clusterProfiler::groupGO(
      gene = sig_gene_list, OrgDb = org, keyType = orgdb_to, ont = "BP", level = go_level))
    ggo_cc <- sm(clusterProfiler::groupGO(
      gene = sig_gene_list, OrgDb = org, keyType = orgdb_to, ont = "CC", level = go_level))
    ego_all_mf <- clusterProfiler::enrichGO(
      gene = sig_gene_list, universe = universe_to, OrgDb = org, ont = "MF", keyType = orgdb_to,
      minGSSize = min_groupsize, pAdjustMethod = padj_type, pvalueCutoff = 1.0)
    mesg("Performing GO over-representation analyses.")
    ego_all_mf_df <- as.data.frame(ego_all_mf, stringsAsFactors = FALSE)
    sig_idx <- ego_all_mf_df[["p.adjust"]] <= pcutoff
    ego_sig_mf_df <- ego_all_mf_df[sig_idx, ]
    ego_all_bp <- clusterProfiler::enrichGO(
      gene = sig_gene_list, universe = universe_to, OrgDb = org, ont = "BP", keyType = orgdb_to,
      minGSSize = min_groupsize, pAdjustMethod = padj_type, pvalueCutoff = 1.0)
    ego_all_bp_df <- as.data.frame(ego_all_bp, stringsAsFactors = FALSE)
    sig_idx <- ego_all_bp_df[["p.adjust"]] <= pcutoff
    ego_sig_bp_df <- ego_all_bp_df[sig_idx, ]
    ego_all_cc <- clusterProfiler::enrichGO(
      gene = sig_gene_list, universe = universe_to, OrgDb = org, ont = "CC", keyType = orgdb_to,
      minGSSize = min_groupsize, pAdjustMethod = padj_type, pvalueCutoff = 1.0)
    ego_all_cc_df <- as.data.frame(ego_all_cc, stringsAsFactors = FALSE)
    sig_idx <- ego_all_cc_df[["p.adjust"]] <= pcutoff
    ego_sig_cc_df <- ego_all_cc_df[sig_idx, ]
    mesg("Found ", nrow(ego_sig_mf_df),
         " MF, ", nrow(ego_sig_bp_df),
         " BP, and ", nrow(ego_sig_cc_df), " CC over-represented hits.")
  }
  de_table_merged <- NULL
  genelist <- vector()
  if (isTRUE(todo[["gse_go"]]) || isTRUE(todo[["gse_kegg"]]) ||
        isTRUE(todo[["gse_reactome"]]) || isTRUE(todo[["gse_dose"]]) ||
        isTRUE(todo[["gse_mesh"]]) || isTRUE(todo[["gse_msigdb"]])) {
    ## Why did I do this? ## Ahh for the GSE analyses, they want ordered gene IDs
    ## Add the entrezIDs to the end
    de_table_merged <- merge(de_table, de_table_namedf, by.x = "row.names", by.y = 1)
    if (updown == "up") {
      new_order <- order(de_table_merged[[gsea_fc_column]], decreasing = TRUE)
      de_table_merged <- de_table_merged[new_order, ]
    } else {
      new_order <- order(de_table_merged[[gsea_fc_column]], decreasing = FALSE)
      de_table_merged <- de_table_merged[new_order, ]
    }
    genelist <- as.vector(de_table_merged[[gsea_fc_column]])
    names(genelist) <- de_table_merged[[orgdb_to]]
    duplicated_names <- duplicated(names(genelist))
    num_duplicates <- sum(duplicated_names)
    if (num_duplicates > 0) {
      mesg("There are ", num_duplicates, " duplicated gene IDs, dropping them.")
      genelist <- genelist[!duplicated_names]
    }
  }
  gse_go <- list()
  gse_go_all_df <- gse_go_sig_df <- data.frame()
  if (isTRUE(todo[["gse_go"]])) {
    ## Hmm this is odd, in the previous calls, I used orgdb_to, but in this set
    ## I am using orgdb_from...
    mesg("Performing GO GSEA of gene lists (this is slow).")
    #else {
    #  genelist <- as.vector(sig_genes[[fc_column]])
    #  names(genelist) <- rownames(sig_genes)
    #}
    ## 2020 04: Adding a pvalue cutoff argument causes an error, I do not know why.
    ## Arguments used by gseGO of interest: exponent, minGSSize/maxGSSize, eps, by(fgsea)
    ## Also, apparently the nperm argument is deprecated.
    gse_go <- try(clusterProfiler::gseGO(
      geneList = genelist, OrgDb = org, keyType = orgdb_to, ont = "ALL", minGSSize = min_groupsize))
    if ("try-error" %in% class(gse_go)) {
      gse_go <- NULL
    } else {
      gse_go_all_df <- as.data.frame(gse_go)
      sig_idx <- gse_go_all_df[["p.adjust"]] <= pcutoff
      gse_go_sig_df <- gse_go_all_df[sig_idx, ]
      mesg("Found ", nrow(gse_go_sig_df), " GO GSE hits.")
    }
  }
  go_data <- list(
    "MF_enrich" = ego_all_mf,
    "BP_enrich" = ego_all_bp,
    "CC_enrich" = ego_all_cc,
    "GO_gse" = gse_go,
    "MF_all" = ego_all_mf_df,
    "MF_sig" = ego_sig_mf_df,
    "BP_all" = ego_all_bp_df,
    "BP_sig" = ego_sig_bp_df,
    "CC_all" = ego_all_cc_df,
    "CC_sig" = ego_sig_cc_df,
    "GO_gse_all" = gse_go_all_df,
    "GO_gse_sig" = gse_go_sig_df)
  all_kegg <- enrich_kegg <- NULL
  all_kegg_df <- sig_kegg_df <- data.frame()
  if (isTRUE(todo[["enrich_kegg"]])) {
    ## Set up to do kegg
    ## Now extract the kegg organism/gene IDs.
    if (is.null(kegg_organism)) {
      org_meta <- AnnotationDbi::metadata(org)
      org_row <- org_meta[["name"]] == "ORGANISM"
      organism <- org_meta[org_row, "value"]
      ## Only grab the first of potentially multiple outputs.
      kegg_organism <- get_kegg_orgn(species = organism)
    }
    if (length(kegg_organism) == 1) {
      mesg("Found organism ID: ", kegg_organism, ".")
    } else if (length(kegg_organism) > 0) {
      kegg_organism <- kegg_organism[[1]]
      if (is.null(kegg_universe)) {
        kegg_universe <- try(KEGGREST::keggConv(kegg_organism, "ncbi-geneid"))
      }
      kegg_sig_intersect <- 0
      kegg_sig_names <- NULL
      if (! "try-error" %in% class(kegg_universe)) {
        kegg_sig_names <- glue("ncbi-geneid:{sig_gene_list}")
        kegg_sig_intersect <- kegg_sig_names %in% names(kegg_universe)
      }
      mesg("Found ", sum(kegg_sig_intersect),
           " matches between the significant gene list and kegg universe.")
      if (sum(kegg_sig_intersect) > 0) {
        all_names <- names(kegg_universe)
        small_universe <- kegg_universe[intersect(kegg_sig_names, all_names)]
        kegg_sig_ids <- unique(as.character(small_universe))
        ##kegg_sig_ids <- unique(as.character(kegg_universe[kegg_sig_intersect]))
        kegg_sig_ids <- gsub(pattern = glue("{kegg_organism}:"),
                             replacement = "", x = kegg_sig_ids)
        mesg("Performing KEGG over-representation analysis.")
      }
      all_kegg <- try(clusterProfiler::enrichKEGG(
        kegg_sig_ids, organism = kegg_organism, keyType = "kegg", pvalueCutoff = 1.0))
      if ("try-error" %in% class(all_kegg)) {
        all_kegg <- NULL
        todo[["gse_kegg"]] <- FALSE
      } else {
        all_kegg_df <- as.data.frame(all_kegg, stringsAsFactors = FALSE)
        sig_idx <- all_kegg_df[["p.adjust"]] <= pcutoff
        sig_kegg_df <- all_kegg_df[sig_idx, ]
      }
    } else { ## We did not get a kegg organism
      todo[["enrich_kegg"]] <- FALSE
      todo[["gse_kegg"]] <- FALSE
      kegg_organism <- NULL
    }
  }
  gse_all_kegg <- gse_sig_kegg <- NULL
  gse_all_kegg_df <- gse_sig_kegg_df <- data.frame()
  if (isTRUE(todo[["gse_kegg"]])) {
    lastcol <- ncol(de_table_merged)
    kegg_genelist <- as.vector(de_table_merged[[fc_column]])
    names(kegg_genelist) <- de_table_merged[[lastcol]]
    kegg_all_names <- glue("ncbi-geneid:{names(kegg_genelist)}")
    kegg_all_intersect <- kegg_all_names %in% names(kegg_universe)
    mesg("Found ", sum(kegg_all_intersect),
         " matches between the gene list and kegg universe.")
    all_names <- names(kegg_universe)
    large_universe <- kegg_universe[intersect(kegg_all_names, names(kegg_universe))]
    kegg_all_ids <- unique(as.character(large_universe))
    kegg_all_ids <- gsub(pattern = glue("{kegg_organism}:"),
                         replacement = "", x = kegg_all_ids)
    names(kegg_genelist) <- kegg_all_ids
    ## Something changed in bioc 3.20 causing some NAs to creep in here.
    na_idx <- is.na(names(kegg_genelist))
    kegg_genelist <- kegg_genelist[!na_idx]
    mesg("Performing KEGG GSEA.")
    gse_all_kegg <- try(
      clusterProfiler::gseKEGG(
        geneList = kegg_genelist, organism = kegg_organism, nPerm = permutations,
        minGSSize = min_groupsize, pvalueCutoff = 1.0, use_internal_data = internal))
    if ("try-error" %in% class(gse_all_kegg)) {
      gse_all_kegg <- NULL
    } else {
      gse_all_kegg_df <- as.data.frame(gse_all_kegg, stringsAsFactors = FALSE)
      sig_idx <- gse_all_kegg_df[["p.adjust"]] <= pcutoff
      gse_sig_kegg_df <- gse_all_kegg_df[sig_idx, ]
      mesg("Found ", nrow(gse_sig_kegg_df), " KEGG GSE hits.")
    }
  }
  kegg_data <- list(
    "kegg_enrich" = all_kegg,
    "kegg_gse" = gse_all_kegg,
    "kegg_all" = all_kegg_df,
    "kegg_sig" = sig_kegg_df,
    "kegg_gse_all" = gse_all_kegg_df,
    "kegg_gse_sig" = gse_sig_kegg_df)
  david_all <- NULL
  david_all_df <- david_sig_df <- data.frame()
  if (isTRUE(todo[["enrich_david"]])) {
    mesg("Attempting DAVID search.")
    david_all <- try(clusterProfiler::enrichDAVID(
      gene = sig_gene_list, minGSSize = min_groupsize,
      idType = david_id, david.user = david_user), silent = TRUE)
    if ("try-error" %in% class(david_all)) {
      todo[["gse_david"]] <- FALSE
    } else {
      david_all_df <- as.data.frame(david_search, stringsAsFactors = FALSE)
      sig_idx <- david_all_df[["p.adjust"]] <= pcutoff
      david_sig_df <- david_all_df[sig_idx, ]
    }
  }

  gse_david <- NULL
  gse_david_all_df <- gse_david_sig_df <- data.frame()
  if (isTRUE(todo[["david_gse"]])) {
    gse_david <- sm(suppressWarnings(clusterProfiler::gseDAVID(
      geneList = genelist, OrgDb = org, keyType = orgdb_to, ont = "ALL", minGSSize = min_groupsize)))
    gse_david_all_df <- as.data.frame(gse_david)
    sig_idx <- gse_david_df[["p.adjust"]] <= pcutoff
    gse_david_sig_df <- gse_david_all_df[sig_idx, ]
  }
  david_data <- list(
    "david_enrich" = david_all,
    "david_all" = david_all_df,
    "david_sig" = david_sig_df,
    "david_gse" = gse_david,
    "david_gse_all" = gse_david_all_df,
    "david_gse_sig" = gse_david_sig_df)

  reactome_all <- NULL
  reacome_all_df <- reactome_sig_df <- data.frame()
  orgdb_name <- orgdb
  if ("OrgDb" %in% class(orgdb)) {
    orgdb_name <- orgdb@.xData[["packageName"]]
  }
  if (orgdb_name == "org.Hs.eg.db") {
    reactome_organism <- "human"
  } else if (orgdb_name == "org.Mm.eg.db") {
    reactome_organism <- "mouse"
  } else {
    todo[["enrich_reactome"]] <- FALSE
    todo[["gse_reactome"]] <- FALSE
  }

  reactome_all <- NULL
  reactome_all_df <- reactome_sig_df <- data.frame()
  if (isTRUE(todo[["enrich_reactome"]])) {
    loaded <- sm(requireNamespace(package = "ReactomePA", quietly = TRUE))
    mesg("Performing reactome over-representation analysis.")
    reactome_all <- try(ReactomePA::enrichPathway(
      gene = sig_gene_list, pvalueCutoff = pcutoff, readable = TRUE,
      pAdjustMethod = padj_type, qvalueCutoff = 1.0, universe = universe_to,
      organism = reactome_organism, minGSSize = min_groupsize, maxGSSize = max_groupsize))
    if ("try-error" %in% class(reactome_all)) {
      todo[["gse_reactome"]] <- FALSE
      reactome_all <- NULL
    } else {
      reactome_all_df <- as.data.frame(reactome_all, stringsAsFactors = FALSE)
      sig_idx <- reactome_all_df[["p.adjust"]] <= pcutoff
      reactome_sig_df <- reactome_all_df[sig_idx, ]
      mesg("Found ", nrow(reactome_sig_df), " reactome over-represented hits.")
    }
  }

  gse_reactome <- NULL
  gse_reactome_all_df <- gse_reactome_sig_df <- data.frame()
  if (isTRUE(todo[["gse_reactome"]])) {
    mesg("Performing reactome GSEA.")
    gse_all_reactome <- try(
      ReactomePA::gsePathway(geneList = genelist, organism = reactome_organism,
                             pvalueCutoff = 1.0, minGSSize = min_groupsize,
                             maxGSSize = max_groupsize))
    if ("try-error" %in% class(gse_all_reactome)) {
      gse_all_reactome <- NULL
    } else {
      gse_reactome_all_df <- as.data.frame(gse_all_reactome, stringsAsFactors = FALSE)
      sig_idx <- gse_reactome_all_df[["p.adjust"]] <= pcutoff
      gse_reactome_sig_df <- reactome_all_df[sig_idx, ]
      mesg("Found ", nrow(gse_reactome_sig_df), " Reactome GSE hits.")
    }
  }
  reactome_data <- list(
    "reactome_enrich" = reactome_all,
    "reactome_all" = reactome_all_df,
    "reactome_sig" = reactome_sig_df,
    "reactome_gse" = gse_reactome,
    "reactome_gse_all" = gse_reactome_all_df,
    "reactome_gse_sig" = gse_reactome_sig_df)

  dose_all <- NULL
  dose_all_df <- dose_sig_df <- data.frame()
  orgn <- "hsa"
  do_db <- "HDO"
  if (organism == "human" || organism == "Homo sapiens") {
    loaded <- sm(requireNamespace(package = "HDO.db", quietly = TRUE))
  } else if (organism == "mouse" || organism == "Mus musculus") {
    orgn <- "mm"
    do_db <- "MPO"
    loaded <- sm(requireNamespace(package = "MPO.db", quietly = TRUE))
  } else {
    todo[["enrich_dose"]] <- FALSE
    todo[["gse_dose"]] <- FALSE
    todo[["enrich_msigdb"]] <- FALSE
    todo[["gse_msigdb"]] <- FALSE
  }
  dose_all <- NULL
  dose_all_df <- dose_sig_df <- data.frame()
  if (isTRUE(todo[["enrich_dose"]])) {
    mesg("Performing DOSE over-representation analysis.")
    dose_all <- try(DOSE::enrichDO(
      gene = sig_gene_list, ont = do_db, organism = orgn,
      pvalueCutoff = pcutoff, pAdjustMethod = padj_type, universe = universe_to,
      minGSSize = min_groupsize, maxGSSize = max_groupsize,
      qvalueCutoff = 1.0, readable = TRUE))
    if ("try-error" %in% class(dose_all)) {
      todo[["gse_dose"]] <- FALSE
      dose_all <- NULL
    } else {
      dose_all_df <- as.data.frame(dose_all, stringsAsFactors = FALSE)
      sig_idx <- dose_all_df[["p.adjust"]] <= pcutoff
      dose_sig_df <- dose_all_df[sig_idx, ]
      mesg("Found ", nrow(dose_sig_df), " DOSE over-representation hits.")
    }
  }

  gse_dose <- NULL
  gse_dose_all_df <- gse_dose_sig_df <- data.frame()
  if (isTRUE(todo[["gse_dose"]])) {
    mesg("Performing DOSE GSEA.")
    gse_dose <- try(
      DOSE::gseDO(geneList = genelist, ont = do_db, organism = orgn,
                  pvalueCutoff = 1.0, minGSSize = min_groupsize,
                  maxGSSize = max_groupsize))
    if ("try-error" %in% class(gse_dose)) {
      gse_dose <- NULL
    } else {
      gse_dose_all_df <- as.data.frame(gse_dose, stringsAsFactors = FALSE)
      sig_idx <- gse_dose_all_df[["p.adjust"]] <= pcutoff
      gse_dose_sig_df <- gse_dose_all_df[sig_idx, ]
      mesg("Found ", nrow(gse_dose_sig_df), " Dose GSE hits.")
    }
  }
  dose_data <- list(
      "dose_enrich" = dose_all,
      "dose_all" = dose_all_df,
      "dose_sig" = dose_sig_df,
      "dose_gse" = gse_dose,
      "dose_gse_all" = gse_dose_all_df,
      "dose_gse_sig" = gse_dose_sig_df)

  mesh_org <- "Homo sapiens"
  mesh_category <- "C"
  mesh_dbname <- "gendoo"
  mesh_all <- NULL
  mesh_all_df <- mesh_sig_df <- data.frame()
  if (isTRUE(todo[["enrich_mesh"]])) {
    if (organism == "human" || organism == "Homo sapiens") {
      loaded <- sm(requireNamespace(package = "MeSH.Hsa.eg.db", quietly = TRUE))
    } else if (organism == "mouse" || organism == "Mus musculus") {
      loaded <- sm(requireNamespace(package = "MeSH.Mm.eg.db", quietly = TRUE))
      mesh_org <- "Mus musculus"
    } else {
      todo[["enrich_mesh"]] <- FALSE
      todo[["gse_mesh"]] <- FALSE
    }
  }
  if (isTRUE(todo[["enrich_mesh"]])) {
    ah <- AnnotationHub::AnnotationHub()
    loaded <- sm(requireNamespace(package = "meshes", quietly = TRUE))
    if (is.null(ah_data)) {
      ah_data <- sm(AnnotationHub::query(ah, c("MeSHDb", mesh_org)))
    }
    orgn_db <- sm(ah_data[[1]])
    if (is.null(mesh_db)) {
      mesh_db <- sm(MeSHDbi::MeSHDb(orgn_db))
    }
    mesh_all <- try(meshes::enrichMeSH(
      gene = sig_gene_list, MeSHDb = mesh_db, database = mesh_dbname,
      pvalueCutoff = 1.0, pAdjustMethod = padj_type, universe = universe_to,
      minGSSize = min_groupsize, maxGSSize = max_groupsize,
      qvalueCutoff = qcutoff), silent = TRUE)
    if ("try-error" %in% class(mesh_all)) {
      todo[["gse_mesh"]] <- FALSE
      mesh_all <- NULL
    } else {
      mesh_all_df <- as.data.frame(mesh_all, stringsAsFactors = FALSE)
      sig_idx <- mesh_all_df[["p.adjust"]] <= pcutoff
      mesh_sig_df <- mesh_all_df[sig_idx, ]
      mesg("Found ", nrow(mesh_sig_df), " MESH over-representation hits.")
    }
  }

  gse_mesh <- NULL
  gse_mesh_all_df <- gse_mesh_sig_df <- data.frame()
  if (isTRUE(todo[["gse_mesh"]])) {
    gse_mesh <- try(
      meshes::gseMeSH(
        geneList = genelist, MeSHDb = mesh_db, database = mesh_dbname,
        pvalueCutoff = 1.0, minGSSize = min_groupsize,
        maxGSSize = max_groupsize))
    if ("try-error" %in% class(gse_mesh)) {
      gse_mesh <- NULL
    } else {
      gse_mesh_all_df <- as.data.frame(gse_mesh, stringsAsFactors = FALSE)
      sig_idx <- gse_mesh_all_df[["p.adjust"]] <= pcutoff
      gse_mesh_sig_df <- gse_mesh_all_df[sig_idx, ]
      mesg("Found ", nrow(gse_mesh_sig_df), " MESH GSE hits.")
    }
  }
  mesh_data <- list(
    "mesh_enrich" = mesh_all,
    "mesh_all" = mesh_all_df,
    "mesh_sig" = mesh_sig_df,
    "mesh_gse" = gse_mesh,
    "mesh_gse_all" = gse_mesh_all_df,
    "mesh_gse_sig" = gse_mesh_sig_df)

  msigdb_all <- NULL
  msigdb_all_df <- msigdb_sig_df <- data.frame()
  if (isTRUE(todo[["enrich_msigdb"]])) {
    ## Currently my msigdb converter only does gene symbols...
    if (is.null(signature_data)) {
      signature_data <- load_gmt_signatures(
        signatures = msig_db, signature_category = msigdb_category, id_type = "entrez")
    }
    if (is.null(signature_df)) {
      signature_df <- signatures_to_df(signature_data)
    }
    expected_term_type <- "ENTREZID"
    current_term_type <- "ENTREZID"
    if (current_term_type != expected_term_type) {
      test_sig_df <- sm(try(clusterProfiler::bitr(
        sig_gene_list, fromType = current_term_type, toType = expected_term_type, OrgDb = org),
        silent = TRUE))
    }
    mesg("Starting msigDB over-representation analysis.")
    msigdb_all <- try(clusterProfiler::enricher(
      sig_gene_list, TERM2GENE = signature_df))
    if ("try-error" %in% class(msigdb_all)) {
      todo[["gse_msigdb"]] <- FALSE
      msigdb_all <- NULL
    } else {
      msigdb_all_df <- as.data.frame(msigdb_all, stringsAsFactors = FALSE)
      sig_idx <- msigdb_all_df[["p.adjust"]] <= pcutoff
      msigdb_sig_df <- msigdb_all_df[sig_idx, ]
      mesg("Found ", nrow(msigdb_sig_df), " KEGG GSE hits.")
    }
  }

  gse_msigdb_all <- NULL
  gse_msigdb_all_df <- gse_msigdb_sig_df <- data.frame()
  if (isTRUE(todo[["gse_msigdb"]])) {
    mesg("Starting msigDB GSEA.")
    gse_msigdb_all <- try(
      clusterProfiler::GSEA(
        genelist, exponent = 1, minGSSize = min_groupsize, maxGSSize = max_groupsize,
        eps = 1e-10, pvalueCutoff = 1.0, pAdjustMethod = padj_type,
        TERM2GENE = signature_df))
    if ("try-error" %in% class(gse_msigdb_all)) {
      gse_msigdb_all <- NULL
    } else {
      gse_msigdb_all_df <- as.data.frame(gse_msigdb_all, stringsAsFactors = FALSE)
      sig_idx <- gse_msigdb_all_df[["p.adjust"]] <= pcutoff
      gse_msigdb_sig_df <- gse_msigdb_all_df[sig_idx, ]
      mesg("Found ", nrow(gse_msigdb_sig_df), " msigDB GSE hits.")
    }
  }
  msigdb_data <- list(
    "msigdb_enrich" = msigdb_all,
    "msigdb_all" = msigdb_all_df,
    "msigdb_sig" = msigdb_sig_df,
    "msigdb_gse" = gse_msigdb_all,
    "msigdb_gse_all" = gse_msigdb_all_df,
    "msigdb_gse_sig" = gse_msigdb_sig_df)

  retlist <- list(
    "all_mappings" = de_table_namedf,
    "go_data" = go_data,
    "kegg_data" = kegg_data,
    "david_data" = david_data,
    "reactome_data" = reactome_data,
    "dose_data" = dose_data,
    "mesh_data" = mesh_data,
    "msigdb_data" = msigdb_data,
    "todo" = todo,
    "orgdb_from" = orgdb_from,
    "orgdb_to" = orgdb_to,
    "kegg_organism" = kegg_organism,
    "kegg_universe" = kegg_universe,
    "reactome_organism" = reactome_organism,
    "mesh_db" = mesh_db,
    "ah_data" = ah_data,
    "signature_data" = signature_data,
    "signature_df" = signature_df,
    "de_table_namedf" = de_table_namedf,
    "sig_genes_namedf" = sig_genes_namedf)
  if (!is.null(excel)) {
    mesg("Writing data to: ", excel, ".")
    excel_ret <- try(write_cp_data(retlist, excel = excel))
    if (class(excel_ret) == "try-error") {
      message("Writing the data to excel failed.")
    } else {
      retlist[["excel"]] <- excel_ret
    }
  }
  class(retlist) <- c("clusterprofiler_result", "list")
  return(retlist)
}

#' Print a clusterprofiler over representation search.
#'
#' @param x Monstrous list of the various results, including but not
#'  limited to plots, go-gene mappings, enrichmed, kegg, david, GO
#'  analyses.
#' @param ... Other args to match the generic.
#' @export
print.clusterprofiler_result <- function(x, ...) {
  mesg("A set of ontologies produced by clusterprofiler.")
  return(invisible(x))
}

#' Set up appropriate option sets for clusterProfiler
#'
#' This hard-sets some defaults for orgdb/kegg databases when using
#' clusterProfiler.
#'
#' @param species  Currently it only works for humans and fruit flies.
cp_options <- function(species) {
  if (species == "dmelanogaster") {
    options <- list(
      orgdb = "org.Dm.eg.db",
      orgdb_from = "FLYBASE",
      orgdb_to = c("ENSEMBL", "SYMBOL", "ENTREZID"),
      kegg_prefix = "Dmel_",
      kegg_organism = "dme",
      kegg_id_column = "FLYBASECG")
  } else if (species == "hsapiens") {
    options <- list(
      orgdb = "org.Hs.eg.db",
      orgdb_from = "ENSEMBL",
      orgdb_to = c("ENSEMBL", "SYMBOL", "ENTREZID"),
      kegg_prefix = "Hsa_",
      kegg_organism = "hsa",
      kegg_id_column = "")
  }
  return(options)
}

#' Generic enrichment using clusterProfiler.
#'
#' culsterProfiler::enricher provides a quick and easy enrichment analysis given
#' a set of siginficant' genes and a data frame which connects each gene to a
#' category.
#'
#' @param sig_genes Set of 'significant' genes as a table.
#' @param de_table All genes from the original analysis.
#' @param db Dataframe of GO->ID matching the gene names of sig_genes to GO
#'  categories.
#' @param current_id Starting ID of the genes.
#' @param needed_id ID type to coerce the data to.
#' @param org DBI to do the conversion.
#' @return Table of 'enriched' categories.
#' @export
simple_cp_enricher <- function(sig_genes, de_table, db, current_id = "ENSEMBL",
                               needed_id = "SYMBOL", org = "org.Hs.eg.db") {
  all_genenames <- rownames(de_table)
  sig_genenames <- rownames(sig_genes)
  if (current_id != needed_id) {
    test_genes <- clusterProfiler::bitr(rownames(sig_genes), fromType = current_id,
                                        toType = needed_id, OrgDb = org)
    test_genes <- test_genes[[needed_id]]
  } else {
    test_genes <- rownames(sig_genes)
  }
  enriched <- clusterProfiler::enricher(test_genes, TERM2GENE = db)
  return(enriched)
}
setGeneric("simple_cp_enricher")

#' Invoke simple_cp_enricher when the input database is a GeneSetCollection.
#'
#' @param sig_genes Set of 'significant' genes as a table.
#' @param de_table All genes from the original analysis.
#' @param db Dataframe of GO->ID matching the gene names of sig_genes to GO
#'  categories.
#' @param current_id Starting ID of the genes.
#' @param needed_id ID type to coerce the data to.
#' @param org DBI to do the conversion.
#' @return Table of 'enriched' categories.
#' @export
setMethod(
  "simple_cp_enricher", signature = signature(db = "GeneSetCollection"),
  definition = function(sig_genes, de_table, db, current_id = "ENSEMBL",
                        needed_id = "SYMBOL", org = "org.Hs.eg.db") {
    db <- signatures_to_df(db)
    simple_cp_enricher(sig_genes, de_table, db = db, current_id = current_id,
                       needed_id = needed_id, org = org)
  })

## EOF
