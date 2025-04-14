## ontology_clusterprofiler.r: Methods to simplify using clusterProfiler.  I
## think clusterprofiler is probably the most complete and comprehensive GSEA
## toolkit.  It is not necessarily as easy to use as I might prefer.  This seeks
## to fill in some corner case.

#' Run simple_clusterprofiler on every table from extract_significant_genes()
#'
#' @param sig Result from extract_significant_genes
#' @param tables Result from combine_de_tables
#' @param according_to Use this result type for the clusterprofiler searches.
#' @param together Concatenate the up/down genes into one set?
#' @param plot_type Choose a plot method as the default.
#' @param ... Arguments to pass to simple_clusterprofiler().
#' @export
all_cprofiler <- function(sig, tables, according_to = "deseq", together = FALSE,
                          orgdb = "org.Hs.eg.db", orgdb_from = NULL, orgdb_to = "ENTREZID",
                          go_level = 3, pcutoff = 0.05, qcutoff = 0.2, fc_column = "logFC",
                          second_fc_column = "deseq_logfc", internal = FALSE, updown = "up",
                          permutations = 1000, min_groupsize = 5, kegg_prefix = NULL,
                          kegg_organism = NULL, do_gsea = TRUE, categories = 12, do_david = FALSE,
                          do_kegg = TRUE, padj_type = "BH", plot_type = "all",
                          do_reactome = TRUE,
                          excel = "excel/all_cp.xlsx", ...) {
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
      ret[[retname_up]] <- simple_clusterprofiler(
        up, table,
        orgdb = orgdb, orgdb_from = orgdb_from, orgdb_to = orgdb_to, go_level = go_level,
        pcutoff = pcutoff, qcutoff = qcutoff, fc_column = fc_column,
        second_fc_column = second_fc_column, internal = internal, updown = updown,
        permutations = permutations, min_groupsize = min_groupsize, kegg_prefix = kegg_prefix,
        kegg_organism = kegg_organism, do_gsea = do_gsea, categories = categories,
        do_david = do_david, do_kegg = do_kegg, padj_type = padj_type,
        do_reactome = do_reactome, excel = chosen_up_xlsx,
        ...)
    } else {
      ret[[retname_up]] <- NULL
    }
    if (down_elements > 0) {
      slept <- Sys.sleep(10)
      chosen_down_xlsx <- file.path(xlsx_dir, glue("{xlsx_base}_{retname_down}.xlsx"))
      ret[[retname_down]] <- simple_clusterprofiler(
        down, table,
        orgdb = orgdb, orgdb_from = orgdb_from, orgdb_to = orgdb_to, go_level = go_level,
        pcutoff = pcutoff, qcutoff = qcutoff, fc_column = fc_column,
        second_fc_column = second_fc_column, internal = internal, updown = updown,
        permutations = permutations, min_groupsize = min_groupsize, kegg_prefix = kegg_prefix,
        kegg_organism = kegg_organism, do_gsea = do_gsea, categories = categories,
        do_david = do_david, do_kegg = do_kegg, padj_type = padj_type,
        do_reactome = do_reactome, excel = chosen_down_xlsx,
        ...)
      #ret[[retname_down]] <- sm(simple_clusterprofiler(down, first_col = fc_col))
    } else {
      ret[[retname_down]] <- NULL
    }
  }
  class(ret) <- "all_cprofiler"
  return(ret)
}

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
  orgdb_sig_from <- orgdb_from
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
#' @return a list
#' @seealso [clusterProfiler] [AnnotationDbi] [KEGGREST]
#' @examples
#' \dontrun{
#'  holyasscrackers <- simple_clusterprofiler(gene_list, all_genes, "org.Dm.eg.db")
#' }
#' @export
simple_clusterprofiler <- function(sig_genes, de_table = NULL, orgdb = "org.Hs.eg.db",
                                   orgdb_from = NULL, orgdb_to = "ENTREZID",
                                   go_level = 3, pcutoff = 0.05, organism = "human",
                                   qcutoff = 0.2, fc_column = "logFC",
                                   second_fc_column = "deseq_logfc", internal = FALSE,
                                   updown = "up", permutations = 1000, min_groupsize = 5,
                                   max_groupsize = 500, kegg_prefix = NULL, kegg_organism = NULL,
                                   do_gsea = TRUE, categories = 12, excel = NULL, do_david = FALSE,
                                   do_kegg = TRUE, david_id = "ENTREZ_GENE_ID", padj_type = "BH",
                                   david_user = "abelew@umd.edu", do_reactome = TRUE,
                                   do_dose = FALSE, do_mesh = FALSE, do_msigdb = FALSE,
                                   mesh_category = "C", mesh_dbname = "gendoo",
                                   msigdb_category = "C2", msig_db = NULL) {
  loaded <- sm(requireNamespace(package = "clusterProfiler", quietly = TRUE))
  loaded <- sm(requireNamespace(package = "DOSE", quietly = TRUE))
  org <- NULL
  ## Start off by figuring out what was given, an OrgDb or the name of one.
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
  }
  sig_genenames <- rownames(sig_genes)
  orgdb_to <- toupper(orgdb_to)

  de_table_namedf <- NULL
  sig_genes_namedf <- NULL
  test_genes_df <- NULL
  test_sig_df <- NULL
  num_sig <- 0
  num_hits <- 0
  orgdb_sig_from <- orgdb_from
  if (is.null(orgdb_from)) {
    mesg("Testing available OrgDb keytypes for the best mapping to entrez.")
    for (k in mapper_keys) {
      test_genes_df <- sm(try(clusterProfiler::bitr(all_genenames, fromType = k,
                                                    toType = orgdb_to, OrgDb = org), silent = TRUE))
      test_sig_df <- sm(try(clusterProfiler::bitr(sig_genenames, fromType = k,
                                                  toType = orgdb_to, OrgDb = org), silent = TRUE))
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

  } else { ## If we do have a column for the OrgDB
    de_table_namedf <- sm(try(clusterProfiler::bitr(all_genenames, fromType = orgdb_from,
                                                    toType = orgdb_to, OrgDb = org), silent = TRUE))
    sig_genes_namedf <- sm(try(clusterProfiler::bitr(sig_genenames, fromType = orgdb_from,
                                                     toType = orgdb_to, OrgDb = org), silent = TRUE))
  }

  if (is.null(sig_genes[[fc_column]]) && is.null(sig_genes[[second_fc_column]])) {
    stop("The fold change column provided no genes, try another column in the data set.")
  } else if (is.null(sig_genes[[fc_column]])) {
    fc_column <- second_fc_column
  }

  gsea_fc_column <- fc_column
  if (is.null(de_table[[gsea_fc_column]]) && is.null(de_table[[second_fc_column]])) {
    message("Unable to find the fold-change column in the de table, not doing gsea.")
    do_gsea <- FALSE
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
  ggo_mf <- ggo_bp <- ggo_cc <- NULL
  ggo_mf <- sm(clusterProfiler::groupGO(gene = sig_gene_list, OrgDb = org,
                                        keyType = orgdb_to,
                                        ont = "MF", level = go_level))
  ggo_bp <- sm(clusterProfiler::groupGO(gene = sig_gene_list, OrgDb = org,
                                        keyType = orgdb_to,
                                        ont = "BP", level = go_level))
  ggo_cc <- sm(clusterProfiler::groupGO(gene = sig_gene_list, OrgDb = org,
                                        keyType = orgdb_to,
                                        ont = "CC", level = go_level))
  ## Recast those groups as dataframes to send back to the user.
  group_go <- list(
    "MF" = as.data.frame(ggo_mf, stringsAsFactors = FALSE),
    "BP" = as.data.frame(ggo_bp, stringsAsFactors = FALSE),
    "CC" = as.data.frame(ggo_cc, stringsAsFactors = FALSE))
  mesg("Found ", nrow(group_go[["MF"]]),
       " MF, ", nrow(group_go[["BP"]]),
       " BP, and ", nrow(group_go[["CC"]]), " CC hits.")
  mesg("Calculating enriched GO groups.")
  ego_all_mf <- clusterProfiler::enrichGO(gene = sig_gene_list, universe = universe_to,
                                          OrgDb = org, ont = "MF", keyType = orgdb_to,
                                          minGSSize = min_groupsize, pAdjustMethod = padj_type,
                                          pvalueCutoff = 1.0)
  ego_sig_mf <- clusterProfiler::enrichGO(gene = sig_gene_list, universe = universe_to,
                                          OrgDb = org, ont = "MF", keyType = orgdb_to,
                                          minGSSize = min_groupsize, pAdjustMethod = padj_type,
                                          pvalueCutoff = pcutoff)
  ego_all_bp <- clusterProfiler::enrichGO(gene = sig_gene_list, universe = universe_to,
                                          OrgDb = org, ont = "BP", keyType = orgdb_to,
                                          minGSSize = min_groupsize, pAdjustMethod = padj_type,
                                          pvalueCutoff = 1.0)
  ## Now extract the p-value significant categories.
  ## This is at least a little bit redundant and should perhaps be revisited/removed.
  ego_sig_bp <- clusterProfiler::enrichGO(gene = sig_gene_list, universe = universe_to,
                                          OrgDb = org, ont = "BP", keyType = orgdb_to,
                                          minGSSize = min_groupsize, pAdjustMethod = padj_type,
                                          pvalueCutoff = pcutoff)
  ego_all_cc <- clusterProfiler::enrichGO(gene = sig_gene_list, universe = universe_to,
                                          OrgDb = org, ont = "CC", keyType = orgdb_to,
                                          minGSSize = min_groupsize, pAdjustMethod = padj_type,
                                          pvalueCutoff = 1.0)
  ego_sig_cc <- clusterProfiler::enrichGO(gene = sig_gene_list, universe = universe_to,
                                          OrgDb = org, ont = "CC", keyType = orgdb_to,
                                          minGSSize = min_groupsize, pAdjustMethod = padj_type,
                                          pvalueCutoff = pcutoff)
  ## Once again, recast the results for the user.
  enrich_go <- list(
    "MF_all" = as.data.frame(ego_all_mf, stringsAsFactors = FALSE),
    "MF_sig" = as.data.frame(ego_sig_mf, stringsAsFactors = FALSE),
    "BP_all" = as.data.frame(ego_all_bp, stringsAsFactors = FALSE),
    "BP_sig" = as.data.frame(ego_sig_bp, stringsAsFactors = FALSE),
    "CC_all" = as.data.frame(ego_all_cc, stringsAsFactors = FALSE),
    "CC_sig" = as.data.frame(ego_sig_cc, stringsAsFactors = FALSE))
  mesg("Found ", nrow(enrich_go[["MF_sig"]]),
       " MF, ", nrow(enrich_go[["BP_sig"]]),
       " BP, and ", nrow(enrich_go[["CC_sig"]]), " CC enriched hits.")
  ## Set up to do GSEA
  gse_go <- list()
  de_table_merged <- NULL
  gse <- list()
  if (isTRUE(do_gsea)) {
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
    ## Hmm this is odd, in the previous calls, I used orgdb_to, but in this set
    ## I am using orgdb_from...
    mesg("Performing GSE analyses of gene lists (this is slow).")
    genelist <- as.vector(de_table_merged[[gsea_fc_column]])
    names(genelist) <- de_table_merged[[orgdb_to]]
    duplicated_names <- duplicated(names(genelist))
    num_duplicates <- sum(duplicated_names)
    if (num_duplicates > 0) {
      mesg("There are ", num_duplicates, " duplicated gene IDs, dropping them.")
      genelist <- genelist[!duplicated_names]
    }
    ## 2020 04: Adding a pvalue cutoff argument causes an error, I do not know why.
    ## Arguments used by gseGO of interest: exponent, minGSSize/maxGSSize, eps, by(fgsea)
    ## Also, apparently the nperm argument is deprecated.
    gse <- suppressWarnings(clusterProfiler::gseGO(geneList = genelist, OrgDb = org,
                                                   keyType = orgdb_to, ont = "ALL",
                                                   minGSSize = min_groupsize))
    gse_go <- as.data.frame(gse)
    mesg("Found ", nrow(gse_go), " enriched hits.")
  } else {
    genelist <- as.vector(sig_genes[[fc_column]])
    names(genelist) <- rownames(sig_genes)
  }
  ## Set up to do kegg
  ## Now extract the kegg organism/gene IDs.
  if (is.null(kegg_organism)) {
    org_meta <- AnnotationDbi::metadata(org)
    org_row <- org_meta[["name"]] == "ORGANISM"
    organism <- org_meta[org_row, "value"]
    ## Only grab the first of potentially multiple outputs.
    kegg_organism <- get_kegg_orgn(species = organism)
    do_kegg <- TRUE
    if (length(kegg_organism) > 0) {
      kegg_organism <- kegg_organism[[1]]
    } else {
      do_kegg <- FALSE
      kegg_organism <- NULL
    }
  }
  all_kegg <- enrich_kegg <- NULL
  if (isTRUE(do_kegg)) {
    kegg_universe <- KEGGREST::keggConv(kegg_organism, "ncbi-geneid")
    kegg_sig_names <- glue("ncbi-geneid:{sig_gene_list}")
    kegg_sig_intersect <- kegg_sig_names %in% names(kegg_universe)
    mesg("Found ", sum(kegg_sig_intersect),
         " matches between the significant gene list and kegg universe.")

    if (sum(kegg_sig_intersect) > 0) {
      all_names <- names(kegg_universe)
      small_universe <- kegg_universe[intersect(kegg_sig_names, all_names)]
      kegg_sig_ids <- unique(as.character(small_universe))
      ##kegg_sig_ids <- unique(as.character(kegg_universe[kegg_sig_intersect]))
      kegg_sig_ids <- gsub(pattern = glue("{kegg_organism}:"),
                           replacement = "", x = kegg_sig_ids)
      mesg("Performing KEGG analyses.")
      all_kegg <- clusterProfiler::enrichKEGG(kegg_sig_ids, organism = kegg_organism,
                                              keyType = "kegg",
                                              pvalueCutoff = 1.0)
      enrich_kegg <- sm(clusterProfiler::enrichKEGG(kegg_sig_ids, organism = kegg_organism,
                                                    keyType = "kegg",
                                                    pvalueCutoff = pcutoff))
    }
  }

  if (is.null(all_kegg)) {
    do_gsea <- FALSE
  }

  gse_all_kegg <- NULL
  gse_sig_kegg <- NULL
  if (isTRUE(do_gsea)) {
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

    gse_all_kegg <- sm(
      clusterProfiler::gseKEGG(geneList = kegg_genelist, organism = kegg_organism,
                               nPerm = permutations, minGSSize = min_groupsize,
                               pvalueCutoff = 1.0, use_internal_data = internal))
    gse_sig_kegg <- sm(
      clusterProfiler::gseKEGG(geneList = kegg_genelist, organism = kegg_organism,
                               nPerm = permutations, minGSSize = min_groupsize,
                               pvalueCutoff = pcutoff, use_internal_data = internal))
  }

  kegg_data <- list(
    "kegg_all" = as.data.frame(all_kegg, stringsAsFactors = FALSE),
    "kegg_sig" = as.data.frame(enrich_kegg, stringsAsFactors = FALSE),
    "kegg_gse_all" = as.data.frame(gse_all_kegg, stringsAsFactors = FALSE),
    "kegg_gse_sig" = as.data.frame(gse_sig_kegg, stringsAsFactors = FALSE))
  mesg("Found ", nrow(kegg_data[["kegg_sig"]]), " KEGG enriched hits.")

  david_data <- NULL
  if (isTRUE(do_david)) {
    mesg("Attempting DAVID search.")
    david_search <- try(clusterProfiler::enrichDAVID(
      gene = sig_gene_list, minGSSize = min_groupsize,
      idType = david_id, david.user = david_user), silent = TRUE)
    if (class(david_search)[[1]] == "try-error") {
      david_data <- NULL
    } else {
      david_data <- as.data.frame(david_search, stringsAsFactors = FALSE)
    }
    mesg("Found ", nrow(david_data), " DAVID hits.")
  }

  reactome_data <- NULL
  reactome_organism <- "human"
  if (isTRUE(do_reactome)) {
    library(ReactomePA)
    if (orgdb == "org.Mm.eg.db") {
      reactome_organism <- "mouse"
    }
    reactome_data <- ReactomePA::enrichPathway(
      gene = sig_gene_list, pvalueCutoff = pcutoff, readable = TRUE,
      pAdjustMethod = padj_type, qvalueCutoff = qcutoff, universe = universe_to,
      organism = reactome_organism,
      minGSSize = min_groupsize, maxGSSize = max_groupsize)
  }

  dose_data <- NULL
  orgn <- "hsa"
  do_db <- "HDO"
  if (isTRUE(do_dose)) {
    if (organism == "human") {
      loaded <- sm(requireNamespace(package = "HDO.db", quietly = TRUE))
    } else if (organism == "mouse") {
      orgn <- "mm"
      do_db <- "MPO"
      loaded <- sm(requireNamespace(package = "MPO.db", quietly = TRUE))
    } else {
      warning("I do not know this DOSE organism, leaving it as human.")
    }
    dose_data <- DOSE::enrichDO(
      gene = sig_gene_list, ont = do_db, organism = orgn,
      pvalueCutoff = pcutoff, pAdjustMethod = padj_type, universe = universe_to,
      minGSSize = min_groupsize, maxGSSize = max_groupsize,
      qvalueCutoff = qcutoff, readable = TRUE)
  }

  mesh_data <- NULL
  mesh_org <- "Homo sapiens"
  mesh_category <- "C"
  mesh_dbname <- "gendoo"
  if (isTRUE(do_mesh)) {
    ah <- AnnotationHub::AnnotationHub()
    loaded <- sm(requireNamespace(package = "meshes", quietly = TRUE))
    if (organism == "human") {
      loaded <- sm(requireNamespace(package = "MeSH.Hsa.eg.db", quietly = TRUE))
    } else if (organism == "mouse") {
      loaded <- sm(requireNamespace(package = "MeSH.Mm.eg.db", quietly = TRUE))
      mesh_org <- "Mus musculus"
    } else {
      warning("I do not know this mesh organism, leaving it as human.")
    }
    ah_data <- AnnotationHub::query(ah, c("MeSHDb", mesh_org))
    orgn_db <- ah_data[[1]]
    mesh_db <- MeSHDbi::MeSHDb(orgn_db)

    mesh_data <- try(meshes::enrichMeSH(
      gene = sig_gene_list, MeSHDb = mesh_db, database = mesh_dbname,
      pvalueCutoff = pcutoff, pAdjustMethod = padj_type, universe = universe_to,
      minGSSize = min_groupsize, maxGSSize = max_groupsize,
      qvalueCutoff = qcutoff))
    if ("try-error" %in% class(mesh_data)) {
      mesh_data <- data.frame()
    }
  }

  msigdb_data <- NULL
  if (isTRUE(do_msigdb)) {
    ## Currently my msigdb converter only does gene symbols...
    signature_data <- load_gmt_signatures(
      signatures = msig_db, signature_category = msigdb_category, id_type = "entrez")
    signature_df <- signatures_to_df(signature_df)
    test_sig_df <- sm(try(clusterProfiler::bitr(sig_gene_list, fromType = "ENTREZID",
                                                toType = "SYMBOL", OrgDb = org), silent = TRUE))
    msigdb_data <- clusterProfiler::enricher(test_sig_df[["SYMBOL"]], TERM2GENE = signature_df)
  }

  mesg("Plotting results, removing most of this.")
  net_sig_mf <- try(
    clusterProfiler::cnetplot(ego_sig_mf, categorySize = "pvalue",
                              color.params = list(foldChange = genelist)), silent = TRUE)
  net_sig_bp <- try(
    clusterProfiler::cnetplot(ego_sig_bp, categorySize = "pvalue",
                              color.params = list(foldChange = genelist)), silent = TRUE)
  net_sig_cc <- try(
    clusterProfiler::cnetplot(ego_sig_cc, categorySize = "pvalue",
                              color.params = list(foldChange = genelist)), silent = TRUE)
  plotlist <- list(
    "net_sig_mf" = net_sig_mf,
    "net_sig_bp" = net_sig_bp,
    "net_sig_cc" = net_sig_cc)
  enrich_objects <- list(
    "MF_all" = ego_all_mf,
    "MF_sig" = ego_sig_mf,
    "BP_all" = ego_all_bp,
    "BP_sig" = ego_sig_bp,
    "CC_all" = ego_all_cc,
    "CC_sig" = ego_sig_cc,
    "gse" = gse,
    "all_kegg" = all_kegg,
    "enrich_kegg" = enrich_kegg,
    "gse_all_kegg" = gse_all_kegg,
    "gse_enrich_kegg" = gse_sig_kegg)
  retlist <- list(
    "all_mappings" = de_table_namedf,
    "sig_mappings" = sig_genes_namedf,
    "group_go" = group_go,
    "enrich_go" = enrich_go,
    "enrich_objects" = enrich_objects,
    "gse_go" = gse_go,
    "kegg_data" = kegg_data,
    "reactome_data" = reactome_data,
    "dose_data" = dose_data,
    "mesh_data" = mesh_data,
    "msigdb_data" = msigdb_data,
    "david_data" = david_data,
    "plots" = plotlist,
    "pvalue_plots" = plotlist)
  class(retlist) <- c("clusterprofiler_result", "list")
  if (!is.null(excel)) {
    mesg("Writing data to: ", excel, ".")
    excel_ret <- try(write_cp_data(retlist, excel = excel))
    if (class(excel_ret) == "try-error") {
      message("Writing the data to excel failed.")
    }
  }
  return(retlist)
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
#' @param go_db Dataframe of GO->ID matching the gene names of sig_genes to GO
#'  categories.
#' @return Table of 'enriched' categories.
simple_cp_enricher <- function(sig_genes, de_table, go_db = NULL) {
  all_genenames <- rownames(de_table)
  sig_genenames <- rownames(sig_genes)
  enriched <- clusterProfiler::enricher(sig_genenames, TERM2GENE = go_db)
  retlist <- list(
    "enriched" = as.data.frame(enriched, stringsAsFactors = FALSE))
  return(retlist)
}

## EOF
