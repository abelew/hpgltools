## ontology_clusterprofiler.r: Methods to simplify using clusterProfiler.  I
## think clusterprofiler is probably the most complete and comprehensive GSEA
## toolkit.  It is not necessarily as easy to use as I might prefer.  This seeks
## to fill in some corner case.

#' @include 01_hpgltools.R
NULL

#' Run simple_clusterprofiler on every table from extract_significant_genes()
#'
#' Most of the options are taken from simple_cprofiler and passed directly down.
#'
#' @param sig Result from extract_significant_genes.
#' @param tables Result from combine_de_tables.
#' @param according_to Use this result type for the clusterprofiler searches.
#' @param together Concatenate the up/down genes into one set?
#' @param todo List of methods to try enrichment/gsea
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
#' @param permutations Run x permutations in clusterProfiler. (deprecated by cp authors)
#' @param min_groupsize Ignore groups with less than x genes.
#' @param max_groupsize Ignore groups with more than x genes.
#' @param kegg_prefix Prefix of this kegg organism ID.
#' @param kegg_organism Full name of this organism when querying KEGG.
#' @param organism Organism for things like DOSE etc, TODO: double check that this is still needed
#' @param categories Plot this number of categories by default.
#' @param padj_type FDR correction method.
#' @param mesh_category Mesh category to query
#' @param mesh_dbname Database in mesh to query
#' @param msigdb_category Category of mSigDB to query
#' @param msig_db Dataframe of mSigDB gene sets.
#' @param kegg_universe Manually set universe of kegg genes/IDs.
#' @param reactome_organism I think not needed any longer.
#' @param mesh_db Location of mesh database/package.
#' @param ah_data Which annotationhub dataset to use with mesh
#' @param signature_data Manually provided list of gene sets.
#' @param signature_df Manually provided df of signatures.
#' @param de_table_namedf Alternate set of gene names.
#' @param sig_genes_namedf Alternate set of gene names for the significant genes.
#' @param excel Output xlsx filename.
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
  ret <- list()
  input_up <- list()
  input_down <- list()
  source <- "significant"

  xlsx_dir <- dirname(excel)
  xlsx_base <- gsub(x = basename(excel), pattern = "\\.[[:alpha:]]{3,}$", replacement = "")
  ## Check if this came from extract_significant_genes or extract_abundant_genes.
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
      Sys.sleep(10)
      chosen_down_xlsx <- file.path(xlsx_dir, glue("{xlsx_base}_{retname_down}.xlsx"))
      args <- list(
        "sig_genes" = down, "de_table" = table, "orgdb" = orgdb, "todo" = todo,
        "orgdb_from" = orgdb_from, "orgdb_to" = orgdb_to, "go_level" = go_level,
        "pcutoff" = pcutoff, "qcutoff" = qcutoff, "fc_column" = fc_column,
        "second_fc_column" = second_fc_column, "internal" = internal, "updown" = updown,
        "permutations" = permutations, "min_groupsize" = min_groupsize, "kegg_prefix" = kegg_prefix,
        "kegg_organism" = kegg_organism, "categories" = categories, "padj_type" = padj_type,
        "excel" = chosen_down_xlsx, "organism" = organism,
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
  class(ret) <- "hpgltools::all_cprofiler"
  return(ret)
}

cp_go_enrich <- function(sig_gene_list, org, orgdb_to, level, universe_to, min_groupsize,
                         padj_type, pcutoff = 0.05) {
  ego_all_mf <- ego_all_bp <- ego_all_cc <- NULL ## GO enrichment data structures
  ego_all_mf_df <- ego_all_bp_df <- ego_all_cc_df <- NULL ## GO enrichment full dfs
  ego_sig_mf_df <- ego_sig_bp_df <- ego_sig_cc_df <- NULL ## pvalue cutoff dfs
  ggo_mf <- try(clusterProfiler::groupGO(
    gene = sig_gene_list, Org = org, keyType = orgdb_to, ont = "MF",
    level = level), silent = TRUE)
  ggo_bp <- try(clusterProfiler::groupGO(
    gene = sig_gene_list, Org = org, keyType = orgdb_to, ont = "BP",
    level = level), silent = TRUE)
  ggo_cc <- try(clusterProfiler::groupGO(
    gene = sig_gene_list, Org = org, keyType = orgdb_to, ont = "CC",
    level = level), silent = TRUE)
  mesg("Performing GO over-representation analyses.")
  ego_all_mf <- clusterProfiler::enrichGO(
    gene = sig_gene_list, universe = universe_to, Org = org, ont = "MF", keyType = orgdb_to,
      minGSSize = min_groupsize, pAdjustMethod = padj_type, pvalueCutoff = 1.0)
    ego_all_mf_df <- as.data.frame(ego_all_mf, stringsAsFactors = FALSE)
    sig_idx <- ego_all_mf_df[["p.adjust"]] <= pcutoff
    ego_sig_mf_df <- ego_all_mf_df[sig_idx, ]
    ego_all_bp <- clusterProfiler::enrichGO(
      gene = sig_gene_list, universe = universe_to, Org = org, ont = "BP", keyType = orgdb_to,
      minGSSize = min_groupsize, pAdjustMethod = padj_type, pvalueCutoff = 1.0)
    ego_all_bp_df <- as.data.frame(ego_all_bp, stringsAsFactors = FALSE)
    sig_idx <- ego_all_bp_df[["p.adjust"]] <= pcutoff
    ego_sig_bp_df <- ego_all_bp_df[sig_idx, ]
    ego_all_cc <- clusterProfiler::enrichGO(
      gene = sig_gene_list, universe = universe_to, Org = org, ont = "CC", keyType = orgdb_to,
      minGSSize = min_groupsize, pAdjustMethod = padj_type, pvalueCutoff = 1.0)
    ego_all_cc_df <- as.data.frame(ego_all_cc, stringsAsFactors = FALSE)
    sig_idx <- ego_all_cc_df[["p.adjust"]] <= pcutoff
    ego_sig_cc_df <- ego_all_cc_df[sig_idx, ]
    mesg("Found ", nrow(ego_sig_mf_df),
         " MF, ", nrow(ego_sig_bp_df),
         " BP, and ", nrow(ego_sig_cc_df), " CC over-represented hits.")
  retlist <- list(
    "ggo_mf" = ggo_mf,
    "ggo_bp" = ggo_bp,
    "ggo_cc" = ggo_cc,
    "ego_all_mf" = ego_all_mf,
    "ego_all_bp" = ego_all_bp,
    "ego_all_cc" = ego_all_cc,
    "ego_all_mf_df" = ego_all_mf_df,
    "ego_all_bp_df" = ego_all_bp_df,
    "ego_all_cc_df" = ego_all_cc_df,
    "ego_sig_mf_df" = ego_sig_mf_df,
    "ego_sig_bp_df" = ego_sig_bp_df,
    "ego_sig_cc_df" = ego_sig_cc_df)
  class(retlist) <- "hpgltools::cp_go_gse"
  return(retlist)
}

cp_go_gsea <- function(genelist, org, orgdb_to, min_groupsize = 5, pcutoff = 0.05) {
  mesg("Performing GO GSEA of gene lists.")
  gse_go <- list()
  gse_go_all_df <- gse_go_sig_df <- data.frame()
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
  retlist <- list(
    "gse_go" = gse_go,
    "gse_go_all_df" = gse_go_all_df,
    "gse_go_sig_df" = gse_go_sig_df)
  class(retlist) <- "hpgltools::cp_go_gsea"
  return(retlist)
}

cp_kegg_organism <- function(organism, orgdb, kegg_organism = NULL) {
  if ("character" == class(orgdb)) {
    orgdb <- get0(orgdb)
  }
  if (is.null(kegg_organism)) {
    org_meta <- AnnotationDbi::metadata(orgdb)
    org_row <- org_meta[["name"]] == "ORGANISM"
    organism <- org_meta[org_row, "value"]
    ## Only grab the first of potentially multiple outputs.
    kegg_organism <- get_kegg_orgn(species = organism)
  }

  if (length(kegg_organism) == 1) {
    mesg("Found organism ID: ", kegg_organism, ".")
  } else if (length(kegg_organism) > 0) {
    kegg_organism <- kegg_organism[[1]]
    mesg("Using kegg organism ID: ", kegg_organism, ".")
  } else {
    kegg_organism <- NULL
  }
  return(kegg_organism)
}

cp_kegg_universe <- function(kegg_organism, sig_gene_list, kegg_universe = NULL) {
  if (is.null(kegg_universe)) {
    kegg_universe <- try(KEGGREST::keggConv(kegg_organism, "ncbi-geneid"))
  }
  kegg_sig_intersect <- 0
  kegg_sig_names <- NULL
  if ("try-error" %in% class(kegg_universe)) {
    return(NULL)
  } else {
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
  } else {
    return(NULL)
  }
  retlist <- list(
    "kegg_universe" = kegg_universe,
    "kegg_sig_ids" = kegg_sig_ids)
  return(retlist)
}

cp_kegg_enrich <- function(kegg_sig_ids, kegg_organism = NULL, pcutoff = 0.05) {
  mesg("Performing KEGG over-representation analysis.")
  all_kegg <- try(clusterProfiler::enrichKEGG(
    kegg_sig_ids, organism = kegg_organism, keyType = "kegg", pvalueCutoff = 1.0))
  if ("try-error" %in% class(all_kegg)) {
    return(NULL)
  } else {
    all_kegg_df <- as.data.frame(all_kegg, stringsAsFactors = FALSE)
    sig_idx <- all_kegg_df[["p.adjust"]] <= pcutoff
    sig_kegg_df <- all_kegg_df[sig_idx, ]
  }
  retlist <- list(
    "all_kegg" = all_kegg,
    "all_kegg_df" = all_kegg_df,
    "sig_kegg_df" = sig_kegg_df)
  class(retlist) <- "hpgltools::cp_kegg_enrich"
  return(retlist)
}

cp_kegg_gsea <- function(de_table_merged, kegg_universe, kegg_organism,
                         min_groupsize = 5, pcutoff = 0.05,
                         internal = FALSE, fc_column = "deseq_logfc") {
  lastcol <- ncol(de_table_merged)
  kegg_genelist <- as.vector(de_table_merged[[fc_column]])
  names(kegg_genelist) <- de_table_merged[[lastcol]]
  kegg_all_names <- glue("ncbi-geneid:{names(kegg_genelist)}")
  kegg_all_intersect <- kegg_all_names %in% names(kegg_universe)
  mesg("Found ", sum(kegg_all_intersect),
       " matches between the gene list and kegg universe.")
  ## all_names <- names(kegg_universe)
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
      geneList = kegg_genelist, organism = kegg_organism,
      minGSSize = min_groupsize, pvalueCutoff = 1.0, use_internal_data = internal))
  if ("try-error" %in% class(gse_all_kegg)) {
    return(NULL)
  } else {
    gse_all_kegg_df <- as.data.frame(gse_all_kegg, stringsAsFactors = FALSE)
    sig_idx <- gse_all_kegg_df[["p.adjust"]] <= pcutoff
    gse_sig_kegg_df <- gse_all_kegg_df[sig_idx, ]
    mesg("Found ", nrow(gse_sig_kegg_df), " KEGG GSE hits.")
  }
  retlist <- list(
    "gse_all_kegg" = gse_all_kegg,
    "gse_all_kegg_df" = gse_all_kegg_df,
    "gse_sig_kegg_df" = gse_sig_kegg_df)
  class(retlist) <- "hpgltools::cp_kegg_gsea"
  return(retlist)
}

cp_david_enrich <- function(sig_gene_list, min_groupsize = 5,
                            david_id = NULL, david_user = NULL, pcutoff = 0.05) {
  if (is.null(david_user)) {
    return(NULL)
  }
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
      david_all_df <- as.data.frame(david_all, stringsAsFactors = FALSE)
      sig_idx <- david_all_df[["p.adjust"]] <= pcutoff
      david_sig_df <- david_all_df[sig_idx, ]
    }
  }
  retlist <- list(
    "david_all" = david_all,
    "david_all_df" = david_all_df,
    "david_sig_df" = david_sig_df)
  class(retlist) <- "hpgltools::cp_david_enrich"
  return(retlist)
}

cp_react_organism <- function(orgdb) {
  orgdb_name <- orgdb
  if ("OrgDb" %in% class(orgdb)) {
    orgdb_name <- orgdb@.xData[["packageName"]]
  }
  if (orgdb_name == "org.Hs.eg.db") {
    reactome_organism <- "human"
  } else if (orgdb_name == "org.Mm.eg.db") {
    reactome_organism <- "mouse"
  } else {
    return(NULL)
  }
  loaded <- sm(try(requireNamespace(package = "ReactomePA", quietly = TRUE)))
  if ("try-error" %in% class(loaded)) {
    return(NULL)
  }
  return(reactome_organism)
}

cp_react_enrich <- function(sig_gene_list, orgdb, universe_to,
                            reactome_organism = NULL, min_groupsize = 5,
                            max_groupsize = 500, padj_type = "BH",
                            pcutoff = 0.05) {
  reactome_all <- NULL
  reactome_all_df <- data.frame()
  reactome_sig_df <- data.frame()
  mesg("Performing reactome over-representation analysis.")
  reactome_all <- try(ReactomePA::enrichPathway(
    gene = sig_gene_list, pvalueCutoff = pcutoff, readable = TRUE,
    pAdjustMethod = padj_type, qvalueCutoff = 1.0, universe = universe_to,
    organism = reactome_organism, minGSSize = min_groupsize, maxGSSize = max_groupsize))
  if ("try-error" %in% class(reactome_all)) {
    return(NULL)
  } else {
    reactome_all_df <- as.data.frame(reactome_all, stringsAsFactors = FALSE)
    sig_idx <- reactome_all_df[["p.adjust"]] <= pcutoff
    reactome_sig_df <- reactome_all_df[sig_idx, ]
    mesg("Found ", nrow(reactome_sig_df), " reactome over-represented hits.")
  }
  retlist <- list(
    "reactome_all" = reactome_all,
    "reactome_all_df" = reactome_all_df,
    "reactome_sig_df" = reactome_sig_df)
  class(retlist) <- "hpgltools::cp_react_enrich"
  return(retlist)
}

cp_react_gsea <- function(genelist, reactome_organism,
                          min_groupsize = 5, max_groupsize = 500, pcutoff = 0.05) {
  gse_all_reactome <- NULL
  gse_reactome_all_df <- data.frame()
  gse_reactome_sig_df <- data.frame()
  mesg("Performing reactome GSEA.")
  gse_all_reactome <- try(
    ReactomePA::gsePathway(geneList = genelist, organism = reactome_organism,
                           pvalueCutoff = 1.0, minGSSize = min_groupsize,
                           maxGSSize = max_groupsize))
  if ("try-error" %in% class(gse_all_reactome)) {
    return(NULL)
  } else {
    gse_reactome_all_df <- as.data.frame(gse_all_reactome, stringsAsFactors = FALSE)
    sig_idx <- gse_reactome_all_df[["p.adjust"]] <= pcutoff
    gse_reactome_sig_df <- gse_reactome_all_df[sig_idx, ]
    mesg("Found ", nrow(gse_reactome_sig_df), " Reactome GSE hits.")
  }
  retlist <- list(
    "gse_all_reactome" = gse_all_reactome,
    "gse_reactome_all_df" = gse_reactome_all_df,
    "gse_reactome_sig_df" = gse_reactome_sig_df)
  class(retlist) <- "hpgltools::cp_react_gsea"
  return(retlist)
}

cp_dose_organism <- function(orgdb) {
  orgdb_name <- orgdb
  dose_info <- list()
  if ("OrgDb" %in% class(orgdb)) {
    orgdb_name <- orgdb@.xData[["packageName"]]
  }
  if (orgdb_name == "org.Hs.eg.db") {
    dose_info[["orgn"]] <- "hsa"
    dose_info[["db"]] <- "HDO"
    loaded <- sm(try(requireNamespace(package = "HDO.db", quietly = TRUE)))
    if ("try-error" %in% class(loaded)) {
      return(NULL)
    }
  } else if (orgdb_name == "org.Mm.eg.db") {
    dose_info[["orgn"]] <- "mm"
    dose_info[["db"]] <- "MPO"
    loaded <- sm(try(requireNamespace(package = "MPO.db", quietly = TRUE)))
    if ("try-error" %in% class(loaded)) {
      return(NULL)
    }
  } else {
    return(NULL)
  }
  return(dose_info)
}

cp_dose_enrich <- function(sig_gene_list, dose_db, dose_orgn,
                           pcutoff = 0.05, padj_type = "BH",
                           universe_to = "ENTREZ", min_groupsize = 5,
                           max_groupsize = 500) {
  dose_all <- NULL
  dose_all_df <- dose_sig_df <- data.frame()
  mesg("Performing DOSE over-representation analysis.")
  dose_all <- try(DOSE::enrichDO(
    gene = sig_gene_list, ont = dose_db, organism = dose_orgn,
    pvalueCutoff = pcutoff, pAdjustMethod = padj_type, universe = universe_to,
    minGSSize = min_groupsize, maxGSSize = max_groupsize,
    qvalueCutoff = 1.0, readable = TRUE))
  if ("try-error" %in% class(dose_all)) {
    return(NULL)
  } else {
    dose_all_df <- as.data.frame(dose_all, stringsAsFactors = FALSE)
    sig_idx <- dose_all_df[["p.adjust"]] <= pcutoff
    dose_sig_df <- dose_all_df[sig_idx, ]
    mesg("Found ", nrow(dose_sig_df), " DOSE over-representation hits.")
  }
  retlist <- list(
    "dose_all" = dose_all,
    "dose_all_df" = dose_all_df,
    "dose_sig_df" = dose_sig_df)
  class(retlist) <- "hpgltools::cp_dose_enrich"
  return(retlist)
}

cp_dose_gsea <- function(genelist, dose_db, dose_orgn,
                         pcutoff = 0.05, padj_type = "BH",
                         min_groupsize = 5, max_groupsize = 500) {
  gse_dose <- NULL
  gse_dose_all_df <- gse_dose_sig_df <- data.frame()
  mesg("Performing DOSE GSEA.")
  gse_dose <- try(
    DOSE::gseDO(geneList = genelist, ont = dose_db, organism = dose_orgn,
                pvalueCutoff = 1.0, minGSSize = min_groupsize,
                maxGSSize = max_groupsize))
  if ("try-error" %in% class(gse_dose)) {
    return(NULL)
  } else {
    gse_dose_all_df <- as.data.frame(gse_dose, stringsAsFactors = FALSE)
    sig_idx <- gse_dose_all_df[["p.adjust"]] <= pcutoff
    gse_dose_sig_df <- gse_dose_all_df[sig_idx, ]
    mesg("Found ", nrow(gse_dose_sig_df), " Dose GSE hits.")
  }
  retlist <- list(
    "gse_dose" = gse_dose,
    "gse_dose_all_df" = gse_dose_all_df,
    "gse_dose_sig_df" = gse_dose_sig_df)
  class(retlist) <- "hpgltools::cp_react_gsea"
  return(retlist)
}

cp_mesh_organism <- function(orgdb) {
  orgdb_name <- orgdb
  mesh_info <- list(
    "mesh_category" = "C",
    "mesh_dbname" = "gendoo")
  if ("OrgDb" %in% class(orgdb)) {
    orgdb_name <- orgdb@.xData[["packageName"]]
  }
  if (orgdb_name == "org.Hs.eg.db") {
    mesh_info[["mesh_org"]] <- "Homo sapiens"
    loaded <- sm(try(requireNamespace(package = "MeSH.Hsa.eg.db", quietly = TRUE)))
    if ("try-error" %in% class(loaded)) {
      return(NULL)
    }
  } else if (orgdb_name == "org.Mm.eg.db") {
    mesh_info[["mesh_org"]] <- "Mus musculus"
    loaded <- sm(try(requireNamespace(package = "MeSH.Mm.eg.db", quietly = TRUE)))
    if ("try-error" %in% class(loaded)) {
      return(NULL)
    }
  } else {
    return(NULL)
  }
  loaded <- sm(try(requireNamespace(package = "AnnotationHub")))
  if ("try-error" %in% class(loaded)) {
    return(NULL)
  }
  loaded <- sm(try(requireNamespace(package = "meshes", quietly = TRUE)))
  if ("try-error" %in% class(loaded)) {
    return(NULL)
  }
  return(mesh_info)
}

cp_mesh_enrich <- function(sig_gene_list, mesh_db, mesh_org,
                           padj_type = "BH", universe_to = "ENTREZID",
                           min_groupsize = 5, max_groupsize = 500, qcutoff = 0.1,
                           pcutoff = 0.05) {
  mesh_all <- try(meshes::enrichMeSH(
    gene = sig_gene_list, MeSHDb = mesh_db, database = mesh_org,
    pvalueCutoff = 1.0, pAdjustMethod = padj_type, universe = universe_to,
    minGSSize = min_groupsize, maxGSSize = max_groupsize,
    qvalueCutoff = qcutoff), silent = TRUE)
  if ("try-error" %in% class(mesh_all)) {
    return(NULL)
  } else {
    mesh_all_df <- as.data.frame(mesh_all, stringsAsFactors = FALSE)
    sig_idx <- mesh_all_df[["p.adjust"]] <= pcutoff
    mesh_sig_df <- mesh_all_df[sig_idx, ]
    mesg("Found ", nrow(mesh_sig_df), " MESH over-representation hits.")
  }
  retlist <- list(
    "mesh_all" = mesh_all,
    "mesh_all_df" = mesh_all_df,
    "mesh_sig_df" = mesh_sig_df)
  class(retlist) <- "hpgltools::cp_mesh_enrich"
  return(retlist)
}

cp_mesh_gsea <- function(genelist, mesh_db, mesh_dbname,
                         pcutoff = pcutoff, min_groupsize = 5, max_groupsize = 500) {
  gse_mesh <- try(
    meshes::gseMeSH(
      geneList = genelist, MeSHDb = mesh_db, database = mesh_dbname,
      pvalueCutoff = 1.0, minGSSize = min_groupsize,
      maxGSSize = max_groupsize))
  if ("try-error" %in% class(gse_mesh)) {
    return(NULL)
  } else {
    gse_mesh_all_df <- as.data.frame(gse_mesh, stringsAsFactors = FALSE)
    sig_idx <- gse_mesh_all_df[["p.adjust"]] <= pcutoff
    gse_mesh_sig_df <- gse_mesh_all_df[sig_idx, ]
    mesg("Found ", nrow(gse_mesh_sig_df), " MESH GSE hits.")
  }
  retlist <- list(
    "gse_mesh" = gse_mesh,
    "gse_mesh_all_df" = gse_mesh_all_df,
    "gse_mesh_sig_df" = gse_mesh_sig_df)
  class(retlist) <- "hpgltools::cp_mesh_gsea"
  return(retlist)
}

cp_msigdb_loaded <- function(signature_data, msig_db, msigdb_category,
                             id_type = "entrez") {
  if (is.null(signature_data)) {
    signature_data <- try(load_gmt_signatures(
      signatures = msig_db, signature_category = msigdb_category, id_type = "entrez"))
    if ("try-error" %in% class(signature_data)) {
      return(NULL)
    }
  }
  if (is.null(signature_df)) {
    signature_df <- try(signatures_to_df(signature_data))
    if ("try-error" %in% class(signature_df)) {
      return(NULL)
    }
  }
  signatures <- list(
    "signature_data" = signature_data,
    "signature_df" = signature_df)
  return(signatures)
}

cp_msigdb_enrich <- function(sig_gene_list, signature_data = NULL, signature_df = NULL,
                             org = NULL, pcutoff = 0.05) {
  msigdb_all <- NULL
  msigdb_all_df <- msigdb_sig_df <- data.frame()

  ## I think I was planning to switch IDs here but never finished the logic.
  ## expected_term_type <- "ENTREZID"
  ## current_term_type <- "ENTREZID"
  ## if (current_term_type != expected_term_type) {
  ##   test_sig_df <- sm(try(clusterProfiler::bitr(
  ##     sig_gene_list, fromType = current_term_type, toType = expected_term_type, OrgDb = org),
  ##     silent = TRUE))
  ## }

  mesg("Starting msigDB over-representation analysis.")
  msigdb_all <- try(clusterProfiler::enricher(
    sig_gene_list, TERM2GENE = signature_df))
  if ("try-error" %in% class(msigdb_all)) {
    return(NULL)
  } else {
    msigdb_all_df <- as.data.frame(msigdb_all, stringsAsFactors = FALSE)
    sig_idx <- msigdb_all_df[["p.adjust"]] <= pcutoff
    msigdb_sig_df <- msigdb_all_df[sig_idx, ]
    mesg("Found ", nrow(msigdb_sig_df), " mSig GSE hits.")
  }
  retlist <- list(
    "msigdb_all" = msigdb_all,
    "msigdb_all_df" = msigdb_all_df,
    "msigdb_sig_df" = msigdb_sig_df)
  class(retlist) <- "hpgltools::cp_msigdb_enrich"
  return(retlist)
}

cp_msigdb_gsea <- function(genelist, signature_data, signature_df, org,
                           min_groupsize = 5, max_groupsize = 500,
                           padj_type = "BH", pcutoff = 0.05) {
  gse_msigdb_all <- NULL
  gse_msigdb_all_df <- gse_msigdb_sig_df <- data.frame()
  mesg("Starting msigDB GSEA.")
  gse_msigdb_all <- try(
    clusterProfiler::GSEA(
      genelist, exponent = 1, minGSSize = min_groupsize, maxGSSize = max_groupsize,
      eps = 1e-10, pvalueCutoff = 1.0, pAdjustMethod = padj_type,
      TERM2GENE = signature_df))
  if ("try-error" %in% class(gse_msigdb_all)) {
    return(NULL)
  } else {
    gse_msigdb_all_df <- as.data.frame(gse_msigdb_all, stringsAsFactors = FALSE)
    sig_idx <- gse_msigdb_all_df[["p.adjust"]] <= pcutoff
    gse_msigdb_sig_df <- gse_msigdb_all_df[sig_idx, ]
    mesg("Found ", nrow(gse_msigdb_sig_df), " msigDB GSEA hits.")
  }
  retlist <- list(
    "gse_msigdb_all" = gse_msigdb_all,
    "gse_msigdb_all_df" = gse_msigdb_all_df,
    "gse_msigdb_sig_df" = gse_msigdb_sig_df)
  class(retlist) <- "hpgltools::cp_msigdb_gsea"
  return(retlist)
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
    "orgdb_sig_from" = orgdb_sig_from)
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
  ## TODO: Separate out these various pieces of functionality to query
  ## stuff like kegg/reactome/etc into separate functions.
  sm(requireNamespace(package = "clusterProfiler", quietly = TRUE))
  sm(requireNamespace(package = "DOSE", quietly = TRUE))
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
  sig_genenames <- rownames(sig_genes)

  mapper_keys <- AnnotationDbi::keytypes(org)
  org_keys <- AnnotationDbi::keys(org)
  all_genenames <- c()
  if (is.null(de_table)) {
    all_genenames <- org_keys
  } else {
    all_genenames <- rownames(de_table)
  }

  if (is.null(orgdb_from)) {
    mesg("Guessing appropriate keytype for the orgdb.")
    key_info <- guess_bitr_keytype(org, orgdb_from, sig_genes = sig_genes, to = orgdb_to)
    ## orgdb_sig_from <- key_info[["orgdb_sig_from"]]
    de_table_namedf <- key_info[["gene_df"]]
    sig_genes_namedf <- key_info[["sig_df"]]
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
  }  ## End of the if cascade starting with testing orgdb_from

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
    warning("No genes were found between the significant genes and the universe.")
    return(NULL)
  }

  ## A final sanity test to see if it is possible to do GSE analyses
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

  ## kegg organism sanity check
  kegg_organism <- cp_kegg_organism(organism, orgdb, kegg_organism = kegg_organism)
  if (is.null(kegg_organism)) {
    todo[["enrich_kegg"]] <- FALSE
    todo[["gse_kegg"]] <- FALSE
  }
  ## Kegg gene list sanity check
  kegg_sig_ids <- NULL
  if (isTRUE(todo[["gse_kegg"]])) {
    kegg_info <- cp_kegg_universe(kegg_organism, sig_gene_list, kegg_universe = kegg_universe)
    if (is.null(kegg_info)) {
      todo[["gse_kegg"]] <- FALSE
    } else {
      kegg_universe <- kegg_info[["kegg_universe"]]
      kegg_sig_ids <- kegg_info[["kegg_sig_ids"]]
    }
  }

  ## Reactome sanity check
  reactome_organism <- cp_react_organism(orgdb)
  if (is.null(reactome_organism)) {
    todo[["enrich_reactome"]] <- FALSE
    todo[["gse_reactome"]] <- FALSE
  }

  ## DOSE sanity check
  dose_orgn <- NULL
  dose_db <- NULL
  dose_db_info <- cp_dose_organism(orgdb)
  if (is.null(dose_db_info)) {
    todo[["enrich_dose"]] <- FALSE
    todo[["gse_dose"]] <- FALSE
  } else {
    dose_orgn <- dose_db_info[["orgn"]]
    dose_db <- dose_db_info[["db"]]
  }

  ## mesh sanity check
  mesh_org <- NULL
  mesh_category <- NULL
  mesh_dbname <- NULL
  mesh_info <- cp_mesh_organism(orgdb)
  if (is.null(mesh_info)) {
    todo[["enrich_mesh"]] <- FALSE
    todo[["gse_mesh"]] <- FALSE
  } else {
    mesh_org <- mesh_info[["mesh_org"]]
    mesh_category <- mesh_info[["mesh_category"]]
    mesh_dbname <- mesh_info[["mesh_dbname"]]
  }

  ## msigdb sanity check
  signature_info <- cp_msigdb_loaded(signature_data, msig_db, msigdb_category)
  if (is.null(signature_info)) {
    todo[["enrich_msigdb"]] <- FALSE
    todo[["gse_msigdb"]] <- FALSE
  } else {
    signature_data <- signature_info[["signature_data"]]
    signature_df <- signature_info[["signature_df"]]
  }

  ## Now we have a universe of geneIDs and significant IDs
  ## Let us perform some analyses...
  go_data <- list()
  if (isTRUE(todo[["enrich_go"]])) {
    go_enrich_info <- cp_go_enrich(sig_gene_list, org, orgdb_to, go_level, universe_to,
                                   min_groupsize, padj_type, pcutoff = pcutoff)
    if (!is.null(go_enrich_info)) {
      go_data[["MF_all"]] <- go_enrich_info[["ego_all_mf_df"]]
      go_data[["MF_sig"]] <- go_enrich_info[["ego_sig_mf_df"]]
      go_data[["BP_all"]] <- go_enrich_info[["ego_all_bp_df"]]
      go_data[["BP_sig"]] <- go_enrich_info[["ego_sig_bp_df"]]
      go_data[["CC_all"]] <- go_enrich_info[["ego_all_cc_df"]]
      go_data[["CC_sig"]] <- go_enrich_info[["ego_sig_cc_df"]]
      go_data[["MF_enrich"]] <- go_enrich_info[["ego_all_mf"]]
      go_data[["BP_enrich"]] <- go_enrich_info[["ego_all_bp"]]
      go_data[["CC_enrich"]] <- go_enrich_info[["ego_all_cc"]]
    }
  }

  if (isTRUE(todo[["gse_go"]])) {
    go_gsea_info <- cp_go_gsea(genelist, org, orgdb_to, min_groupsize = min_groupsize,
                               pcutoff = pcutoff)
    if (!is.null(go_gsea_info)) {
      go_data[["GO_gse"]] <- go_gsea_info[["gse_go"]]
      go_data[["GO_gse_all"]] <- go_gsea_info[["gse_go_all_df"]]
      go_data[["GO_gse_sig"]] <- go_gsea_info[["gse_go_sig_df"]]
    }
  }
  ## This is foolishly written I know, but I am going to come back through and clean it asap
  ## until then I want to be able to mess with these on a per-instance basis

  kegg_data <- list()
  if (isTRUE(todo[["enrich_kegg"]])) {
    kegg_info <- cp_kegg_enrich(kegg_sig_ids, kegg_organism = kegg_organism, pcutoff = pcutoff)
    if (!is.null(kegg_info)) {
      kegg_data <- kegg_info
    }
  }

  if (isTRUE(todo[["gse_kegg"]])) {
    gse_kegg_data <- cp_kegg_gsea(de_table_merged, kegg_universe, kegg_organism,
                                  min_groupsize = min_groupsize, pcutoff = pcutoff,
                                  internal = internal, fc_column = fc_column)
    if (!is.null(gse_kegg_data)) {
      kegg_data[["gse_all_kegg"]] <- gse_kegg_data[["gse_all_kegg"]]
      kegg_data[["gse_all_kegg_df"]] <- gse_kegg_data[["gse_all_kegg_df"]]
      kegg_data[["gse_sig_kegg_df"]] <- gse_kegg_data[["gse_sig_kegg_df"]]
    }
  }

  ## DAVID requires a java interface which is prone to segmentation faults.
  ## I may revisit this in the future, so I will leave it here, but disabled.
  todo[["enrich_david"]] <- FALSE
  david_data <- list()
  if (isTRUE(todo[["enrich_david"]])) {
    david_data <- cp_david_enrich(sig_gene_list, min_groupsize = min_groupsize,
                                  david_id = david_id, david_user = david_user)
  }

  reactome_data <- list()
  if (isTRUE(todo[["enrich_reactome"]])) {
    reactome_info <- cp_react_enrich(sig_gene_list, orgdb, universe_to,
                                     reactome_organism = reactome_organism,
                                     min_groupsize = min_groupsize, max_groupsize = max_groupsize,
                                     padj_type = padj_type, pcutoff = pcutoff)
    if (!is.null(reactome_info)) {
      reactome_data[["reactome_all"]] <- reactome_info[["reactome_all"]]
      reactome_data[["reactome_all_df"]] <- reactome_info[["reactome_all_df"]]
      reactome_data[["reactome_sig_df"]] <- reactome_info[["reactome_sig_df"]]
    }
  }

  if (isTRUE(todo[["gse_reactome"]])) {
    reactome_info <- cp_react_gsea(genelist, reactome_organism,
                                   min_groupsize = min_groupsize, max_groupsize = max_groupsize,
                                   pcutoff = pcutoff)
    if (!is.null(reactome_info)) {
      reactome_data[["gse_all_reactome"]] <- reactome_info[["gse_all_reactome"]]
      reactome_data[["gse_reactome_all_df"]] <- reactome_info[["gse_reactome_all_df"]]
      reactome_data[["gse_reactome_sig_df"]] <- reactome_info[["gse_reactome_sig_df"]]
    }
  }


  dose_data <- list()
  if (isTRUE(todo[["enrich_dose"]])) {
    dose_info <- cp_dose_enrich(sig_gene_list, dose_db, dose_orgn,
                                pcutoff = pcutoff, padj_type = padj_type,
                                universe_to = universe_to, min_groupsize = min_groupsize,
                                max_groupsize = max_groupsize)
    if (!is.null(dose_info)) {
      dose_data[["dose_all"]] <- dose_info[["dose_all"]]
      dose_data[["dose_all_df"]] <- dose_info[["dose_all_df"]]
      dose_data[["dose_sig_df"]] <- dose_info[["dose_sig_df"]]
    }
  }

  if (isTRUE(todo[["gse_dose"]])) {
    mesg("Performing DOSE GSEA.")
    dose_info <- cp_dose_gsea(genelist, dose_db, dose_orgn,
                              pcutoff = pcutoff, padj_type = padj_type,
                              min_groupsize = min_groupsize, max_groupsize = max_groupsize)
    if (!is.null(dose_info)) {
      dose_data[["gse_dose"]] <- dose_info[["gse_dose"]]
      dose_data[["gse_dose_all_df"]] <- dose_info[["gse_dose_all_df"]]
      dose_data[["gse_dose_sig_df"]] <- dose_info[["gse_dose_sig_df"]]
    }
  }

  mesh_data <- list()
  if (isTRUE(todo[["enrich_mesh"]])) {
    mesh_info <- cp_mesh_enrich(sig_gene_list, mesh_db, mesh_org,
                                padj_type = padj_type, universe_to = universe_to,
                                min_groupsize = min_groupsize, max_groupsize,
                                qcutoff = qcutoff, pcutoff = pcutoff)
    if (!is.null(mesh_info)) {
      mesh_data[["mesh_all"]] <- mesh_info[["mesh_all"]]
      mesh_data[["mesh_all_df"]] <- mesh_info[["mesh_all_df"]]
      mesh_data[["mesh_sig_df"]] <- mesh_info[["mesh_sig_df"]]
    }
  }

  if (isTRUE(todo[["gse_mesh"]])) {
    mesh_info <- cp_mesh_gsea(genelist, mesh_db, mesh_dbname,
                              pcutoff = pcutoff, min_groupsize = min_groupsize,
                              max_groupsize = max_groupsize)
    if (!is.null(mesh_info)) {
      mesh_data[["gse_mesh"]] <- mesh_info[["gse_mesh"]]
      mesh_data[["gse_mesh_all_df"]] <- mesh_info[["gse_mesh_all_df"]]
      mesh_data[["gse_mesh_sig_df"]] <- mesh_info[["gse_mesh_sig_df"]]
    }
  }

  msigdb_data <- list()
  if (isTRUE(todo[["enrich_msigdb"]])) {
    msigdb_info <- cp_msigdb_enrich(sig_gene_list, signature_data, signature_df, org,
                                    pcutoff = pcutoff)
    if (!is.null(msigdb_info)) {
      msigdb_data[["msigdb_all"]] <- msigdb_info[["msigdb_all"]]
      msigdb_data[["msigdb_all_df"]] <- msigdb_info[["msigdb_all_df"]]
      msigdb_data[["msigdb_sig_df"]] <- msigdb_info[["msigdb_sig_df"]]
    }
  }

  if (isTRUE(todo[["gse_msigdb"]])) {
    msigdb_info <- cp_msigdb_gsea(genelist, signature_data, signature_df, org,
                                  min_groupsize = min_groupsize, max_groupsize = max_groupsize,
                                  padj_type = padj_type, pcutoff = pcutoff)
    if (!is.null(msigdb_info)) {
      msigdb_data[["gse_msigdb_all"]] <- msigdb_info[["gse_msigdb_all"]]
      msigdb_data[["gse_msigdb_all_df"]] <- msigdb_info[["gse_msigdb_all_df"]]
      msigdb_data[["gse_msigdb_sig_df"]] <- msigdb_info[["gse_msigdb_sig_df"]]
    }
  }

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
  class(retlist) <- "hpgltools::simple_clusterprofiler"
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
  ## all_genenames <- rownames(de_table)
  ## sig_genenames <- rownames(sig_genes)
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
