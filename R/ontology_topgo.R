
#' Coerce a df of gene IDs and GO categories to the topGO list
#'
#' The read_mappings function in topGO produces a list of x element
#' where x is named by the gene IDs and the value of each element
#' is a character vector comprised of the categories associated with it.
#' @param go_db The df used by goseq and most commonly by me.
df_to_mappings <- function(go_db) {
  gene_ids <- unique(go_db[["ID"]])
  go_lst <- list()
  for (id in seq_along(gene_ids)) {
    id_name <- gene_ids[id]
    id_vector <- c()
    categories <- go_db[["ID"]] == id_name
    category_names <- go_db[categories, "GO"]
    if (sum(categories) > 0) {
      vector <- go_db[categories, "GO"]
      na <- is.na(vector)
      vector[!na]
      empty <- vector == ""
      vector <- vector[!empty]
      if (length(vector) == 0) {
        next
      }
      go_lst[[id_name]] <- vector
    }
  } ## End iterating over genes
  return(go_lst)
}

#' Perform a simplified topgo analysis.
#'
#' This will attempt to make it easier to run topgo on a set of genes.
#' The way I organized these data structures is completely stupid. I
#' want to convert the data from topgo to clusterprofiler for ease of
#' plotting, but because of the terrible way I organized everything
#' that is likely to be difficult.
#'
#' @param sig_genes Data frame of differentially expressed genes, containing IDs
#'  any other columns.
#' @param goid_map File containing mappings of genes to goids in the format
#'  expected by topgo.
#' @param go_db Data frame of the goids which may be used to make the goid_map.
#' @param pvals Set of pvalues in the DE data which may be used to improve the
#'  topgo results.
#' @param limitby Test to index the results by.
#' @param limit Ontology pvalue to use as the lower limit.
#' @param signodes I don't remember right now.
#' @param sigforall Provide the significance for all nodes?
#' @param numchar Character limit for the table of results.
#' @param selector Function name for choosing genes to include.
#' @param pval_column Column from which to acquire scores.
#' @param overwrite Yeah I do not remember this one either.
#' @param densities Densities, yeah, the densities...
#' @param pval_plots Include pvalue plots of the results a la clusterprofiler?
#' @param parallel Use doParallel?
#' @param excel  Print the results to an excel file?
#' @param ... Other options which I do not remember right now!
#' @return Big list including the various outputs from topgo
#' @seealso [topGO]
#' @export
simple_topgo <- function(sig_genes, goid_map = "id2go.map", go_db = NULL,
                         pvals = NULL, limitby = "fisher", limit = 0.1, signodes = 100,
                         sigforall = TRUE, numchar = 300, selector = "topDiffGenes",
                         pval_column = "deseq_adjp", overwrite = FALSE, densities = TRUE,
                         pval_plots = TRUE, parallel = TRUE, excel = NULL, ...) {
  ## Some neat ideas from the topGO documentation:
  ## geneList <- getPvalues(exprs(eset), classlabel = y, alternative = "greater")
  ## A variant of these operations make it possible to give topGO scores so that
  ## a larger array of tests may be performed
  ## x <- topDiffGenes(geneList)
  ## sum(x) ## the number of selected genes
  ## If we do something like above to give scores to all the 'DEgenes', then we
  ## set up the GOdata object like this: mf_GOdata = new("topGOdata",
  ## description = "something", ontology = "BP", allGenes = entire_geneList,
  ## geneSel = topDiffGenes, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 2)

  ## What I should do: create the data structure returned by readMappings()
  ## Creating this file denovo is dumb
  if (isTRUE(overwrite) && file.exists(goid_map)) {
    removed <- file.remove(goid_map)
  }
  #gomap_info <- make_id2gomap(goid_map = goid_map, go_db = go_db, overwrite = overwrite)
  #if (is.null(gomap_info)) {
  #  warning("There appears to have been a problem generating the gomap.")
  #}
  #geneID2GO <- topGO::readMappings(file = goid_map)
  geneID2GO <- df_to_mappings(go_db)
  annotated_genes <- names(geneID2GO)
  if (is.null(go_db)) {
    go_db <- reshape2::melt(geneID2GO)
    go_db <- go_db[, c("L1", "value")]
    colnames(go_db) <- c("ID", "GO")
  }

  if (is.null(sig_genes[["ID"]])) {
    sig_genes[["ID"]] <- make.names(rownames(sig_genes), unique = TRUE)
  }
  fisher_interesting_genes <- as.factor(
      as.integer(annotated_genes %in% sig_genes[["ID"]]))
  names(fisher_interesting_genes) <- annotated_genes
  ks_interesting_genes <- as.integer(!annotated_genes %in% sig_genes[["ID"]])
  if (!is.null(sig_genes[[pval_column]])) {
    ## I think this needs to include the entire gene universe, not only the set
    ## of x differentially expressed genes, making this an explicit
    ## as.vector(as.numeric()) because it turns out the values from DESeq are
    ## characters.
    pvals <- as.vector(as.numeric(sig_genes[[pval_column]]))
    names(pvals) <- rownames(sig_genes)
    for (p in seq_along(pvals)) {
      name <- names(pvals)[p]
      ks_interesting_genes[[name]] <- pvals[p]
    }
  }
  ks_interesting_genes <- as.vector(ks_interesting_genes)
  names(ks_interesting_genes) <- annotated_genes

  ## Instead of invoking library(topGO), I can requireNamespace && attachNamespace.
  ## "GOMFTerm not found"
  ## Ergo, requireNamespace() is insufficient!

  ontology_result <- list(
      "BP" = list(),
      "MF" = list(),
      "CC" = list())
  returns <- list()
  methods <- c("fisher", "KS", "EL", "weight")
  if (isTRUE(parallel)) {
    cl <- parallel::makeCluster(4)
    doParallel::registerDoParallel(cl)
    tt <- sm(requireNamespace("parallel"))
    tt <- sm(requireNamespace("doParallel"))
    tt <- sm(requireNamespace("iterators"))
    tt <- sm(requireNamespace("foreach"))
    res <- foreach(c = seq_along(methods),
                 ## .combine = "c", .multicombine = TRUE,
                 .packages = c("hpgltools", "Hmisc", "topGO")) %dopar% {
                   type <- methods[c]
                   returns[[type]] <- do_topgo(type, go_map = geneID2GO,
                                               fisher_genes = fisher_interesting_genes,
                                               ks_genes = ks_interesting_genes)
    }
    stopped <- parallel::stopCluster(cl)
    if (class(stopped)[1] == "try-error") {
      warning("There was a problem stopping the parallel cluster.")
    }
    for (r in seq_along(methods)) {
      a_result <- res[[r]]
      type <- a_result[["MF"]][["type"]]
      ontology_result[["MF"]][[type]] <- a_result[["MF"]]
      ontology_result[["BP"]][[type]] <- a_result[["BP"]]
      ontology_result[["CC"]][[type]] <- a_result[["CC"]]
    }
  } else {
    for (r in seq_along(methods)) {
      type <- methods[r]
      message("Starting ", type, ".")
      a_result <- do_topgo(type, go_map = geneID2GO,
                           fisher_genes = fisher_interesting_genes,
                           ks_genes = ks_interesting_genes)
      ontology_result[["MF"]][[type]] <- a_result[["MF"]]
      ontology_result[["BP"]][[type]] <- a_result[["BP"]]
      ontology_result[["CC"]][[type]] <- a_result[["CC"]]
    }
  }

  included_methods <- names(ontology_result[["MF"]])
  p_dists <- list()
  for (o in c("BP", "MF", "CC")) {
    for (m in included_methods) {
      name <- glue::glue("{tolower(o)}_{limitby}")
      p_dists[[name]] <- try(plot_histogram(
          ontology_result[[o]][[m]][["test_result"]]@score,
          bins = 20))
    }
  }

  godata <- list()
  results <- list()
  if (!is.null(ontology_result[["MF"]][["fisher"]])) {
    godata[["fisher_bp"]] = ontology_result[["BP"]][["fisher"]][["godata"]]
    godata[["fisher_mf"]] = ontology_result[["MF"]][["fisher"]][["godata"]]
    godata[["fisher_cc"]] = ontology_result[["CC"]][["fisher"]][["godata"]]
    results[["fisher_bp"]] = ontology_result[["BP"]][["fisher"]][["test_result"]]
    results[["fisher_mf"]] = ontology_result[["MF"]][["fisher"]][["test_result"]]
    results[["fisher_cc"]] = ontology_result[["CC"]][["fisher"]][["test_result"]]
  }

  if (!is.null(ontology_result[["MF"]][["KS"]])) {
    godata[["ks_bp"]] = ontology_result[["BP"]][["KS"]][["godata"]]
    godata[["ks_mf"]] = ontology_result[["MF"]][["KS"]][["godata"]]
    godata[["ks_cc"]] = ontology_result[["CC"]][["KS"]][["godata"]]
    results[["ks_bp"]] = ontology_result[["BP"]][["KS"]][["test_result"]]
    results[["ks_mf"]] = ontology_result[["MF"]][["KS"]][["test_result"]]
    results[["ks_cc"]] = ontology_result[["CC"]][["KS"]][["test_result"]]
    results[["el_bp"]] = ontology_result[["BP"]][["EL"]][["test_result"]]
  }

  if (!is.null(ontology_result[["MF"]][["EL"]])) {
    godata[["el_bp"]] = ontology_result[["BP"]][["EL"]][["godata"]]
    godata[["el_mf"]] = ontology_result[["MF"]][["EL"]][["godata"]]
    godata[["el_cc"]] = ontology_result[["CC"]][["EL"]][["godata"]]
    results[["el_mf"]] = ontology_result[["MF"]][["EL"]][["test_result"]]
    results[["el_bp"]] = ontology_result[["BP"]][["EL"]][["test_result"]]
    results[["el_cc"]] = ontology_result[["CC"]][["EL"]][["test_result"]]
  }

  if (!is.null(ontology_result[["MF"]][["weight"]])) {
    godata[["weight_bp"]] = ontology_result[["BP"]][["weight"]][["godata"]]
    godata[["weight_mf"]] = ontology_result[["MF"]][["weight"]][["godata"]]
    godata[["weight_cc"]] = ontology_result[["CC"]][["weight"]][["godata"]]
    results[["weight_bp"]] = ontology_result[["BP"]][["weight"]][["test_result"]]
    results[["weight_mf"]] = ontology_result[["MF"]][["weight"]][["test_result"]]
    results[["weight_cc"]] = ontology_result[["CC"]][["weight"]][["test_result"]]
  }

  tables <- try(topgo_tables(results, godata, limitby = limitby, limit = limit,
                             orderby = limitby, ranksof = limitby))
  if (class(tables)[1] == "try-error") {
    tables <- NULL
  }

  mf_densities <- bp_densities <- cc_densities <- list()
  if (isTRUE(densities)) {
    bp_densities <- sm(
      plot_topgo_densities(godata[["fisher_bp"]], tables[["bp_interesting"]]))
    mf_densities <- sm(
      plot_topgo_densities(godata[["fisher_mf"]], tables[["mf_interesting"]]))
    cc_densities <- sm(
      plot_topgo_densities(godata[["fisher_cc"]], tables[["cc_interesting"]]))
  } else {
    message("simple_topgo(): Set densities = TRUE for ontology density plots.")
  }

  retlist <- list(
      "go_db" = go_db,
      "input" = sig_genes,
      "godata" = godata,
      "results" = results,
      "tables" = tables,
      "mf_densities" = mf_densities,
      "bp_densities" = bp_densities,
      "cc_densities" = cc_densities,
      "pdists" = p_dists)

  enrich_results <- list()
  for (ont in c("bp", "mf", "cc")) {
    message("Getting enrichResult for ontology: ", ont, ".")
    enrich_results[[ont]] <- topgo2enrich(retlist, ontology = ont, pval = limit,
                                          column = limitby)
  }
  retlist[["enrich_results"]] <- enrich_results

  pval_histograms <- list()
  fisher_ps <- c(retlist[["tables"]][["mf_subset"]][["fisher"]],
                 retlist[["tables"]][["bp_subset"]][["fisher"]],
                 retlist[["tables"]][["cc_subset"]][["fisher"]])
  pval_histograms[["fisher"]] <- sm(try(plot_histogram(fisher_ps, bins = 50))) +
    ggplot2::ylab("Number of ontologies observed.") +
    ggplot2::xlab("Fisher exact test score.")
  ks_ps <- c(retlist[["tables"]][["mf_subset"]][["KS"]],
             retlist[["tables"]][["bp_subset"]][["KS"]],
             retlist[["tables"]][["cc_subset"]][["KS"]])
  pval_histograms[["KS"]] <- sm(try(plot_histogram(ks_ps, bins = 50))) +
    ggplot2::ylab("Number of ontologies observed.") +
    ggplot2::xlab("KS test score.")
  el_ps <- c(retlist[["tables"]][["mf_subset"]][["EL"]],
             retlist[["tables"]][["bp_subset"]][["EL"]],
             retlist[["tables"]][["cc_subset"]][["EL"]])
  pval_histograms[["EL"]] <- sm(try(plot_histogram(el_ps, bins = 50))) +
    ggplot2::ylab("Number of ontologies observed.") +
    ggplot2::xlab("EL test score.")
  weight_ps <- c(retlist[["tables"]][["mf_subset"]][["weight"]],
                 retlist[["tables"]][["bp_subset"]][["weight"]],
                 retlist[["tables"]][["cc_subset"]][["weight"]])
  pval_histograms[["weight"]] <- sm(try(plot_histogram(weight_ps, bins = 50))) +
    ggplot2::ylab("Number of ontologies observed.") +
    ggplot2::xlab("Weighted test score.")
  retlist[["pvalue_histograms"]] <- pval_histograms

  if (!is.null(excel)) {
    message("Writing data to: ", excel, ".")
    excel_ret <- sm(try(write_topgo_data(retlist, excel = excel)))
    retlist[["excel"]] <- excel_ret
  }
  retlist[["goid_map"]] <- goid_map
  class(retlist) <- "topgo_result"
  return(retlist)
}
setGeneric("simple_topgo")

#' Coerce simple_topgo to accept a vector of gene IDs instead of a real dataframe of significance.
#'
#' Doing this voids the topgo warantee.
#'
#' @param sig_genes Data frame of differentially expressed genes, containing IDs
#'  any other columns.
#' @param goid_map File containing mappings of genes to goids in the format
#'  expected by topgo.
#' @param go_db Data frame of the goids which may be used to make the goid_map.
#' @param pvals Set of pvalues in the DE data which may be used to improve the
#'  topgo results.
#' @param limitby Test to index the results by.
#' @param limit Ontology pvalue to use as the lower limit.
#' @param signodes I don't remember right now.
#' @param sigforall Provide the significance for all nodes?
#' @param numchar Character limit for the table of results.
#' @param selector Function name for choosing genes to include.
#' @param pval_column Column from which to acquire scores.
#' @param overwrite Yeah I do not remember this one either.
#' @param densities Densities, yeah, the densities...
#' @param pval_plots Include pvalue plots of the results a la clusterprofiler?
#' @param parallel Use doParallel?
#' @param excel  Print the results to an excel file?
#' @param ... Other options which I do not remember right now!
setMethod(
  "simple_topgo", signature = signature(sig_genes = "character"),
  definition = function(sig_genes, goid_map = "id2go.map", go_db = NULL,
                        pvals = NULL, limitby = "fisher", limit = 0.1,
                        signodes = 100, sigforall = TRUE, numchar = 300,
                        selector = "topDiffGenes", pval_column = "deseq_adjp",
                        overwrite = FALSE, densities = FALSE,
                        pval_plots = TRUE, excel = NULL, ...) {
    fake_df <- data.frame(row.names = sig_genes)
    fake_df[["ID"]] <- rownames(fake_df)
    fake_df[[pval_column]] <- 0.01
    warning("Faking a dataframe with significance of every gene as 0.01 because this was given a vector of gene IDs.")
    simple_topgo(fake_df, goid_map = goid_map, go_db = go_db,
                 pvals = pvals, limitby = limitby, limit = limit,
                 signodes = signodes, sigforall = sigforall,
                 numchar = numchar, selector = selector,
                 pval_column = pval_column, overwrite = overwrite,
                 densities = densities, pval_plots = pval_plots, excel = excel,
                 ...)
  })

#' An attempt to make topgo invocations a bit more standard.
#'
#' My function 'simple_topgo()' was excessively long and a morass of copy/pasted
#' fragments.  This attempts to simplify that and converge on a single piece
#' of code for all the methodologies provided by topgo.
#'
#' @param type Type of topgo search to perform: fisher, KS, EL, or weight.
#' @param go_map Mappings of gene and GO IDs.
#' @param fisher_genes List of genes used for fisher analyses.
#' @param ks_genes List of genes used for KS analyses.
#' @param selector Function to use when selecting genes.
#' @param sigforall Provide significance metrics for all ontologies observed,
#'  not only the ones deemed statistically significant.
#' @param numchar A limit on characters printed when printing topgo tables
#'  (used?)
#' @param pval_column Column from which to extract DE p-values.
#' @param overwrite Overwrite an existing gene ID/GO mapping?
#' @param cutoff Define 'significant'?
#' @param densities Perform gene density plots by ontology?
#' @param pval_plots Print p-values plots as per clusterProfiler?
#' @return List of results from the various tests in topGO.
#' @seealso [topGO]
#' @export
do_topgo <- function(type, go_map = NULL, fisher_genes = NULL, ks_genes = NULL,
                     selector = "topDiffGenes", sigforall = TRUE, numchar = 300,
                     pval_column = "adj.P.Val", overwrite = FALSE,
                     cutoff = 0.05, densities = FALSE, pval_plots = TRUE) {
  tt <- try(sm(requireNamespace("topGO")), silent = TRUE)
  tt <- try(sm(attachNamespace("topGO")), silent = TRUE)
  retlist <- list(
      "BP" = list(),
      "MF" = list(),
      "CC" = list())
  for (ont in names(retlist)) {
    switchret <- switch(
        type,
        "fisher" = {
          retlist[[ont]][["type"]] <- "fisher"
          retlist[[ont]][["godata"]] <- new(
              "topGOdata",
              description = ont,
              ontology = ont,
              allGenes = fisher_genes,
              annot = topGO::annFUN.gene2GO,
              gene2GO = go_map)
          retlist[[ont]][["test_stat"]] <- new(
              "classicCount",
              testStatistic = topGO::GOFisherTest,
              name = "Fisher test")
          retlist[[ont]][["test_result"]] <- topGO::getSigGroups(
                                                        retlist[[ont]][["godata"]],
                                                        retlist[[ont]][["test_stat"]])
        },
        "KS" = {
          retlist[[ont]][["type"]] <- "KS"
          retlist[[ont]][["godata"]] <- new(
              "topGOdata",
              description = ont,
              ontology = ont,
              allGenes = ks_genes,
              geneSel = get(selector),
              annot = topGO::annFUN.gene2GO,
              gene2GO = go_map)
          retlist[[ont]][["test_stat"]] <- new(
              "classicScore",
              testStatistic = topGO::GOKSTest,
              name = "KS test")
          retlist[[ont]][["test_result"]] <- topGO::getSigGroups(
                                                        retlist[[ont]][["godata"]],
                                                        retlist[[ont]][["test_stat"]])
        },
        "EL" = {
          retlist[[ont]][["type"]] <- "EL"
          retlist[[ont]][["godata"]] <- new(
              "topGOdata",
              description = ont,
              ontology = ont,
              allGenes = ks_genes,
              geneSel = get(selector),
              annot = topGO::annFUN.gene2GO,
              gene2GO = go_map)
          retlist[[ont]][["test_stat"]] <- new(
              "elimScore",
              testStatistic = topGO::GOKSTest,
              name = "KS test",
              cutOff = cutoff)
          retlist[[ont]][["test_result"]] <- topGO::getSigGroups(
                                                        retlist[[ont]][["godata"]],
                                                        retlist[[ont]][["test_stat"]])
        },
        "weight" = {
          retlist[[ont]][["type"]] <- "weight"
          retlist[[ont]][["godata"]] <- new(
              "topGOdata",
              description = ont,
              ontology = ont,
              allGenes = fisher_genes,
              annot = topGO::annFUN.gene2GO,
              gene2GO = go_map)
          retlist[[ont]][["test_stat"]] <- new(
              "weightCount",
              testStatistic = topGO::GOFisherTest,
              name = "Fisher test",
              sigRatio = "ratio")
          retlist[[ont]][["test_result"]] <- topGO::getSigGroups(
                                                        retlist[[ont]][["godata"]],
                                                        retlist[[ont]][["test_stat"]])
        },
        {
          message("I do not know that type.")
        }
    ) ## End of the switch
  }  ## End of for BP/MF/CC
  ##retlist[["pdists"]] <- try(plot_histogram(retlist[["test_result"]]@score, bins = 20))
  class(retlist) <- c("topgo_result", "list")
  return(retlist)
}

#' Print a topgo over representation search.
#'
#' @param x List of the various over/under representation analyses
#'  provided by topGO, the associated plots, and coerced enrichResults.
#' @param ... Other args to match the generic.
#' @export
print.topgo_result <- function(x, ...) {
  bp_entries <- nrow(x[["tables"]][["bp_over_enriched"]])
  mf_entries <- nrow(x[["tables"]][["mf_over_enriched"]])
  cc_entries <- nrow(x[["tables"]][["cc_over_enriched"]])
  summary_string <- glue("TopGO found {bp_entries} BP categories, {mf_entries} MF categories, and \\
{cc_entries} CC categories.")
  message(summary_string)
  enrichplot::dotplot(x[["enrich_results"]][["bp"]])
  return(invisible(x))
}

#' Convert a simple_topgo() result to an enrichResult.
#'
#' Same idea as goseq2enrich.
#'
#' @param retlist result from simple_topgo()
#' @param ontology Ontology subtree to act upon.
#' @param pval Cutoff, hmm I think I need to standardize these.
#' @param organism org name/data.
#' @param column Table column to export.
#' @param padjust_method Use this method for the pvalues for the enrich result.
#' @return enrichResult object ready to pass to things like dotplot.
topgo2enrich <- function(retlist, ontology = "mf", pval = 0.05, organism = NULL,
                         column = "fisher", padjust_method = "BH") {
  result_name <- paste0(column, "_", tolower(ontology))

  term_column <- paste0("Term_", column)
  sig_column <- paste0("Significant_", column)
  annot_column <- paste0("Annotated_", column)
  q_column <- paste0("padj_", column)
  interesting_name <- paste0(tolower(ontology), "_interesting")
  godata <- retlist[["godata"]][[result_name]]
  result_data <- retlist[["results"]][[result_name]]
  interesting <- retlist[["tables"]][[interesting_name]]
  bg_genes <- godata@allGenes
  scores <- interesting[[column]]
  adjusted <- p.adjust(scores)
  input_type <- class(retlist[["input"]])
  if (input_type == "data.frame") {
    sig_genes <- rownames(retlist[["input"]])
  } else {
    sig_genes <- retlist[["input"]]
  }
  mesg("Gather genes per category, this is slow.")
  genes_per_category <- gather_ontology_genes(retlist, ontology = ontology, pval = pval,
                                              column = column, include_all = FALSE)
  ## One of the biggest oddities of enrichResult objects: the scores
  ## are explicitly ratio _strings_, thus 0.05 is bizarrely '5/100' and fails otherwise.
  category_genes <- gsub(pattern = ", ", replacement = "/", x = genes_per_category[["sig"]])
  names(category_genes) <- rownames(genes_per_category)
  keepers <- names(category_genes) %in% rownames(interesting)
  category_genes_kept <- category_genes[keepers]
  interesting[["gene_ids"]] <- as.character(category_genes_kept)
  interesting[["tmp"]] <- length(bg_genes)
  ## FIXME: This is _definitely_ wrong for BgRatio
  representation_df <- data.frame(
      "ID" = rownames(interesting),
      "Description" = interesting[[term_column]],
      ## The following two lines are ridiculous, but required for the enrichplots to work.
      "GeneRatio" = paste0(interesting[[sig_column]], "/", interesting[[annot_column]]),
      "BgRatio" = paste0(interesting[[sig_column]], "/", interesting[["tmp"]]),
      "pvalue" = interesting[[column]],
      "p.adjust" = adjusted,
      "qvalue" = interesting[[q_column]],
      "geneID" = interesting[["gene_ids"]],
      "Count" = interesting[[sig_column]],
      stringsAsFactors = FALSE)
  rownames(representation_df) <- representation_df[["ID"]]
  if (is.null(organism)) {
    organism <- "UNKNOWN"
  }
  ret <- new("enrichResult",
             result = representation_df,
             pvalueCutoff = pval,
             pAdjustMethod = padjust_method,
             qvalueCutoff = pval,
             gene = sig_genes,
             universe = godata@graph@nodes,
             ## universe = extID,
             geneSets = list(up=sig_genes),
             ## geneSets = geneSets,
             organism = organism,
             keytype = "UNKNOWN",
             ontology = ontology,
             readable = FALSE)
  return(ret)
}

#' Make pretty tables out of topGO data
#'
#' The topgo function GenTable is neat, but it needs some simplification to not
#' be obnoxious.
#'
#' @param result Topgo result.
#' @param godata The ontology database.
#' @param limit Pvalue limit defining 'significant'.
#' @param limitby Type of test to perform.
#' @param numchar How many characters to allow in the description?
#' @param orderby Which of the available columns to order the table by?
#' @param ranksof Which of the available columns are used to rank the data?
#' @return prettier tables
#' @seealso [topGO]
#' @export
topgo_tables <- function(results, godata, limit = 0.1, limitby = "fisher",
                         numchar = 300, orderby = "fisher", ranksof = "fisher") {
  ## The following if statement could be replaced by get(limitby)
  ## But I am leaving it as a way to ensure that no shenanigans ensue
  if (limitby == "fisher") {
    mf_siglist <- names(which(results[["fisher_mf"]]@score <= limit))
    bp_siglist <- names(which(results[["fisher_bp"]]@score <= limit))
    cc_siglist <- names(which(results[["fisher_bp"]]@score <= limit))
  } else if (limitby == "KS") {
    mf_siglist <- names(which(results[["ks_mf"]]@score <= limit))
    bp_siglist <- names(which(results[["ks_bp"]]@score <= limit))
    cc_siglist <- names(which(results[["ks_bp"]]@score <= limit))
  } else if (limitby == "EL") {
    mf_siglist <- names(which(results[["el_mf"]]@score <= limit))
    bp_siglist <- names(which(results[["el_bp"]]@score <= limit))
    cc_siglist <- names(which(results[["el_bp"]]@score <= limit))
  } else if (limitby == "weight") {
    mf_siglist <- names(which(results[["weight_mf"]]@score <= limit))
    bp_siglist <- names(which(results[["weight_bp"]]@score <= limit))
    cc_siglist <- names(which(results[["weight_bp"]]@score <= limit))
  } else {
    stop("I can only limit by: fisher, KS, EL, or weight.")
  }
  siglist <- list(
      "mf" = mf_siglist,
      "bp" = bp_siglist,
      "cc" = cc_siglist)
  topnode_list <- list(
      "mf" = length(siglist[["mf"]]),
      "bp" = length(siglist[["bp"]]),
      "cc" = length(siglist[["cc"]]))
  interest_lst <- list()
  all_lst <- list()
  for (ont in c("mf", "bp", "cc")) {
    fisher_name <- glue::glue("fisher_{ont}")
    ks_name <- glue::glue("ks_{ont}")
    el_name <- glue::glue("el_{ont}")
    weight_name <- glue::glue("weight_{ont}")
    if (topnode_list[[ont]] == 0) {
      next
    }
    fisher <- try(topGO::GenTable(godata[[fisher_name]], results[[fisher_name]],
                                  topNodes = topnode_list[[ont]], numChar = numchar))
    if ("try-error" %in% class(fisher)) {
      interest_lst[[ont]] <- data.frame()
      all_lst[[ont]] <- data.frame()
      next
    }
    colnames(fisher) <- paste0(colnames(fisher), "_fisher")
    ks <- try(topGO::GenTable(godata[[ks_name]], results[[ks_name]],
                              topNodes = topnode_list[[ont]], numChar = numchar))
    if ("try-error" %in% class(ks)) {
      interest_lst[[ont]] <- data.frame()
      all_lst[[ont]] <- data.frame()
      next
    }
    colnames(ks) <- paste0(colnames(ks), "_ks")
    el <- try(topGO::GenTable(godata[[el_name]], results[[el_name]],
                              topNodes = topnode_list[[ont]], numChar = numchar))
    if ("try-error" %in% class(el)) {
      interest_lst[[ont]] <- data.frame()
      all_lst[[ont]] <- data.frame()
      next
    }
    colnames(el) <- paste0(colnames(el), "_el")
    weight <- try(topGO::GenTable(godata[[weight_name]], results[[weight_name]],
                                  topNodes = topnode_list[[ont]], numChar = numchar))
    if ("try-error" %in% class(weight)) {
      interest_lst[[ont]] <- data.frame()
      all_lst[[ont]] <- data.frame()
      next
    }
    colnames(weight) <- paste0(colnames(weight), "_weight")
    rownames(fisher) <- fisher[["GO.ID_fisher"]]
    fisher[["GO.ID_fisher"]] <- NULL
    rownames(ks) <- ks[["GO.ID_ks"]]
    ks[["GO.ID_ks"]] <- NULL
    rownames(el) <- el[["GO.ID_el"]]
    el[["GO.ID_el"]] <- NULL
    rownames(weight) <- weight[["GO.ID_weight"]]
    weight[["GO.ID_weight"]] <- NULL
    all <- merge(fisher, ks, by = "row.names", all = TRUE)
    rownames(all) <- all[["Row.names"]]
    all[["Row.names"]] <- NULL
    all <- merge(all, el, by = "row.names", all = TRUE)
    rownames(all) <- all[["Row.names"]]
    all[["Row.names"]] <- NULL
    all <- merge(all, weight, by = "row.names", all = TRUE)
    rownames(all) <- all[["Row.names"]]
    all[["Row.names"]] <- NULL
    colnames(all) <- gsub(x = colnames(all), pattern = "^result1_", replacement = "")
    all[["fisher"]] <- gsub(x = all[["fisher"]], pattern = "^< ", replacement = "")
    all[["fisher"]] <- as.numeric(all[["fisher"]])
    all[["ks"]] <- gsub(x = all[["ks"]], pattern = "^< ", replacement = "")
    all[["ks"]] <- as.numeric(all[["ks"]])
    all[["el"]] <- gsub(x = all[["el"]], pattern = "^< ", replacement = "")
    all[["el"]] <- as.numeric(all[["el"]])
    all[["weight"]] <- gsub(x = all[["weight"]], pattern = "^< ", replacement = "")
    all[["weight"]] <- as.numeric(all[["weight"]])
    all_lst[[ont]] <- all

    if ("data.frame" %in% class(all)) {
      all[["padj_fisher"]] <- p.adjust(all[["fisher"]])
      all[["padj_ks"]] <- p.adjust(all[["ks"]])
      all[["padj_el"]] <- p.adjust(all[["el"]])
      all[["padj_weight"]] <- p.adjust(all[["weight"]])
      interest_idx <- all[[limitby]] <= limit
      na_idx <- is.na(interest_idx)
      interest_idx[na_idx] <- FALSE
      interesting <- all[interest_idx, ]
      interesting[["ont"]] <- ont
      interest_lst[[ont]] <- interesting
    }
  }

  tables <- list(
      "mf_subset" = all_lst[["mf"]],
      "bp_subset" = all_lst[["bp"]],
      "cc_subset" = all_lst[["cc"]],
      "mf_interesting" = interest_lst[["mf"]],
      "bp_interesting" = interest_lst[["bp"]],
      "cc_interesting" = interest_lst[["cc"]])
  return(tables)
}

#' Make a go mapping from IDs in a format suitable for topGO.
#'
#' When using a non-supported organism, one must write out mappings in the
#' format expected by topgo.  This handles that process and gives a summary of
#' the new table.
#'
#' @param goid_map TopGO mapping file.
#' @param go_db If there is no goid_map, create it with this data frame.
#' @param overwrite Rewrite the mapping file?
#' @return Summary of the new goid table.
#' @seealso [topGO]
#' @export
make_id2gomap <- function(goid_map = "reference/go/id2go.map", go_db = NULL,
                          overwrite = FALSE) {
  id2go_test <- file.info(goid_map)
  goids_dir <- dirname(goid_map)
  new_go <- NULL
  if (!file.exists(goids_dir)) {
    dir.create(goids_dir, recursive = TRUE)
  }
  if (isTRUE(overwrite)) {
    if (is.null(go_db)) {
      stop("There is neither a id2go file nor a data frame of goids.")
    } else {
      message("Generating a id2go file in the format expected by topGO.")
      new_go <- reshape2::dcast(go_db, ID ~ ., value.var = "GO",
                                fun.aggregate = paste, collapse = ",")

      write.table(new_go, file = goid_map, sep = "\t",
                  row.names = FALSE, quote = FALSE, col.names = FALSE)
      rm(id2go_test)
    }
  } else {
    ## overwrite is not true
    if (is.na(id2go_test[["size"]])) {
      if (is.null(go_db)) {
        stop("There is neither a id2go file nor a data frame of goids.")
      } else {
        message("Attempting to generate a id2go file in the format expected by topGO.")
        new_go <- reshape2::dcast(go_db, ID~., value.var = "GO",
                                  fun.aggregate = paste, collapse = ",")
        write.table(new_go, file = goid_map, sep = "\t",
                    row.names = FALSE, quote = FALSE, col.names = FALSE)
        id2go_test <- file.info(goid_map)
      }
    } else {
      ## There already exists a file, so return its stats
      ## new_go = id2go_test
      new_go <- goid_map
    }
  }

  ## Pass back an easier to handle go database.
  return(new_go)
}

hpgl_topdiffgenes <- function(scores, df = get0("sig_genes"), direction = "up") {
  ## Testing parameters
  ##scores = pvals
  ##df = epi_cl14clbr_high
  ## Here is the original topDiffGenes
  ## topDiffGenes <- function(allScore) {
  ##   return(allScore < 0.01)
  ##}
  ## my version of this will expect a limma result table from which I will
  ## extract the entries with low p-values and logFCs which are high or low
  quartiles <- summary(df)
}

#' A very simple selector of strong scoring genes (by p-value)
#'
#' This function was provided in the topGO documentation, but not defined.
#' It was copied/pasted here.  I have ideas for including up/down expression
#' but have so far deemed them not needed because I am feeding topGO
#' already explicit lists of genes which are up/down/whatever.
#' But it still is likely to be useful to be able to further subset the data.
#'
#' @param allScore The scores of the genes
#' @export
topDiffGenes <- function(allScore) { return(allScore < 0.01) }

#' Plot the ontology DAG.
#'
#' This function was stolen from topgo in order to figure out where it was failing.
#'
#' @param graph Graph from topGO
#' @return Weights!
#' @export
getEdgeWeights <- function(graph) {
  weightsList <- graph::edgeWeights(graph)
  to <- lapply(weightsList, names)
  from <- graph::nodes(graph)

  if (any(is.na(unlist(to))) || any(is.na(from))) {
    stop("Edge names do not match node names.")
  }

  edge.names <- paste(rep(from, Biobase::listLen(to)), unlist(to), sep = "~")
  edge.weights <- unlist(weightsList)
  names(edge.weights) <- edge.names
  return(edge.weights)
}

#' A minor hack of the topGO GOplot function.
#'
#' This allows me to change the line widths from the default.
#'
#' @param dag DAG tree of ontologies.
#' @param sigNodes Set of significant ontologies (with p-values).
#' @param dag.name Name for the graph.
#' @param edgeTypes Types of the edges for graphviz.
#' @param nodeShape.type Shapes on the tree.
#' @param genNodes Generate the nodes?
#' @param wantedNodes Subset of the ontologies to plot.
#' @param showEdges Show the arrows?
#' @param useFullNames Full names of the ontologies (they can get long).
#' @param oldSigNodes I dunno.
#' @param nodeInfo Hmm.
#' @param maxchars Maximum characters per line inside the shapes.
#' @return Topgo plot!
#' @seealso [topGO]
#' @export
hpgl_GOplot <- function(dag, sigNodes, dag.name = "GO terms", edgeTypes = TRUE,
                        nodeShape.type = c("box", "circle", "ellipse", "plaintext")[3],
                        genNodes = NULL, wantedNodes = NULL, showEdges = TRUE,
                        useFullNames = TRUE, oldSigNodes = NULL, nodeInfo = NULL,
                        maxchars = 30) {
  ## Original function definition had nodeInfo = nodeInfo
  if(!missing(sigNodes)) {
    sigNodeInd <- TRUE
  } else {
    sigNodeInd <- FALSE
  }

  ## we set the global Graphviz attributes
  ## graphAttrs <- getDefaultAttrs(layoutType = 'dot')
  graphAttrs <- Rgraphviz::getDefaultAttrs(layoutType = 'dot')
  graphAttrs[["cluster"]] <- NULL
  graphAttrs[["edge"]][["arrowsize"]] <- "0.4"
  graphAttrs[["edge"]][["weight"]] <- "0.01"

  ##graphAttrs[["graph"]][["splines"]] <- FALSE
  graphAttrs[["graph"]][["size"]] <- "12.0,12.0"
  graphAttrs[["graph"]][["margin"]] <- "0.0,0.0"
  ##  graphAttrs[["graph"]][["ranksep"]] <- "0.02"
  ##  graphAttrs[["graph"]][["nodesep"]] <- "0.30"

  ## set the node shape
  graphAttrs[["node"]][["shape"]] <- nodeShape.type
  ##graphAttrs[["node"]][["fixedsize"]] <- FALSE
  ## set the fontsize for the nodes labels
  graphAttrs[["node"]][["fontsize"]] <- "20.0"
  graphAttrs[["node"]][["height"]] <- "2.0"
  graphAttrs[["node"]][["width"]] <- "3.0"
  graphAttrs[["graph"]][["size"]] <- "12,12"
  graphAttrs[["node"]][["color"]] <- "lightblue"
  graphAttrs[["node"]][["fontname"]] <- "arial"
  graphAttrs[["node"]][["style"]] <- "invis"

  ## set the local attributes lists
  nodeAttrs <- list()
  edgeAttrs <- list()

  ## try to use adaptive node size
  ##nodeAttrs[["fixedsize"]][nodes(dag)] <- rep(FALSE, numNodes(dag))

  if(is.null(nodeInfo)) {
    nodeInfo <- character(graph::numNodes(dag))
    names(nodeInfo) <- graph::nodes(dag)
  } else {
    ##        nodeInfo <- paste('\\\n', nodeInfo, sep = '')
    nodeInfo <- gsub(pattern = "(\\w.{18}).*(\\\\\\n)",
                     replacement = "\\1\\2", x = nodeInfo, perl = TRUE)
    nodeInfo <- glue("\\\n{nodeInfo}")
  }

  ## a good idea is to use xxxxxxx instead of GO:xxxxxxx as node labes
  node.names <- graph::nodes(dag)
  if(!useFullNames) {
    nodeAttrs[["label"]] <- sapply(node.names,
                                   function(x) {
                                     return(paste(substr(x, 4, nchar(node.names[1])),
                                                  nodeInfo[x], sep = ''))
                                   })
  } else {
    nodeAttrs[["label"]] <- glue("{node.names}{nodeInfo}")
    names(nodeAttrs[["label"]]) <- node.names
  }

  ## we will change the shape and the color of the nodes that generated the dag
  if(!is.null(wantedNodes)) {
    diffNodes <- setdiff(wantedNodes, genNodes)
    if(length(diffNodes) > 0) {
      nodeAttrs[["color"]][diffNodes] <- rep("lightblue", .ln <- length(diffNodes))
      nodeAttrs[["shape"]][diffNodes] <- rep("circle", .ln)
      nodeAttrs[["height"]][diffNodes] <- rep("0.45", .ln)
      ## nodeAttrs[["width"]][diffNodes] <- rep("0.6", .ln)
      ## nodeAttrs[["fixedsize"]][wantedNodes] <- rep(TRUE, .ln)
    }
  }

  ## we will change the shape and the color of the nodes we want back
  if(!is.null(genNodes)) {
    nodeAttrs[["color"]][genNodes] <- rep("lightblue", .ln <- length(genNodes))
    nodeAttrs[["shape"]][genNodes] <- rep("box", .ln)
    ## nodeAttrs[["fixedsize"]][genNodes] <- rep(FALSE, .ln)
  }

  ## we will use different fillcolors for the nodes
  if(sigNodeInd) {
    if(!is.null(oldSigNodes)) {
      old.logSigNodes <- log10(sort(oldSigNodes[graph::nodes(dag)]))
      old.range <- range(old.logSigNodes)
      logSigNodes <- log10(sort(sigNodes))
      logSigNodes[logSigNodes < old.range[1]] <- old.range[1]
      logSigNodes[logSigNodes > old.range[2]] <- old.range[2]
      ## debug:  old.range == range(logSigNodes)
      ## if(!identical(all.equal(old.range, range(logSigNodes)), TRUE)) {
      ## stop('some stupid error here :)')
      ## }
    } else {
      old.logSigNodes <- logSigNodes <- log10(sort(sigNodes))
    }

    sigColor <- round(logSigNodes - range(logSigNodes)[1] + 1)
    old.sigColor <- round(old.logSigNodes - range(old.logSigNodes)[1] + 1)

    mm <- max(sigColor, old.sigColor)
    sigColor <- sigColor + (mm - max(sigColor))

    colorMap <- heat.colors(mm)
    nodeAttrs[["fillcolor"]] <- unlist(
        lapply(sigColor, function(x) return(colorMap[x])))
  }

  if(!showEdges) {
    graphAttrs[["edge"]][["color"]] <- "white"
  } else {
    ## if we want to differentiate between 'part-of' and 'is-a' edges
    if (edgeTypes) {
      ## 0 for a is_a relation,  1 for a part_of relation
      ## edgeAttrs[["color"]] <- ifelse(getEdgeWeights(dag) == 0, 'black', 'red')
      edgeAttrs[["color"]] <- ifelse(
          hpgltools::getEdgeWeights(dag) == 0, "black", "black")
    }
  }

  final_dag <- Rgraphviz::agopen(graph = dag, name = dag.name, attrs = graphAttrs,
                                 nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs)
  return(final_dag)
}

#' A hack of topGO's groupDensity()
#'
#' This just adds a couple wrappers to avoid errors in groupDensity.
#'
#' @param object TopGO enrichment object.
#' @param whichGO Individual ontology group to compare against.
#' @param ranks Rank order the set of ontologies?
#' @param rm.one Remove pvalue = 1 groups?
#' @return plot of group densities.
hpgl_GroupDensity <- function(object, whichGO, ranks = TRUE, rm.one = FALSE) {
  ##   groupMembers <- try(topGO::genesInTerm(object, whichGO)[[1]])
  groupMembers <- try(topGO::genesInTerm(object))
  if (class(groupMembers)[1] == "try-error") {
    return(NULL)
  }
  allS <- topGO::geneScore(object, use.names = TRUE)
  if (rm.one) {
    allS <- allS[allS < 0.99]
  }
  xlab <- "Gene' score"
  if (ranks) {
    allS <- BiocGenerics::rank(allS, ties.method = "random")
    xlab <- "Gene's rank"
  }
  group <- as.integer(names(allS) %in% groupMembers)
  xx <- data.frame(
      "score" = allS,
      "group" = factor(
          group, labels = paste0(c('complementary', whichGO), "  (", table(group), ")")))
  plot <- lattice::densityplot(
                       ~ score | group, data = xx, layout = c(1, 2), xlab = xlab)
  return(plot)
}

## EOF
