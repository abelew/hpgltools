#' Run simple_gprofiler on every table from extract_significant_genes()
#'
#' @param sig Result from extract_significant_genes
#' @param according_to Use this result type for the gprofiler searches.
#' @param together Concatenate the up/down genes into one set?
#' @param plot_type Choose a plot method as the default.
#' @param sleep Give the gProfiler servers a break between queries.
#' @param ... Arguments to pass to simple_gprofiler().
#' @export
all_gprofiler <- function(sig, according_to = "deseq", together = FALSE,
                          sleep = 7, plot_type = "dotplot", ...) {
  ret <- list()
  input_up <- list()
  input_down <- list()
  source <- "significant"
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
    slept <- Sys.sleep(sleep)
    name <- sig_names[i]
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
    mesg("Starting ", name, ".")

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
      ret[[retname_up]] <- simple_gprofiler2(up, first_col = fc_col,
                                             plot_type = plot_type,
                                             ...)
      #ret[[retname_up]] <- sm(simple_gprofiler(up, first_col = fc_col))
    } else {
      ret[[retname_up]] <- NULL
    }
    if (down_elements > 0) {
      slept <- Sys.sleep(10)
      ret[[retname_down]] <- simple_gprofiler2(down, first_col = fc_col,
                                               plot_type = plot_type,
                                               ...)
      #ret[[retname_down]] <- sm(simple_gprofiler(down, first_col = fc_col))
    } else {
      ret[[retname_down]] <- NULL
    }
  }
  class(ret) <- "all_gprofiler"
  return(ret)
}

#' Run searches against the web service g:Profiler.
#'
#' This is the beginning of a reimplementation to use gprofiler2.  However,
#' AFAICT gprofiler2 does not yet actually work for anything other than their GO
#' data.
#'
#' @param sig_genes Guess!  The set of differentially expressed/interesting
#'  genes.
#' @param species Organism supported by gprofiler.
#' @param convert Use gProfileR's conversion utility?
#' @param first_col First place used to define the order of 'significant'.
#' @param second_col If that fails, try a second column.
#' @param do_go Perform GO search?
#' @param do_kegg Perform KEGG search?
#' @param do_reactome Perform reactome search?
#' @param do_mi Do miRNA search?
#' @param do_tf Search for transcription factors?
#' @param do_corum Do corum search?
#' @param do_hp Do the hp search?
#' @param do_hpa Do the hpa search?
#' @param do_wp Do the wp search?
#' @param significant Only return the statistically significant hits?
#' @param exclude_iea Passed directly to gprofiler2.
#' @param do_under Perform under-representation search?
#' @param evcodes Get the set of evcodes in the data?  This makes it take
#'  longer.
#' @param threshold p-value 'significance' threshold.
#' @param adjp Method to adjust p-values.
#' @param domain_scope Passed to gprofiler2.
#' @param bg Background genes.
#' @param ordered Is the data in a ranked order by significance?
#' @param id_col Which column in the table should be used for gene ID
#'  crossreferencing?  gProfiler uses Ensembl ids.  So if you have a table of
#'  entrez or whatever, translate it!
#' @param plot_type Use this plot type for images.
#' @param excel Print the results to an excel file?
#' @param enrich_id_column Column from which to extract more readable gene IDs when
#'  creating a clusterProfiler-compatible enrich object.
#' @param ... Primarily for changing options when writing a xlsx output.
#' @return a list of results for go, kegg, reactome, and a few more.
#' @seealso [gProfiler]
#' @examples
#' \dontrun{
#'  gprofiler_is_nice_and_easy <- simple_gprofiler(genes, species='mmusculus')
#' }
#' @export
simple_gprofiler2 <- function(sig_genes, species = "hsapiens", convert = TRUE,
                              first_col = "deseq_logfc", second_col = "logfc", do_mf = TRUE,
                              do_bp = TRUE, do_cc = FALSE, do_kegg = TRUE, do_reactome = TRUE,
                              do_mi = TRUE, do_tf = TRUE, do_corum = TRUE, do_hp = TRUE,
                              do_hpa = TRUE, do_wp = TRUE, significant = TRUE,
                              exclude_iea = FALSE, do_under = FALSE, evcodes = TRUE,
                              threshold = 0.05, adjp = "g_SCS", domain_scope = "annotated",
                              bg = NULL, min_genes = 10, ordered = TRUE, id_col = "row.names",
                              plot_type = "dotplot", excel = NULL, min_go_level = 3, ...) {

  gene_list <- NULL
  num_genes <- 0
  gene_ids <- NULL
  if (is.null(id_col)) {
    id_col <- "row.names"
  }
  if ("data.frame" %in% class(sig_genes)) {
    if (id_col == "row.names") {
      gene_ids <- rownames(sig_genes)
    } else {
      gene_ids <- sig_genes[[id_col]]
    }
  } else {
    if (is.null(names(sig_genes))) {
      gene_ids <- as.character(sig_genes)
    } else {
      gene_ids <- names(sig_genes)
    }
  }
  ## gProfiler is somewhat focused on human data, turn some searches off if
  ## we are not querying human data.
  if (species != "hsapiens") {
    do_hpa <- FALSE
  }

  retlst <- list()
  sources <- c()
  type_names <- c()
  ## I just noticed that you can pass GO:MF GO:BP and GO:CC
  ## As a result I perhaps should change this.
  if (isTRUE(do_bp)) {
    retlst[["BP"]] <- data.frame()
    sources <- c(sources, "GO:BP")
    type_names <- c(type_names, "BP")
  }
  if (isTRUE(do_cc)) {
    retlst[["CC"]] <- data.frame()
    sources <- c(sources, "GO:CC")
    type_names <- c(type_names, "CC")
  }
  if (isTRUE(do_corum)) {
    retlst[["CORUM"]] <- data.frame()
    sources <- c(sources, "CORUM")
    type_names <- c(type_names, "CORUM")
  }
  if (isTRUE(do_hp)) {
    retlst[["HP"]] <- data.frame()
    sources <- c(sources, "HP")
    type_names <- c(type_names, "HP")
  }
  if (isTRUE(do_hpa)) {
    retlst[["HPA"]] <- data.frame()
    sources <- c(sources, "HPA")
    type_names <- c(type_names, "HPA")
  }
  if (isTRUE(do_kegg)) {
    retlst[["KEGG"]] <- data.frame()
    sources <- c(sources, "KEGG")
    type_names <- c(type_names, "KEGG")
  }
  if (isTRUE(do_mi)) {
    retlst[["MIRNA"]] <- data.frame()
    sources <- c(sources, "MIRNA")
    type_names <- c(type_names, "MIRNA")
  }
  if (isTRUE(do_mf)) {
    retlst[["MF"]] <- data.frame()
    sources <- c(sources, "GO:MF")
    type_names <- c(type_names, "MF")
  }
  if (isTRUE(do_reactome)) {
    retlst[["REAC"]] <- data.frame()
    sources <- c(sources, "REAC")
    type_names <- c(type_names, "REAC")
  }
  if (isTRUE(do_tf)) {
    retlst[["TF"]] <- data.frame()
    sources <- c(sources, "TF")
    type_names <- c(type_names, "TF")
  }
  if (isTRUE(do_wp)) {
    retlst[["WP"]] <- data.frame()
    sources <- c(sources, "WP")
    type_names <- c(type_names, "WP")
  }

  if (sum(grepl(pattern = "gene:", x = gene_ids)) > 0) {
    warning("Hey, it looks like you forgot to strip off the htseq prefix for the gene IDs.")
    gene_ids <- gsub(x = gene_ids, pattern = "gene:", replacement = "")
  }

  retlst[["input"]] <- sig_genes
  interactive_plots <- list()
  gost_plots <- list()
  gost_links <- list()
  sig_tables <- list()
  ## Set up a vector to count the number of results observed.
  num_hits <- rep(0, length(type_names))
  names(num_hits) <- type_names
  num_genes <- length(gene_ids)
  if (num_genes <= min_genes) {
    message("There are only, ", num_genes, " returning null.")
    return(NULL)
  }
  for (t in seq_along(type_names)) {
    source <- sources[t]
    type <- type_names[t]
    mesg("Performing gProfiler ", type, " search of ",
         num_genes, " genes against ", species, ".")
    Sys.sleep(3)
    ## To avoid the error: "'names' attribute [14] must be the same length as
    ## the vector [1]"
    gene_ids <- as.vector(gene_ids)
    a_result <- sm(try(gprofiler2::gost(
      query = gene_ids, organism = species, evcodes = evcodes, significant = significant,
      ordered_query = ordered, user_threshold = threshold, correction_method = adjp,
      domain_scope = domain_scope, custom_bg = bg, sources = source)))
    a_df <- data.frame(stringsAsFactors = FALSE)
    if ("try-error" %in% class(a_result)) {
      mesg("The ", type, " method failed for this organism.")
    } else if (is.null(a_result)) {
      mesg("There was no result for the ", type, " search.")
    } else {
      a_df <- a_result[["result"]]
      sig_idx <- a_df[["p_value"]] <= threshold
      sig_df <- a_df[sig_idx, ]
      mesg(type, " search found ", nrow(sig_df), " hits.")
      num_hits[[type]] <- nrow(sig_df)
      sig_tables[[type]] <- sig_df
      gost_links[[type]] <- sm(gprofiler2::gost(
        query = gene_ids, organism = species, evcodes = evcodes, significant = significant,
        ordered_query = ordered, user_threshold = threshold, correction_method = adjp,
        domain_scope = domain_scope, custom_bg = bg, sources = source, as_short_link = TRUE))
      interactive_plots[[type]] <- sm(try(
        gprofiler2::gostplot(a_result, capped = TRUE, interactive = TRUE), silent = TRUE))
      gost_plots[[type]] <- sm(try(
        gprofiler2::gostplot(a_result, capped = FALSE, interactive = FALSE), silent = TRUE))
    }
    enrich_name <- paste0(type, "_enrich")
    ## Note to self, now that I think about it I think gprofiler2 provides its own p-adjustment.
    retlst[[type]] <- a_df
    retlst[[enrich_name]] <- gprofiler2enrich(retlst, ontology = type,
                                              cutoff = threshold,
                                              min_go_level = min_go_level)
    ##    retlst[[enrich_name]] <- gprofiler2enrich(retlst, ontology = type,
    ##                                              cutoff = threshold, enrich_ids = enrich_ids,
    ##                                              min_go_level = min_go_level)
  } ## End iterating over the set of default sources.
  retlst[["num_genes"]] <- num_genes
  retlst[["interactive_plots"]] <- interactive_plots
  retlst[["num_hits"]] <- num_hits
  retlst[["gost_plots"]] <- gost_plots
  retlst[["gost_links"]] <- gost_links
  retlst[["significant"]] <- sig_tables

  if (!is.null(excel)) {
    mesg("Writing data to: ", excel, ".")
    excel_ret <- sm(try(write_gprofiler_data(retlst, excel = excel, ...)))
    retlst[["excel"]] <- excel_ret
    mesg("Finished writing data.")
  }

  if (plot_type == "barplot") {
    retlst[["pvalue_plots"]] <- try(plot_gprofiler2_pval(retlst))
  } else if (plot_type == "dotplot") {
    retlst[["pvalue_plots"]] <- try(plot_gprofiler2_pval(retlst))
    mesg("Add a little logic here to use enrichplot::dotplot().")
  } else {
    retlst[["pvalue_plots"]] <- list()
  }

  retlst[["species"]] <- species
  retlst[["threshold"]] <- threshold
  class(retlst) <- c("gprofiler_result", "list")
  return(retlst)
}

#' Redirect users to simple_gprofiler2
#'
#' @param ... Arguments passed to simple_gprofiler2()
#' @export
simple_gprofiler <- function(...) {
  simple_gprofiler2(...)
}

#' Recast gProfiler data to the output class produced by clusterProfiler.
#'
#' I would like to use the various clusterProfiler plots more easily.
#' Therefore I figured it would be advantageous to coerce the various
#' outputs from gprofiler and friends into the datastructure produced by
#' clusterProfiler.
#'
#' @param retlst Output from simple_gprofiler()
#' @param ontology Category type to extract, currently only GO?
#' @param cutoff Use a p-value cutoff to get only the significant
#'  categories?
#' @param organism Set the orgdb organism name?
#' @param padjust_method what it says on the tin.
#' @return The same 'enrich' datastructure produced by clusterProfiler.
#' @export
gprofiler2enrich <- function(retlst, ontology = "MF", cutoff = 1,
                             organism = NULL, padjust_method = "BH",
                             ## enrich_ids = NULL,
                             min_go_level = 3) {
  is_go <- FALSE
  if (ontology == "MF" || ontology == "BP" || ontology == "CC") {
    is_go <- TRUE
  }
  interesting <- retlst[[ontology]]
  sig_genes <- c()
  sig_genes_input <- retlst[["input"]]
  ##if (!is.null(enrich_ids)) {
  ##  rownames(sig_genes_input) <- make.names(enrich_ids[["enrich"]], unique = TRUE)
  ##}
  if (class(sig_genes_input)[1] == "character") {
    sig_genes <- sig_genes_input
  } else if ("data.frame" %in% class(sig_genes_input)) {
    sig_genes <- rownames(retlst[["input"]])
  } else {
    stop("I do not know this input data type when extracting the input genes.")
  }
  if (is.null(interesting)) {
    return(NULL)
  }
  if (nrow(interesting) == 0) {
    return(NULL)
  }

  bg_genes <- sum(!duplicated(sort(interesting[["term_id"]])))
  interesting[["tmp"]] <- bg_genes
  interesting[["adjusted"]] <- p.adjust(interesting[["p_value"]], method = padjust_method)
  ##if (!is.null(enrich_ids)) {
  ##  for (row in seq_len(nrow(enrich_ids))) {
  ##    from_id <- enrich_ids[row, "id"]
  ##    to_id <- enrich_ids[row, "enrich"]
  ##    interesting[["intersection"]] <- gsub(x = interesting[["intersection"]], pattern = from_id,
  ##                                          replacement = to_id)
  ##  }
  ##}

  genes_per_category <- interesting[, c("term_id", "intersection")]
  category_genes <- gsub(pattern = ",\\s*", replacement = "/",
                         x = genes_per_category[["intersection"]])

  ## Right now the cutoff is 1.0, which is not particularly interesting/useful.
  interesting_cutoff_idx <- interesting[["p_value"]] <= cutoff
  interesting_cutoff <- interesting[interesting_cutoff_idx, ]

  ## Note that for the moment I am repeating the pvalue/p.adjust/qvalue because
  ## I am reasonably certain that gprofiler2 does its own adjustment.
  representation_df <- data.frame(
    "ID" = interesting[["term_id"]],
    "Description" = interesting[["term_name"]],
    ## The following two lines are ridiculous, but required for the enrichplots to work.
    "GeneRatio" = paste0(interesting[["intersection_size"]], "/", interesting[["term_size"]]),
    "BgRatio" = paste0(interesting[["term_size"]], "/", interesting[["query_size"]]),
    "pvalue" = interesting[["p_value"]],
    "p.adjust" = interesting[["p_value"]],
    "qvalue" = interesting[["p_value"]],
    "geneID" = category_genes,
    "Count" = interesting[["intersection_size"]],
    stringsAsFactors = FALSE)
  rownames(representation_df) <- representation_df[["ID"]]
  if (is.null(organism)) {
    organism <- "UNKNOWN"
  }
  ret <- new("enrichResult",
             result = representation_df,
             pvalueCutoff = cutoff,
             pAdjustMethod = padjust_method,
             qvalueCutoff = cutoff,
             gene = sig_genes,
             ## universe = extID,
             geneSets = list(up = sig_genes),
             ## geneSets = geneSets,
             organism = organism,
             keytype = "UNKNOWN",
             ontology = ontology,
             readable = FALSE)
  if (is.numeric(min_go_level) && isTRUE(is_go)) {
    dropped <- seq_len(min_go_level)
    for (drop in dropped) {
      ret <- clusterProfiler::dropGO(ret, level = drop)
    }
  }
  return(ret)
}

## EOF
