#' A temporary alias to subset_genes
#'
#' @param ... Parameters passed to subset_genes().
#' @export
exclude_genes <- function(...) {
  message("Note, I renamed this to subset_genes().")
  subset_genes(...)
}

#' Exclude some genes given a pattern match
#'
#' Because I am too lazy to remember that expressionsets use matrix subsets for
#' gene and sample.  Also those methods lead to shenanigans when I want to know
#' what happened to the data over the course of the subset.
#'
#' @param input Expressionset containing exp object.
#' @param column fData column to use for subsetting.
#' @param method Either remove explicit rows, or keep them.
#' @param ids Specific IDs to exclude.
#' @param warning_cutoff Print the sample IDs for anything which has less than this percent left.
#' @param meta_column Save the amount of data lost to this metadata column when not null.
#' @param patterns Character list of patterns to remove/keep
#' @param ... Extra arguments are passed to arglist, currently unused.
#' @return A smaller expt
#' @seealso [create_expt()] [Biobase]
#' @examples
#'  \dontrun{
#'   all_expt <- create_expt(metadata)
#'   ## This assumes a column in the metadata named 'txtype' containing the
#'   ## information telling us what type of transcript each gene is.
#'   no_ribosomes <- exclude_genes_expt(all_expt, column = "txtype",
#'                                      patterns = c("snRNA", "tRNA", "rRNA"))
#'   i_hate_these_genes <- exclude_genes_expt(all_expt, ids = c("gene1", "gene2"))
#'   only_ribosomes <- exclude_genes_expt(all_expt, method = "keep")
#' }
#' @export
subset_genes <- function(input, column = "txtype", method = "remove", ids = NULL,
                         warning_cutoff = 90, meta_column = NULL,
                         patterns = c("snRNA", "tRNA", "rRNA"), ...) {
  message("Generic gene subset function.")
}

#' FIXME: Redo this with S4. Subset a SE bu genes.
#'
#' @param se SummarizedExperiment from which to remove some genes..
#' @param column fData column to use for subsetting.
#' @param method Either remove explicit rows, or keep them.
#' @param ids Specific IDs to exclude.
#' @param warning_cutoff Print the sample IDs for anything which has less than this percent left.
#' @param meta_column Save the amount of data lost to this metadata column when not null.
#' @param patterns Character list of patterns to remove/keep
#' @param ... Extra arguments are passed to arglist, currently unused.#'
subset_genes_se <- function(se, column = "txtype", method = "remove", ids = NULL,
                            warning_cutoff = 90, meta_column = NULL,
                            patterns = c("snRNA", "tRNA", "rRNA"), ...) {
  arglist <- list(...)
  annotations <- rowData(se)
  expression <- assay(se)
  if (is.null(ids) && is.null(annotations[[column]])) {
    message("The ", column, " column is null, doing nothing.")
    return(se)
  }

  orig <- se
  pattern_string <- ""
  for (pat in patterns) {
    pattern_string <- glue("{pattern_string}{pat}|")
  }
  silly_string <- gsub(pattern = "\\|$", replacement = "", x = pattern_string)
  idx <- rep(x = TRUE, times = nrow(annotations))
  if (is.null(ids)) {
    idx <- grepl(pattern = silly_string, x = annotations[[column]], perl = TRUE)
  } else if (is.logical(ids) || is.numeric(ids)) {
    idx <- ids
  } else {
    idx <- rownames(expression) %in% ids
  }
  kept <- NULL
  removed <- NULL
  kept_sums <- NULL
  removed_sums <- NULL
  if (method == "remove") {
    kept <- se[!idx, ]
    removed <- se[idx, ]
  } else {
    kept <- se[idx, ]
    removed <- se[!idx, ]
  }
  tximport_info <- txinfo(se)
  if (!is.null(tximport_info)) {
    df <- tximport_info[["raw"]][["abundance"]][idx, ]
    tximport_info[["raw"]][["abundance"]] <- df
    df <- tximport_info[["raw"]][["counts"]][idx, ]
    tximport_info[["raw"]][["counts"]] <- df
    df <- tximport_info[["raw"]][["length"]][idx, ]
    tximport_info[["raw"]][["length"]] <- df
    df <- tximport_info[["scaled"]][["abundance"]][idx, ]
    tximport_info[["scaled"]][["abundance"]] <- df
    df <- tximport_info[["scaled"]][["counts"]][idx, ]
    tximport_info[["scaled"]][["counts"]] <- df
    df <- tximport_info[["scaled"]][["length"]][idx, ]
    tximport_info[["scaled"]][["length"]] <- df
    txinfo(se) <- tximport_info
  }
  message("subset_genes(), before removal, there were ",
          nrow(rowData(orig)), " genes, now there are ",
          nrow(rowData(kept)), ".")
  all_tables <- assay(orig)
  all_sums <- colSums(all_tables)
  kept_tables <- assay(kept)
  kept_sums <- colSums(kept_tables)
  removed_tables <- assay(removed)
  removed_sums <- colSums(removed_tables)
  pct_kept <- (kept_sums / all_sums) * 100.0
  pct_na <- is.na(pct_kept)
  pct_kept[pct_na] <- 0
  pct_removed <- (removed_sums / all_sums) * 100.0
  pct_na <- is.na(pct_removed)
  pct_removed[pct_na] <- 0
  summary_table <- rbind(kept_sums, removed_sums, all_sums,
                         pct_kept, pct_removed)
  rownames(summary_table) <- c("kept_sums", "removed_sums", "all_sums",
                               "pct_kept", "pct_removed")
  mesg("Average percent of the counts kept after filtering: ",
       toString(sprintf(fmt = "%.3f", mean(pct_kept))))
  ## FIXME: This should be handled by dispatch
  S4Vectors::metadata(se)[["summary_table"]] <- summary_table
  if (!is.null(meta_column)) {
    colData(kept)[meta_column] <- summary_table["pct_removed", ]
  }

  warning_idx <- summary_table["pct_kept", ] < warning_cutoff
  if (sum(warning_idx) > 0) {
    message("There are ", sum(warning_idx), " samples which kept less than ",
            warning_cutoff, " percent counts.")
    print(summary_table["pct_kept", warning_idx])
  }
  return(kept)
}
setGeneric("subset_genes")

setMethod(
  "subset_genes", signature = signature(input = "SummarizedExperiment"),
  definition = function(input, column = "txtype", method = "remove", ids = NULL,
                        warning_cutoff = 90, meta_column = NULL,
                        patterns = c("snRNA", "tRNA", "rRNA"), ...) {
    subset_genes_se(input, column = column, method = method, ids = ids,
                    warning_cutoff = warning_cutoff, meta_column = meta_column,
                    patterns = patterns, ...)
  })

## EOF
