## model_varpartition.r: Use VariancePartition to analyze the measureability of
## variance in a data set given different experimental models.  Variance
## Partition is fun, but weird.  This file seeks to simplify and standardize
## these methods.

#' A shortcut for replotting the percent plots from variancePartition.
#'
#' In case I wish to look at different numbers of genes from variancePartition
#' and/or different columns to sort from.
#'
#' @param varpart_output List returned by varpart()
#' @param n How many genes to plot.
#' @param column The df column to use for sorting.
#' @param decreasing high->low or vice versa?
#' @return The percent variance bar plots from variancePartition!
#' @seealso [variancePartition]
#' @export
replot_varpart_percent <- function(varpart_output, n = 30, column = NULL, decreasing = TRUE) {
  sorted <- varpart_output[["sorted_df"]]
  if (!is.null(column)) {
    if (column %in% colnames(sorted)) {
      sorted <- sorted[order(sorted[[column]], decreasing = decreasing), ]
    } else {
      message("The column ", column,
              "is not in the sorted data frame returned by varpart(). ",
              "Leaving the data frame alone.")
    }
  }
  new_plot <- variancePartition::plotPercentBars(sorted[1:n, ])
  retlist <- list(
      "resorted" = sorted,
      "plot" = new_plot)
  class(retlist) <- "reordered_varpart"
  return(retlist)
}

#' Print the result of a reordered variance partition analysis.
#'
#' @param x List of a resorted variance partition analysis and its plot.
#' @param ... Other args to match the generic.
#' @export
print.reordered_varpart <- function(x, ...) {
  plot(x[["plot"]])
  return(invisible(x))
}

#' Use variancePartition to try and understand where the variance lies in a data set.
#'
#' The arguments and usage of variancePartition are a bit opaque.  This function
#' attempts to fill in reasonable values and simplify its invocation.
#'
#' @param input Some data
#' @param fstring Formula string describing the factors to query.
#' @param do_fit Perform a fitting using variancePartition?
#' @param cor_gene Provide a set of genes to look at the correlations, defaults
#'  to the first gene.
#' @param cpus Number cpus to use
#' @param genes Number of genes to count.
#' @param parallel Use doParallel?
#' @param strict_filter Perform a strict filtering of the results via median_by_factor and dropping
#'  any genes with a 0.
#' @param modify_input Add annotation columns with the variance/factor?
#' @return List of plots and variance data frames
#' @seealso [variancePartition] DOI:10.1186/s12859-016-1323-z.
#' @export
simple_varpart <- function(input, fstring = "~ condition + batch",
                           do_fit = FALSE, cor_gene = 1,
                           cpus = NULL, genes = 40, parallel = TRUE,
                           strict_filter = TRUE, modify_input = TRUE) {
  cl <- NULL
  para <- NULL
  lib_result <- sm(requireNamespace("variancePartition"))
  att_result <- sm(try(attachNamespace("variancePartition"), silent = TRUE))
  lib_result <- sm(requireNamespace("BiocParallel"))
  att_result <- sm(try(attachNamespace("BiocParallel"), silent = TRUE))
  if (isTRUE(parallel)) {
    cl <- NULL
    if (is.null(cpus)) {
      cpus <- parallel::detectCores() - 4
      if (cpus < 1) {
        cpus <- 1
      }
      cl <- parallel::makeCluster(cpus)
    } else {
      cl <- parallel::makeCluster(cpus)
    }
    para <- doParallel::registerDoParallel(cl)
    ## multi <- BiocParallel::MulticoreParam()
  }
  design <- colData(input)
  rank_test <- test_design_model_rank(design, fstring)
  fctrs <- rank_test[["factors"]]
  condition_fctr <- rank_test[["factors"]][1]
  for (f in rank_test[["factors"]]) {
    design[[f]] <- droplevels(as.factor(design[[f]]))
  }

  test_formula <- as.formula(fstring)
  ## I think the simple filter is insufficient and I need there to be
  ## no genes with 0 counts in any one condition.
  norm <- sm(normalize(input, filter = "simple"))
  if (isTRUE(strict_filter)) {
    test <- sm(median_by_factor(norm, fact = condition_fctr, fun = "mean"))
    all_condition_gt_zero_idx <- rowSums(test[["medians"]] == 0) == 0
    kept_gt <- rownames(assay(norm))[all_condition_gt_zero_idx]
    norm <- norm[kept_gt, ]
  }
  data <- assay(norm)

  mesg("Fitting the expressionset to the model, this is slow.")
  fit_extract <- try(variancePartition::fitExtractVarPartModel(data, test_formula, design))
  ## my_extract <- try(variancePartition::fitVarPartModel(data, my_model, design))
  if ("try-error" %in% class(fit_extract)) {
    mesg("A couple of common errors:
An error like 'vtv downdated' may be because there are too many 0s, filter the data and rerun.
An error like 'number of levels of each grouping factor must be < number of observations' means
that the factor used is not appropriate for the analysis - it really only works for factors
which are shared among multiple samples.")
    message("Retrying with only condition in the model.")
    test_formula <- as.formula(glue("~ {condition_fctr}"))
    fit_extract <- try(variancePartition::fitExtractVarPartModel(data, test_formula, design))
    if ("try-error" %in% class(fit_extract)) {
      message("Attempting again with only condition failed.")
      stop()
    }
  }

  ## A new dataset has some NAs!
  na_idx <- is.na(fit_extract)
  if (sum(na_idx) > 0) {
    warning("There are ", sum(na_idx), " NAs in this data, something may be wrong.")
    message("There are ", sum(na_idx), " NAs in this data, something may be wrong.")
    message("Converting NAs to 0.")
    fit_extract[na_idx] <- 0
  }

  sorted_fit <- sortCols(fit_extract)
  order_idx <- order(sorted_fit[[condition_fctr]], decreasing = TRUE)
  sorted_fit <- sorted_fit[order_idx, ]
  ## Recent error noticed when checking that variances sum to 1
  ## This is because sometimes we have smaller data sets
  if (genes > ncol(sorted_fit)) {
    genes <- ncol(sorted_fit)
  }

  percent_plot <- variancePartition::plotPercentBars(sorted_fit[1:genes, ])
  partition_plot <- variancePartition::plotVarPart(sorted_fit)

  fitting <- NULL
  stratify_batch_plot <- NULL
  stratify_condition_plot <- NULL
  if (isTRUE(do_fit)) {
    ## Try fitting with lmer4
    fitting <- variancePartition::fitVarPartModel(exprObj = data,
                                                  formula = test_formula, data = design)
    last_fact <- fctrs[length(fctrs)]
    idx <- order(design[[condition_fctr]], design[[last_fact]])
    ##first <- variancePartition::plotCorrStructure(fitting, reorder = idx)
    test_strat <- data.frame(Expression = data[3, ],
                             first = design[[condition_fctr]],
                             last = design[[last_fact]])
    batch_expression <- as.formula(glue("Expression ~ first"))
    cond_expression <- as.formula(glue("Expression ~ last"))
    stratify_batch_plot <- variancePartition::plotStratify(batch_expression, test_strat)
    stratify_condition_plot <- variancePartition::plotStratify(cond_expression, test_strat)
  }

  if (isTRUE(parallel)) {
    para <- parallel::stopCluster(cl)
  }

  ret <- list(
    "fstring" = fstring,
    "model_used" = test_formula,
    "percent_plot" = percent_plot,
    "partition_plot" = partition_plot,
    "sorted_df" = sorted_fit,
    "fitted_df" = fit_extract,
    "fitting" = fitting,
    "stratify_batch_plot" = stratify_batch_plot,
    "stratify_condition_plot" = stratify_condition_plot)
  if (isTRUE(modify_input) && nrow(rowData(input)) > 0) {
    new_input <- input
    tmp_annot <- rowData(new_input)
    tmp_annot[["Row.names"]] <- NULL
    added_data <- sorted_fit
    colnames(added_data) <- glue("variance_{colnames(added_data)}")
    ## Note that we are getting these variance numbers from data which was filtered
    ## thus we need all.x to get the IDs to match up.
    tmp_annot <- merge(tmp_annot, added_data, by = "row.names", all.x = TRUE)
    rownames(tmp_annot) <- tmp_annot[["Row.names"]]
    tmp_annot[["Row.names"]] <- NULL
    annot_order <- rownames(assay(new_input))
    tmp_annot <- tmp_annot[annot_order, ]
    ## Make it possible to use a generic expressionset, though maybe this is
    ## impossible for this function.
    rowData(new_input) <- tmp_annot
    ret[["modified_input"]] <- new_input
  }
  class(ret) <- "varpart"
  return(ret)
}

#' Print variance partition results.
#'
#' @param x List of results from variancePartition including the model
#'  information, percent/partition plots, dataframes of the
#'  fitted/sorted data by variance, etc.
#' @param ... Other args to match the generic.
#' @export
print.varpart <- function(x, ...) {
  summary_string <- glue("The result of using variancePartition with the model:
{x[['model_string']]}")
  message(summary_string)
  plot(x[["partition_plot"]])
  return(invisible(x))
}

#' Attempt to use variancePartition's fitVarPartModel() function.
#'
#' Note the word 'attempt'.  This function is so ungodly slow that it probably
#' will never be used.
#'
#' @param exp Input expressionset.
#' @param factors Set of factors to query
#' @param cpus Number of cpus to use in doParallel.
#' @return Summaries of the new model,  in theory this would be a nicely
#'  batch-corrected data set.
#' @seealso [variancePartition]
varpart_summaries <- function(exp, factors = c("condition", "batch"), cpus = 6) {
  cl <- parallel::makeCluster(cpus)
  doParallel::registerDoParallel(cl)
  model_string <- "~ "
  for (fact in factors) {
    model_string <- glue("{model_string} (1|{fact}) + ")
  }
  model_string <- gsub(pattern = "\\+ $", replacement = "", x = model_string)
  my_model <- as.formula(model_string)
  norm <- sm(normalize_exp(exp, filter = TRUE))
  data <- assay(norm)
  design <- colData(exp)
  summaries <- variancePartition::fitVarPartModel(data, my_model, design, fxn = summary)
  return(summaries)
}

## EOF
