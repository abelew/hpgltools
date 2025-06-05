## de_ebseq.r: Simplify the invocation of EBSeq, the Bayesian statistical method
## of performing differential expression analysis.  EBSeq provides a third
## family of DE method: DESeq/EdgeR provides explicitly negative binomial
## distribution aware methods, while limma assumes microarrayish.  Both these
## families use various statistical models.  In contrast, EBSeq uses prior
## probabilities and sampling to gather its values for each gene.

#' Set up model matrices contrasts and do pairwise comparisons of all conditions
#' using EBSeq.
#'
#' Invoking EBSeq is confusing, this should help.
#'
#' @param input Dataframe/vector or expt class containing data, normalization
#'  state, etc.
#' @param patterns Set of expression patterns to query.
#' @param model_fstring Formula string describing the model of interest.
#' @param null_fstring Formula string describing the null model.
#' @param model_svs Matrix of SVs or character describing how to find them.
#' @param keepers Perform a specific set of contrasts instead of all?
#' @param ng_vector I think this is for isoform quantification, but am not yet
#'  certain.
#' @param rounds Number of iterations for doing the multi-test
#' @param target_fdr Definition of 'significant'
#' @param method The default ebseq methodology is to create the set of all
#'  possible 'patterns' in the data; for data sets which are more than
#'  trivially complex, this is not tenable, so this defaults to subsetting the
#'  data into pairs of conditions.
#' @param norm Normalization method to use.
#' @param force Force ebseq to accept bad data (notably NA containing
#'  stuff from proteomics.
#' @param keep_underscore Sanitize away underscores?
#' @param ... Extra arguments currently unused.
#' @return List containing tables from ebseq, the conditions tested, and the
#'  ebseq table of conditions.
#' @seealso [limma_pairwise()] [deseq_pairwise()] [edger_pairwise()] [basic_pairwise()]
#' @examples
#'  \dontrun{
#'   expt <- create_expt(metadata = "sample_sheet.xlsx", gene_info = annotations)
#'   ebseq_de <- ebseq_pairwise(input = expt)
#' }
#' @export
ebseq_pairwise <- function(input = NULL, patterns = NULL,
                           model_fstring = "~ 0 + condition + batch",
                           null_fstring = "~", model_svs = NULL,
                           keepers = NULL, ng_vector = NULL, rounds = 20,
                           target_fdr = 0.05, method = "pairwise_subset",
                           norm = "median", force = FALSE, keep_underscore = TRUE,
                           ...) {
  arglist <- list(...)

  mesg("Starting EBSeq pairwise comparisons.")
  current <- state(input)
  filteredp <- current[["filter"]]
  if ((is.null(filteredp) || filteredp == "raw") &&
        isTRUE(filter)) {
    input <- sm(normalize(input, filter = filter))
  }
  fctrs <- get_formula_factors(model_fstring)
  factors <- fctrs[["factors"]]
  condition_column <- factors[1]
  input <- sanitize_expt(input, keep_underscore = keep_underscore, factors = factors)
  input_data <- choose_binom_dataset(input, force = force)
  design <- pData(input)
  conditions <- droplevels(as.factor(design[[condition_column]]))
  data <- as.matrix(input_data[["data"]])
  condition_table <- table(conditions)
  condition_levels <- levels(conditions)
  numerators <- denominators <- contrasts_performed <- c()

  if (method == "pairwise_subset") {
    result <- ebseq_pairwise_subset(input, model_fstring = model_fstring,
                                    ng_vector = ng_vector, rounds = rounds,
                                    target_fdr = target_fdr, norm = norm, force = force,
                                    keepers = keepers, keep_underscore = keep_underscore,
                                    ...)
    numerators <- result[["numerators"]]
    denominators <- result[["denominators"]]
    contrasts_performed <- result[["contrasts_performed"]]
  } else {
    mesg("Starting single EBSeq invocation.")
    multi <- FALSE
    if (length(condition_levels) < 2) {
      stop("You have fewer than 2 conditions.")
    } else if (length(condition_levels) == 2) {
      mesg("Invoking ebseq with 2-condition parameters.")
      result <- ebseq_two(input, ng_vector = ng_vector, rounds = rounds,
                          target_fdr = target_fdr, norm = norm)
      numerators <- result[["numerator"]]
      denominators <- result[["denominators"]]
    } else if (length(condition_levels) > 5) {
      stop("Beyond 5 conditions generates too many patterns, ",
           "please provide a pattern matrix, or 'all_same'.")
    } else {
      mesg("Invoking ebseq with parameters suitable for a few conditions.")
      result <- ebseq_few(input, model_fstring = model_fstring,
                          conditions, patterns = patterns,
                          ng_vector = ng_vector, rounds = rounds,
                          target_fdr = target_fdr, norm = norm)
    }
  }

  retlist <- list(
      "all_tables" = result,
      "conditions" = conditions,
      "conditions_table" = condition_table,
      "contrasts_performed" = contrasts_performed,
      "denominators" = denominators,
      "method" = "ebseq",
      "numerators" = numerators)
  class(retlist) <- c("ebseq_pairwise", "list")
  return(retlist)
}

#' @export
print.ebseq_pairwise <- function(x, ...) {
  summary_string <- glue("The results from the EBSeq pairwise analysis.")
  message(summary_string)
  return(invisible(x))
}

#' Perform pairwise comparisons with ebseq, one at a time.
#'
#' This uses the same logic as in the various *_pairwise functions to invoke
#' the 'normal' ebseq pairwise comparison for each pair of conditions in an
#' expressionset.  It therefore avoids the strange logic inherent in the ebseq
#' multitest function.
#'
#' @param input Expressionset/expt to perform de upon.
#' @param model_fstring Formula string describing the model of interest.
#' @param ng_vector Passed on to ebseq, I forget what this does.
#' @param rounds Passed on to ebseq, I think it defines how many iterations to
#'  perform before return the de estimates
#' @param target_fdr If we reach this fdr before iterating rounds
#'  times, return.
#' @param keepers Specify a set of contrasts to perform here.#'
#' @param conditions Factor of conditions in the data, used to define the
#'  contrasts.
#' @param norm EBseq normalization method to apply to the data.
#' @param force Flag used to force inappropriate data into the various methods.
#' @param keep_underscore Sanitize away underscores?
#' @param ... Extra arguments passed downstream.
#' @return A pairwise comparison of the various conditions in the data.
#' @seealso [ebseq_pairwise()]
ebseq_pairwise_subset <- function(input, model_fstring = "~ 0 + condition + batch",
                                  ng_vector = NULL, rounds = 10, target_fdr = 0.05,
                                  keepers = NULL, conditions = NULL, norm = "median",
                                  force = FALSE, keep_underscore = TRUE,
                                  ...) {
  mesg("Starting EBSeq pairwise subset.")
  ## Now that I understand pData a bit more, I should probably remove the
  ## conditions/batches slots from my expt classes.
  design <- pData(input)
  fctrs <- get_formula_factors(model_fstring)
  condition_column <- fctrs[["factors"]][1]
  design <- pData(input)
  conditions <- droplevels(as.factor(design[[condition_column]]))
  data <- exprs(input)
  condition_table <- table(conditions)
  condition_levels <- levels(conditions)

  model_mtrx <- model.matrix(as.formula(model_fstring), data = design)
  apc <- make_pairwise_contrasts(model_mtrx, conditions, do_identities = FALSE,
                                 do_extras = FALSE, keepers = keepers,
                                 keep_underscore = keep_underscore,
                                 ...)
  contrasts_performed <- c()
  numerators <- denominators <- c()
  retlst <- list()
  for (c in seq_along(apc[["names"]])) {
    name  <- apc[["names"]][[c]]
    contrasts_performed <- c(name, contrasts_performed)
    a_name <- gsub(pattern = "^(.*)_vs_(.*)$", replacement = "\\1", x = name)
    b_name <- gsub(pattern = "^(.*)_vs_(.*)$", replacement = "\\2", x = name)
    if (! a_name %in% conditions) {
      message("The contrast ", a_name, " is not in the results: ", conditions)
      message("If this is not an extra contrast, then this is an error.")
      next
    }
    idx <- conditions == b_name | conditions == a_name
    pair <- input[, idx]
    pair_data <- exprs(pair)
    pair_conditions <- droplevels(conditions[idx])
    a_result <- ebseq_two(pair_data, pair_conditions, numerator = b_name, denominator = a_name,
                          ng_vector = ng_vector, rounds = rounds, target_fdr = target_fdr,
                          norm = norm, force = force)
    numerators <- c(numerators, b_name)
    denominators <- c(denominators, a_name)
    retlst[[name]] <- a_result
  }
  retlst[["numerators"]] <- numerators
  retlst[["denominators"]] <- denominators
  retlst[["contrasts_performed"]] <- contrasts_performed
  return(retlst)
}

#' Choose the ebseq normalization method to apply to the data.
#'
#' EBSeq provides three normaliation methods.  Median, Quantile, and Rank.
#' Choose among them here.
#'
#' @param data_mtrx This is exprs(expressionset)
#' @param norm The method to pass along.
#' @return a new matrix using the ebseq specific method of choice.
#' @seealso [EBSeq]
ebseq_size_factors <- function(data_mtrx, norm = NULL) {
  ## Set up a null normalization vector
  normalized <- rep(x = 1, times = ncol(data_mtrx))
  names(normalized) <- colnames(data_mtrx)
  ## If the parameter passed matches an ebseq function, use it, if it doesn't,
  ## we still have the null.
  if (norm == "median") {
    normalized <- EBSeq::MedianNorm(data_mtrx, alternative = TRUE)
  } else if (norm == "quantile") {
    normalized <- EBSeq::QuantileNorm(data_mtrx)
  } else if (norm == "rank") {
    normalized <- EBSeq::RankNorm(data_mtrx)
  } else {
    message("I do not know the norm method: ", norm, ", using median.")
    normalized <- EBSeq::MedianNorm(data_mtrx, alternative = TRUE)
  }
  return(normalized)
}

#' Invoke EBMultiTest() when we do not have too many conditions to deal with.
#'
#' Starting at approximately 5 conditions, ebseq becomes too unwieldy to use
#' effectively. But, its results until then are pretty neat.
#'
#' @param data Expressionset/matrix
#' @param conditions Factor of conditions in the data to compare.
#' @param model_fstring Formula string describing the model of interest.
#' @param patterns Set of patterns as described in the ebseq documentation to query.
#' @param ng_vector Passed along to ebmultitest().
#' @param rounds Passed to ebseq.
#' @param target_fdr Passed to ebseq.
#' @param norm Normalization method to apply to the data.
#' @seealso [ebseq_pairwise()]
ebseq_few <- function(data, conditions, model_fstring = "~ 0 + condition + batch",
                      patterns = NULL, ng_vector = NULL, rounds = 10,
                      target_fdr = 0.05, norm = "median") {

  ## Reminder about the meanings of 'patterns':
  ## Each row is a set of mean expression levels
  ## Each column is a condition
  ## Therefore, if a row has all three columns with a '1', then this pattern
  ## signifies a scenario in which all conditions have the same mean
  ## expression.
  ## If a row is '1', '1', '2', '2' then there are two sets of the same
  ## expression, 1 and 2.
  ## etc etc.
  if (is.null(patterns)) {
    patterns <- EBSeq::GetPatterns(conditions)
  } else if (patterns == "all_same") {
    patterns <- data.frame(row.names = "Pattern1")
    for (i in conditions) {
      patterns[1, i] <- 1
    }
  }

  normalized <- ebseq_size_factors(data, norm)
  eb_output <- EBSeq::EBMultiTest(
                          Data = data, NgVector = ng_vector,
                          Conditions = conditions, AllParti = patterns,
                          sizeFactors = normalized, maxround = rounds)
  posteriors <- EBSeq::GetMultiPP(eb_output)
  fold_changes <- EBSeq::GetMultiFC(eb_output)

  pp_df <- as.data.frame(eb_output[["PPMat"]])
  ## Drop the pattern which looks for all the same
  interesting_patterns <- as.data.frame(patterns[-1, ])

  table_lst <- list()
  ## FIXME: Change this to use a portion of make_pairwise_contrasts()
  ## so that numerators and denominators do not get flipped.
  for (i in seq_len(ncol(fold_changes[["FCMat"]]))) {
    column <- colnames(fold_changes[["FCMat"]])[i]
    contrast <- gsub(pattern = "Over", replacement = "_vs_", x = column)
    numerator <- gsub(pattern = "^(.*)_vs_(.*)$", replacement = "\\1", x = contrast)
    denominator <- gsub(pattern = "^(.*)_vs_(.*)$", replacement = "\\2", x = contrast)
    chosen_pattern_idx <- interesting_patterns[[numerator]] ==
      interesting_patterns[[denominator]]
    ## This should provide the name of the pattern which, when I query the
    ## PPMatrix will provide the likelihood that the two are the same
    ## I can therefore take the 1-that to get a likelihood different.
    chosen_pattern <- rownames(interesting_patterns)[chosen_pattern_idx]

    table <- data.frame(row.names = rownames(fold_changes[["FCMat"]]))
    table[["ebseq_FC"]] <- fold_changes[["FCMat"]][, i]
    table[["logFC"]] <- fold_changes[["Log2FCMat"]][, i]
    table[["ebseq_postfc"]] <- fold_changes[["Log2PostFCMat"]][, i]
    table[["ebseq_mean"]] <- fold_changes[["CondMeans"]][, i]
    table[["pp_all_same"]] <- pp_df[["Pattern1"]]
    table[["pp_pair_diff"]] <- pp_df[[chosen_pattern]]
    table[["putative_p_value"]] <- 1 - table[["pp_pair_diff"]]
    table_lst[[contrast]] <- table
  }

  retlst <- list(
      "all_tables" = table_lst,
      "conditions" = conditions,
      "denominator" = denominator,
      "method" = "ebseq",
      "numerator" = numerator)
  return(retlst)
}

#' The primary function used in my EBSeq implementation.
#'
#' Most of the time, my invocation of ebseq will fall into this function.
#'
#' @param pair_data Matrix containing the samples comprising two experimental
#'  factors of interest.
#' @param conditions Factor of conditions in the data.
#' @param numerator Which factor has the numerator in the data.
#' @param denominator Which factor has the denominator in the data.
#' @param fast The EBSeq fast argument.
#' @param ng_vector Passed to ebseq.
#' @param rounds Passed to ebseq.
#' @param Alpha The ebseq alpha parameter.
#' @param Beta The ebseq beta parameter.
#' @param Qtrm Ibid.
#' @param QtrmCut Ibid.
#' @param step1 Ibid.
#' @param step2 Ibid.
#' @param thre Ibid.
#' @param sthre Ibid.
#' @param filter Ibid.
#' @param stopthre Ibid.
#' @param target_fdr Passed to ebseq.
#' @param norm Normalization method of ebseq to apply.
#' @param force Force inappropriate data into ebseq?
#' @return EBSeq result table with some extra formatting.
#' @seealso [ebseq_pairwise()]
ebseq_two <- function(pair_data, conditions,
                      numerator = 2, denominator = 1, fast = TRUE,
                      ng_vector = NULL, rounds = 20, Alpha = NULL,
                      Beta = NULL, Qtrm = 1, QtrmCut = 0, step1 = 1e-06,
                      step2 = 0.01, thre = log(2), sthre = 0,
                      filter = 10, stopthre = 1e-04,
                      target_fdr = 0.05, norm = "median",
                      force = FALSE) {
  normalized <- ebseq_size_factors(pair_data, norm = norm)
  ## I think this should be removed in lieu of the imputation functions
  if (isTRUE(force)) {
    mesg("Forcing out NA values by putting in the mean of all data.")
    ## Put NA values (proteomics) to the mean of the existing values in the hopes
    ## they will not mess anything up too badly.
    na_idx <- is.na(pair_data)
    pair_data[na_idx] <- mean(pair_data, na.rm = TRUE)
  }
  eb_output <- sm(EBSeq::EBTest(
    Data = pair_data, NgVector = NULL, Conditions = conditions,
    sizeFactors = normalized, maxround = rounds, fast = fast,
    Alpha = Alpha, Beta = Beta, Qtrm = Qtrm, QtrmCut = QtrmCut, step1 = step1,
    step2 = step2, thre = thre, sthre = sthre, filter = filter, stopthre = stopthre))
  posteriors <- EBSeq::GetPPMat(eb_output)
  fold_changes <- EBSeq::PostFC(eb_output)
  eb_result <- EBSeq::GetDEResults(eb_output, FDR = target_fdr)
  mean_df <- as.data.frame(eb_output[["Mean"]])
  colnames(mean_df) <- c("numerator_mean", "denominator_mean")
  meanlist_df <- as.data.frame(eb_output[["MeanList"]])
  colnames(meanlist_df) <- c("ebseq_mean")
  varlist_df <- as.data.frame(eb_output[["VarList"]])
  colnames(varlist_df) <- c("ebseq_var")
  p_df <- as.data.frame(eb_result[["PPMat"]])
  table <- data.frame(row.names = rownames(posteriors))
  table[["ebseq_FC"]] <- fold_changes[["RealFC"]]
  table[["logFC"]] <- log2(table[["ebseq_FC"]])
  table[["ebseq_c1mean"]] <- as.numeric(mean_df[["numerator_mean"]])
  table[["ebseq_c2mean"]] <- as.numeric(mean_df[["denominator_mean"]])
  table <- merge(table, meanlist_df, by = "row.names", all.x = TRUE)
  rownames(table) <- table[["Row.names"]]
  table[["Row.names"]] <- NULL
  table <- merge(table, varlist_df, by = "row.names", all.x = TRUE)
  rownames(table) <- table[["Row.names"]]
  table[["Row.names"]] <- NULL
  table[["ebseq_postfc"]] <- fold_changes[["PostFC"]]
  table <- merge(table, p_df, by = "row.names", all.x = TRUE)
  rownames(table) <- table[["Row.names"]]
  table[["Row.names"]] <- NULL
  ## This is incorrect I think, but being used as a placeholder until I figure out how to
  ## properly adjust a set prior probabilities.
  table[["ebseq_adjp"]] <- table[["PPEE"]]
  table[["PPEE"]] <- NULL

  ## Finally, make sure the 'direction' matches my conception of numerator/denominator.
  eb_direction <- fold_changes[["Direction"]]
  eb_numerator <- gsub(pattern = "^(.*) Over (.*)$", replacement = "\\1", x = eb_direction)
  eb_denominator <- gsub(pattern = "^(.*) Over (.*)$", replacement = "\\2", x = eb_direction)
  if (! eb_numerator == denominator) {
    table[["ebseq_FC"]] <- 1 / table[["ebseq_FC"]]
    table[["logFC"]] <- -1 * table[["logFC"]]
    table[["ebseq_postfc"]] <- 1 / table[["ebseq_postfc"]]
  }
  table[["ebseq_FC"]] <- NULL
  return(table)
}

## EOF
