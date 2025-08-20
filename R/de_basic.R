## de_basic.r: An implementation of a simplified, statistical model unaware
## differential expression method.  This is intended essentially as a negative
## control to get a sense of how 'intrusive' other methods need to be in order
## to get their various results when performing a differential expression
## analysis.

#' The simplest possible differential expression method.
#'
#' Perform a pairwise comparison among conditions which takes
#' nothing into account.  It _only_ takes the conditions, a mean value/variance
#' among them, divides by condition, and returns the result.  No fancy
#' nomalizations, no statistical models, no nothing.  It should be the very
#' worst method possible. But, it should also provide a baseline to compare the
#' other tools against, they should all do better than this, always.
#'
#' Tested in test_27de_basic.R
#' This function was written after the corresponding functions in de_deseq.R,
#' de_edger.R, and de_limma.R.  Like those, it performs the full set of pairwise
#' comparisons and returns a list of the results.  As mentioned above, unlike
#' those, it is purposefully stupid.
#'
#' @param input Count table by sample.
#' @param model_fstring Formula string which describes the experimental model.
#' @param null_fstring Formula string describing the null hypothesis (not used).
#' @param model_svs Method to extract surrogate variables (not used).
#' @param annot_df Extra annotation dataframe.
#' @param keepers Set of specific contrasts to perform instead of all.
#' @param fx What function to use for mean/median?
#' @param keep_underscore Sanitize model underscores?
#' @param ... Extra options passed to arglist.
#' @return Df of pseudo-logFC, p-values, numerators, and denominators.
#' @seealso [deseq_pairwise()] [limma_pairwise()] [edger_pairwise()] [ebseq_pairwise()]
#' @examples
#' \dontrun{
#'  expt <- create_expt(metadata = "sample_sheet.xlsx", gene_info = "annotations")
#'  basic_de <- basic_pairwise(expt)
#'  basic_tables <- combine_de_tables(basic_de)
#' }
#' @export
basic_pairwise <- function(input = NULL, model_fstring = "~ 0 + condition + batch",
                           null_fstring = "~", model_svs = NULL,
                           annot_df = NULL, keepers = NULL,
                           fx = "mean", keep_underscore = TRUE, ...) {
  arglist <- list(...)
  mesg("Starting basic pairwise comparisons.")
  input <- sanitize_expt(input, keep_underscore = keep_underscore)
  fctrs <- get_formula_factors(model_fstring)
  contrast_factor <- fctrs[["factors"]][1]
  input_data <- choose_basic_dataset(input, force = force)
  design <- pData(input)

  conditions <- droplevels(as.factor(design[[contrast_factor]]))
  batches <- droplevels(as.factor(design[["batch"]]))
  data <- as.matrix(input_data[["data"]])
  conditions_table <- table(conditions)
  batches_table <- table(batches)
  condition_levels <- levels(conditions)
  model_mtrx <- model.matrix(as.formula(model_fstring), data = design)
  num_conds <- length(condition_levels)
  ## These will be filled with num_conds columns and numRows(input) rows.
  median_table <- data.frame()
  variance_table <- data.frame()
  ## First use conditions to rbind a table of medians by condition.
  mesg("Basic step 1/3: Creating ", fx, " and variance tables.")
  median_colnames <- c()
  for (d in seq_len(num_conds)) {
    condition_name <- condition_levels[d]
    median_colnames <- append(median_colnames, condition_name)
    columns <- which(conditions == condition_name)
    if (length(columns) == 1) {
      med <- data.frame(data[, columns], stringsAsFactors = FALSE)
      var <- as.data.frame(matrix(NA, ncol = 1, nrow = nrow(med)),
                           stringsAsFactors = FALSE)
    } else {
      med_input <- as.matrix(data[, columns])
      if (fx == "mean") {
        med <- data.frame(matrixStats::rowMeans2(x = med_input, na.rm = TRUE))
      } else {
        med <- data.frame(Biobase::rowMedians(as.matrix(med_input)))
      }
      colnames(med) <- c(condition_name)
      var <- as.data.frame(genefilter::rowVars(as.matrix(med_input)))
      colnames(var) <- c(condition_name)
    }
    if (d == 1) {
      median_table <- med
      variance_table <- var
    } else {
      median_table <- cbind(median_table, med)
      variance_table <- cbind(variance_table, var)
    }
  } ## end creation of median and variance tables.
  colnames(median_table) <- median_colnames
  colnames(variance_table) <- median_colnames
  rownames(median_table) <- rownames(data)
  rownames(variance_table) <- rownames(data)
  ## We have tables of the median values by condition
  ## Now perform the pairwise comparisons
  comparisons <- data.frame()
  tvalues <- data.frame()
  pvalues <- data.frame()
  num_done <- 0
  column_list <- c()
  total_contrasts <- length(levels(as.factor(conditions)))
  if (is.null(keepers)) {
    total_contrasts <- (total_contrasts * (total_contrasts + 1)) / 2
  } else if ("list" %in% class(keepers)) {
    total_contrasts <- length(keepers)
  } else {
    total_contrasts <- (total_contrasts * (total_contrasts + 1)) / 2
  }
  mesg("Basic step 2/3: Performing ", total_contrasts, " comparisons.")
  apc <- make_pairwise_contrasts(model_mtrx, conditions, contrast_factor = contrast_factor,
                                 do_identities = FALSE, do_extras = FALSE,
                                 keepers = keepers, keep_underscore = keep_underscore,
                                 ...)
  contrasts_performed <- c()
  for (d in seq_along(apc[["names"]])) {
    num_done <- num_done + 1
    name  <- apc[["names"]][[d]]
    c_name <- gsub(pattern = "^(.*)_vs_(.*)$", replacement = "\\1", x = name)
    d_name <- gsub(pattern = "^(.*)_vs_(.*)$", replacement = "\\2", x = name)
    contrasts_performed <- append(name, contrasts_performed)
    if (! c_name %in% colnames(median_table)) {
      message("The contrast ", name, " is not in the results.")
      message("If this is not an extra contrast, then this is an error.")
      next
    }
    division <- data.frame(
        median_table[, c_name] - median_table[, d_name])
    column_list <- append(column_list, name)
    colnames(division) <- name
    ## Lets see if I can make a dirty p-value
    xcols <- which(conditions == c_name)
    ycols <- which(conditions == d_name)
    xdata <- as.data.frame(data[, xcols])
    ydata <- as.data.frame(data[, ycols])

    t_data <- vector("list", nrow(xdata))
    p_data <- vector("list", nrow(xdata))
    for (j in seq_len(nrow(xdata))) {
      test_result <- try(suppressWarnings(
        stats::wilcox.test(x = as.numeric(xdata[j, ]), y = as.numeric(ydata[j, ]),
                           exact = FALSE, alternative = "two.sided")),
        silent = TRUE)
      if (class(test_result) == "htest") {
        t_data[[j]] <- test_result[[1]]
        p_data[[j]] <- test_result[[3]]
      } else {
        t_data[[j]] <- 0
        p_data[[j]] <- 1
      }
    } ## Done calculating cheapo p-values

    if (d == 1) {
      comparisons <- division
      tvalues <- t_data
      pvalues <- p_data
    } else {
      comparisons <- cbind(comparisons, division)
      tvalues <- cbind(tvalues, t_data)
      pvalues <- cbind(pvalues, p_data)
    }
  } ## End for each contrast

  ## Because of the way I made tvalues/pvalues into a list
  ## If only 1 comparison was performed, the resulting data structure never gets coerced into a
  ## data frame.  Therefore I am performing this check which, if a single comparison was done, adds
  ## a second column, performs the coercion, then strips it away.  This is a stupid way
  ## of doing what I want.
  if (num_done == 1) {
    tvalues <- cbind(tvalues, t_data)
    pvalues <- cbind(pvalues, p_data)
    tvalues <- as.data.frame(tvalues)
    pvalues <- as.data.frame(pvalues)
    tvalues <- tvalues[-1]
    pvalues <- pvalues[-1]
  }
  comparisons[is.na(comparisons)] <- 0
  tvalues[is.na(tvalues)] <- 0
  pvalues[is.na(pvalues)] <- 1
  rownames(comparisons) <- rownames(data)
  rownames(tvalues) <- rownames(data)
  rownames(pvalues) <- rownames(data)
  all_tables <- list()

  mesg("Basic step 3/3: Creating faux DE Tables.")
  for (e in seq_along(colnames(comparisons))) {
    colname <- colnames(comparisons)[[e]]
    fc_column <- comparisons[, e]
    t_column <- as.numeric(tvalues[, e])
    p_column <- as.numeric(pvalues[, e])
    fc_column[mapply(is.infinite, fc_column)] <- 0
    numer_denom <- strsplit(x = colname, split = "_vs_")[[1]]
    numerator <- numer_denom[1]
    denominator <- numer_denom[2]
    num_col <- paste0("numerator_", fx)
    den_col <- paste0("denominator_", fx)
    fc_table <- data.frame(
        "numerator_var" = variance_table[[numerator]],
        "denominator_var" = variance_table[[denominator]],
        "t" = t_column,
        "p" = p_column,
        "logFC" = fc_column)
    fc_table[[num_col]] <- median_table[[numerator]]
    fc_table[[den_col]] <- median_table[[denominator]]
    fc_table <- fc_table[, c(num_col, den_col, "numerator_var",
                             "denominator_var", "t", "p", "logFC")]
    fc_table[["adjp"]] <- stats::p.adjust(as.numeric(fc_table[["p"]]), method = "BH")

    fc_table[[num_col]] <- signif(
        x = fc_table[[num_col]], digits = 4)
    fc_table[[den_col]] <- signif(
        x = fc_table[[den_col]], digits = 4)
    ## I am thinking to change my mind about this formatting, since
    ## it recasts the numbers as characters, and that is dumb.
    fc_table[["t"]] <- signif(x = fc_table[["t"]], digits = 4)
    fc_table[["logFC"]] <- signif(x = fc_table[["logFC"]], digits = 4)
    rownames(fc_table) <- rownames(data)
    all_tables[[e]] <- fc_table
  }
  names(all_tables) <- colnames(comparisons)

  retlist <- list(
      "all_pairwise" = comparisons,
      "all_tables" = all_tables,
      "conditions_table" = table(conditions),
      "conditions" = conditions,
      "contrasts_performed" = contrasts_performed,
      "input_data" = input,
      "medians" = median_table,
      "method" = "basic",
      "num_contrasts" = total_contrasts,
      "variances" = variance_table)
  class(retlist) <- c("basic_pairwise", "list")
  if (!is.null(arglist[["basic_excel"]])) {
    retlist[["basic_excel"]] <- write_basic(retlist, excel = arglist[["basic_excel"]])
  }
  return(retlist)
}

#' @export
print.basic_pairwise <- function(x, ...) {
  summary_string <- glue("The results from the Basic pairwise analysis.")
  message(summary_string)
  return(invisible(x))
}

#' Attempt to ensure that input data to basic_pairwise() is suitable.
#'
#' basic_pairwise() assumes log2 data as input, use this to ensure that is true.
#'
#' @param input An expressionset containing expt to test and/or modify.
#' @param force If we want to try out other distributed data sets, force it in using me.
#' @param ... future options, I think currently unused.
#' @return data ready for basic_pairwise()
#' @seealso [Biobase] [choose_dataset()] [normalize()]
#' @examples
#' \dontrun{
#'  ready <- choose_basic_dataset(expt)
#' }
choose_basic_dataset <- function(input, force = FALSE, ...) {
  ## arglist <- list(...)
  conditions <- conditions(input)
  batches <- batches(input)
  data <- as.data.frame(exprs(input))
  state <- state(input)
  tran_state <- state[["transform"]]
  libsize <- NULL
  if (is.null(tran_state)) {
    tran_state <- "raw"
  }
  conv_state <- state[["conversion"]]
  ## Note that voom takes care of this for us.
  if (is.null(conv_state)) {
    conv_state <- "raw"
  }
  norm_state <- state[["normalization"]]
  if (is.null(norm_state)) {
    norm_state <- "raw"
  }
  filt_state <- state[["filter"]]
  if (is.null(filt_state)) {
    filt_state <- "raw"
  }

  ready <- input
  if (isTRUE(force)) {
    message("Leaving the data alone, regardless of normalization state.")
  } else {
    if (filt_state == "raw") {
      message("Basic step 0/3: Filtering data.")
      ready <- sm(normalize(ready, filter = TRUE))
    }
    if (norm_state == "raw") {
      message("Basic step 0/3: Normalizing data.")
      ready <- sm(normalize(ready, norm = "quant"))
    }
    if (conv_state == "raw") {
      message("Basic step 0/3: Converting data.")
      ready <- sm(normalize(ready, convert = "cbcbcpm"))
    }
  }
  ## No matter what we do, it must be logged.
  message("I think this is failing? ", class(ready))
  if (tran_state == "raw") {
    message("Basic step 0/3: Transforming data.")
    ready <- normalize(ready, transform = "log2")
  }
  data <- exprs(ready)
  libsize <- colSums(data)
  rm(ready)
  retlist <- list(
      "libsize" = libsize,
      "conditions" = conditions,
      "batches" = batches,
      "data" = data)
  return(retlist)
}

#' Writes out the results of a basic search using write_de_table()
#'
#' Looking to provide a single interface for writing tables from basic and friends.
#'
#' Tested in test_26basic.R
#'
#' @param data Output from basic_pairwise()
#' @param ... Options for writing the xlsx file.
#' @seealso [basic_pairwise()] [write_de_table()]
#' @examples
#' \dontrun{
#'  finished_comparison <- basic_pairwise(expressionset)
#'  data_list <- write_basic(finished_comparison)
#' }
#' @export
write_basic <- function(data, ...) {
  result <- write_de_table(data, type = "basic", ...)
  return(result)
}

## EOF
