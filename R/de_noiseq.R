## NOISeq looks to me like it does not use a statistical model in the same way as
## limma/deseq/edger, but instead implements its own (un)guided batch correction method
## and expects the user to use that data as input for its set of comparisons.

## Thus, I think setting model_batch to TRUE is likely the place to use their default function
## 'ARSyNseq()' and if we use svaseq, simply send the modified counts to the various noiseq()
## functions which follow.

#' Perform pairwise comparisons using noiseq.
#'
#' @param input Expressionset to compare.
#' @param model_fstring Formula string describing the model of interest.
#' @param null_fstring Formula string describing the null model.
#' @param model_svs Surrogate matrix or how to find them.
#' @param extra_contrasts Extra contrasts beyond the default set to seek.
#' @param annot_df Extra annotations.
#' @param force Force the data even if it violates the method's assumptions.
#' @param keepers Perform the comparison only over these specific
#'  contrasts instead of all.
#' @param batch Noiseq batch factor.
#' @param logtransf noiseq log transformer.
#' @param k Taken from the noiseq docs.
#' @param norm Normalization method (noiseq oddly defaults to rpkm).
#' @param factor Metadata factor over which to iterate.
#' @param lc taken from the noiseq docs.
#' @param r taken from the noiseq docs.
#' @param adj taken from the noiseq docs.
#' @param a0per taken from the noiseq docs.
#' @param filter Filter the data?
#' @param keep_underscore Sanitize out underscores?
#' @param ... Extra arguments.
#' @return List similar to deseq_pairwise/edger_pairwise/etc.
#' @seealso DOI:10.1093/nar/gkv711
#' @export
noiseq_pairwise <- function(input = NULL, model_fstring = "~ 0 + condition + batch",
                            null_fstring = "~", model_svs = NULL,
                            extra_contrasts = NULL, annot_df = NULL,
                            force = FALSE, keepers = NULL, batch = FALSE,
                            logtransf = FALSE, k = 0.5, norm = "tmm",
                            factor = NULL, lc = 1, r = 20, adj = 1.5,
                            a0per = 0.9, filter = 1,
                            keep_underscore = TRUE, ...) {
  arglist <- list(...)
  if (isTRUE(filter)) {
    filter <- 1
  }

  mesg("Starting noiseq pairwise comparisons.")
  current <- state(input)
  filteredp <- current[["filter"]]
  if ((is.null(filteredp) || filteredp == "raw") &&
        isTRUE(filter)) {
    input <- sm(normalize(input, filter = filter))
  }
  fctrs <- get_formula_factors(model_fstring)
  factors <- fctrs[["factors"]]
  condition_column <- factors[1]
  if (is.null(factor)) {
    factor <- condition_column
  }
  input <- sanitize_se(input, keep_underscore = keep_underscore, factors = factors)
  input_data <- choose_binom_dataset(input, force = force)
  count_mtrx <- input_data[["data"]]
  design <- pData(input)
  conditions <- droplevels(as.factor(design[[condition_column]]))
  batches <- droplevels(as.factor(design[["batch"]]))
  condition_table <- table(conditions)
  batch_table <- table(batches)
  condition_levels <- levels(conditions)
  design[[condition_column]] <- conditions
  if (length(levels(batches)) > 1) {
    design[["batch"]] <- batches
  } else {
    design[["batch"]] <- NULL
    batch <- FALSE
  }
  mesg("This noiseq pairwise comparison should compare across:")
  print(condition_table)
  noiseq_input <- NOISeq::readData(input_data[["data"]], factors = design)
  if (is.null(design[["batch"]])) {
    norm_input <- try(NOISeq::ARSyNseq(noiseq_input, factor = NULL, batch = TRUE,
                                       norm = norm, logtransf = logtransf), silent = TRUE)
  } else {
    norm_input <- try(NOISeq::ARSyNseq(noiseq_input, factor = "batch", batch = TRUE,
                                       norm = norm, logtransf = FALSE), silent = TRUE)
  }
  if ("try-error" %in% class(norm_input)) {
    norm_input <- noiseq_input
  }
  model_mtrx <- model.matrix(as.formula(model_fstring), data = design)
  apc <- make_pairwise_contrasts(model_mtrx, conditions, keepers = keepers,
                                 keep_underscore = keep_underscore)
  contrast_list <- list()
  result_list <- list()
  lrt_list <- list()
  sc <- vector("list", length(apc[["names"]]))
  end <- length(apc[["names"]])
  coefficient_df <- data.frame()
  final_coef_colnames <- c()
  density_theta_plots <- list()
  for (con in seq_along(apc[["names"]])) {
    name <- apc[["names"]][[con]]
    numerator <- apc[["numerators"]][[name]]
    denominator <- apc[["denominators"]][[name]]
    tmp_file <- tmpmd5file(pattern = "noiseq_density_theta", fileext = ".png")
    this_plot <- png(filename = tmp_file)
    controlled <- dev.control("enable")
    noiseq_table <- sm(NOISeq::noiseqbio(
      norm_input, k = k, norm = norm, factor = condition_column, lc = lc,
      r = r, adj = adj, plot = TRUE, a0per = a0per,
      conditions = c(numerator, denominator)))
    if (class(noiseq_table)[1] != "try-error") {
      density_theta_plots[[name]] <- grDevices::recordPlot()
    }
    plotted <- dev.off()
    removed <- file.remove(tmp_file)
    invert <- FALSE
    actual_comparison <- strsplit(noiseq_table@comparison, " - ")[[1]]
    actual_numerator <- actual_comparison[1]
    actual_denominator <- actual_comparison[2]
    if (actual_numerator != numerator) {
      invert <- TRUE
    }
    noiseq_result <- noiseq_table@results[[1]]
    rename_col <- colnames(noiseq_result) == "log2FC"
    colnames(noiseq_result)[rename_col] <- "logFC"
    noiseq_result[["p"]] <- 1.0 - noiseq_result[["prob"]]
    noiseq_result[["adjp"]] <- p.adjust(noiseq_result[["p"]])
    coefficient_names <- colnames(noiseq_result)[c(1, 2)]
    current_colnames <- colnames(coefficient_df)
    tmp_coef_df <- noiseq_result[, c(1, 2)]
    kept <- ! coefficient_names %in% current_colnames
    final_coef_colnames <- c(final_coef_colnames, coefficient_names[kept])
    if (ncol(coefficient_df) == 0) {
      coefficient_df <- data.frame(row.names = rownames(tmp_coef_df))
    }
    if (sum(kept) > 0) {
      for (co in colnames(tmp_coef_df)) {
        coefficient_df <- cbind.data.frame(coefficient_df, log2(tmp_coef_df[[co]]))
      }
    }
    colnames(noiseq_result) <- c("num_mean", "den_mean", "theta", "prob", "logFC", "p", "adjp")
    if (isTRUE(invert)) {
      colnames(noiseq_result) <- c("den_mean", "num_mean", "theta", "prob", "logFC", "p", "adjp")
      noiseq_result <- noiseq_result[, c("num_mean", "den_mean", "theta", "prob", "logFC", "p", "adjp")]
      noiseq_result[["logFC"]] <- -1.0 * noiseq_result[["logFC"]]
    }
    result_list[[name]] <- noiseq_result
  }
  final_coef_colnames <- gsub(x = final_coef_colnames, pattern = "_mean", replacement = "")
  colnames(coefficient_df) <- final_coef_colnames
  coef_na_idx <- is.na(coefficient_df)
  coefficient_df[coef_na_idx] <- 0

  retlist <- list(
      "all_tables" = result_list,
      "batches" = batches,
      "batch_table" = batch_table,
      "coefficients" = coefficient_df,
      "conditions" = conditions,
      "condition_table" = condition_table,
      "contrast_list" = contrast_list,
      "contrasts" = apc,
      "contrasts_performed" = apc[["names"]],
      "density_theta_plots" = density_theta_plots,
      "input_data" = input,
      "method" = "noiseq",
      "model" = model_mtrx,
      "model_fstring" = model_fstring,
      "norm_input" = norm_input)
  class(retlist) <- c("noiseq_pairwise", "list")
  return(retlist)
}

#' Print the results from noiseq_pairwise().
#'
#' @param x result to print.
#' @param ... arbitrary arguments.
#' @export
print.noiseq_pairwise <- function(x, ...) {
  summary_string <- glue("The results from the Noiseq pairwise analysis.")
  message(summary_string)
  return(invisible(x))
}

#' Passes noiseq result to write_de_table.
#'
#' @param data noiseq_pairwise() result.
#' @param ... arbitrary arguments.
#' @export
write_noiseq <- function(data, ...) {
  result <- write_de_table(data, type = "noiseq", ...)
  return(result)
}

## EOF
