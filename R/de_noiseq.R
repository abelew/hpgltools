## NOISeq looks to me like it does not use a statistical model in the same way as
## limma/deseq/edger, but instead implements its own (un)guided batch correction method
## and expects the user to use that data as input for its set of comparisons.

## Thus, I think setting model_batch to TRUE is likely the place to use their default function
## 'ARSyNseq()' and if we use svaseq, simply send the modified counts to the various noiseq()
## functions which follow.

#' Perform pairwise comparisons using noiseq.
#'
#' @param input Expressionset to compare.
#' @param conditions Set of conditions to query
#' @param batches known batches in the data, or a surrogate estimator.
#' @param model_cond Add condition to the model?
#' @param model_batch Add batch to the model, noiseq has its own combat-like method,
#'  so maybe not necessary?
#' @param annot_df Extra annotations.
#' @param k Taken from the noiseq docs.
#' @param norm Normalization method (noiseq oddly defaults to rpkm).
#' @param factor Metadata factor over which to iterate.
#' @param lc taken from the noiseq docs.
#' @param r taken from the noiseq docs.
#' @param adj taken from the noiseq docs.
#' @param a0per taken from the noiseq docs.
#' @param filter Filter the data?
#' @param keepers Perform the comparison only over these specific contrasts instead of all.
#' @param ... Extra arguments.
#' @return List similar to deseq_pairwise/edger_pairwise/etc.
#' @seealso DOI:10.1093/nar/gkv711
#' @export
noiseq_pairwise <- function(input = NULL, conditions = NULL,
                            batches = NULL, model_cond = TRUE,
                            model_batch = TRUE, model_sv = NULL,
                            model_intercept = FALSE, alt_model = NULL,
                            annot_df = NULL,
                            k = 0.5, norm = "tmm", factor = "condition",
                            lc = 1, r = 20, adj = 1.5, a0per = 0.9, filter = 1,
                            keepers = NULL, keep_underscore = TRUE, ...) {
  arglist <- list(...)

  mesg("Starting noiseq pairwise comparisons.")
  input <- sanitize_expt(input, keep_underscore = keep_underscore)
  input_data <- choose_binom_dataset(input, force = force)
  design <- pData(input)
  conditions <- design[["condition"]]
  conditions_table <- table(conditions)
  batches <- design[["batch"]]
  batches_table <- table(batches)
  data <- input_data[["data"]]
  conditions <- droplevels(as.factor(conditions))
  batches <- droplevels(as.factor(batches))

  noiseq_input <- NOISeq::readData(input_data[["data"]], factors = pData(input))
  norm_input <- NOISeq::ARSyNseq(noiseq_input, factor = "condition", batch = FALSE,
                                 norm = norm, logtransf = FALSE)

  ## Yes I know NOISeq doesn't use models in the same way as other
  ## methods I have applied, but this will make it easier to set up
  ## the contrasts.
  model_choice <- choose_model(norm_input, conditions = conditions, batches = batches,
                               model_batch = model_batch, model_cond = model_cond,
                               model_intercept = model_intercept, model_sv = model_sv,
                               alt_model = alt_model, keep_underscore = keep_underscore,
                               ...)
  model_including <- model_choice[["including"]]
  if (class(model_choice[["model_batch"]])[1] == "matrix") {
    model_batch <- model_choice[["model_batch"]]
  }
  model_data <- model_choice[["chosen_model"]]
  model_string <- model_choice[["chosen_string"]]
  apc <- make_pairwise_contrasts(model_data, conditions, keepers = keepers,
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
    ## Noiseq uses the levels of the factor to define numerator/denominator.
    pData(norm_input)[[factor]] <- as.factor(pData(norm_input)[[factor]])
    pData(norm_input)[[factor]] <- relevel(pData(norm_input)[[factor]], numerator)
    ## Noiseq recasts the pData() as a factor and blows away my levels!
    tmp_file <- tmpmd5file(pattern = "noiseq_density_theta", fileext = ".png")
    this_plot <- png(filename = tmp_file)
    controlled <- dev.control("enable")
    noiseq_table <- sm(NOISeq::noiseqbio(
      norm_input, k = k, norm = norm, factor = factor, lc = lc,
      r = r, adj = adj, plot = TRUE, a0per = a0per, filter = filter,
      conditions = c(numerator, denominator)))
    if (class(noiseq_table)[1] != "try-error") {
      density_theta_plots[[name]] <- grDevices::recordPlot()
    }
    dev.off()
    removed <- file.remove(tmp_file)
    invert <- FALSE
    actual_comparison <- strsplit(noiseq_table@comparison, " - ")[[1]]
    actual_numerator <- actual_comparison[1]
    actual_denominator <- actual_comparison[2]
    if (actual_numerator != numerator) {
      invert <- TRUE
    }
    noiseq_result <- noiseq_table@results[[1]]
    ## It looks to me like noiseq flips the logFC compared to other methods.
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
      "batches_table" = batches_table,
      "coefficients" = coefficient_df,
      "conditions" = conditions,
      "conditions_table" = conditions_table,
      "contrast_list" = contrast_list,
      "contrasts" = apc,
      "contrasts_performed" = apc[["names"]],
      "density_theta_plots" = density_theta_plots,
      "input_data" = input,
      "method" = "noiseq",
      "model" = model_data,
      "model_string" = model_string,
      "norm_input" = norm_input)
  class(retlist) <- c("noiseq_pairwise", "list")
  return(retlist)
}

## EOF
