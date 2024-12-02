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
                            model_batch = TRUE, annot_df = NULL,
                            k = 0.5, norm = "rpkm", factor = "condition",
                            lc = 1, r = 20, adj = 1.5, a0per = 0.9, filter = 1,
                            keepers = NULL, keep_underscore = FALSE, ...) {
  arglist <- list(...)

  message("Starting noiseq pairwise comparisons.")
  input <- sanitize_expt(input, keep_underscore = keep_underscore)
  input_data <- choose_binom_dataset(input, force = force)
  design <- pData(input)
  conditions <- design[["condition"]]
  conditions_table <- table(conditions)
  batches <- design[["batch"]]
  batches_table <- table(batches)
  data <- input_data[["data"]]
  conditions <- as.factor(conditions)
  batches <- as.factor(batches)

  noiseq_input <- NOISeq::readData(input_data[["data"]], factors = pData(input))
  norm_input <- NOISeq::ARSyNseq(noiseq_input, factor = "condition", batch = FALSE,
                                 norm = norm, logtransf = FALSE)

  ## Yes I know NOISeq doesn't use models in the same way as other
  ## methods I have applied, but this will make it easier to set up
  ## the contrasts.
  model_choice <- choose_model(norm_input, conditions = conditions,
                               batches = batches,
                               model_batch = model_batch,
                               model_cond = model_cond)
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
  for (con in seq_along(apc[["names"]])) {
    name <- apc[["names"]][[con]]
    numerator <- apc[["numerators"]][[name]]
    denominator <- apc[["denominators"]][[name]]
    ## Noiseq uses the levels of the factor to define numerator/denominator.
    pData(norm_input)[[factor]] <- as.factor(pData(norm_input)[[factor]])
    relevel(pData(norm_input)[[factor]], numerator)
    ## Noiseq recasts the pData() as a factor and blows away my levels!
    noiseq_table <- sm(NOISeq::noiseqbio(
      norm_input, k = 0.5, norm = norm, factor = factor, lc = 1,
      r = 20, adj = 1.5, plot = TRUE, a0per = 0.9, filter = 1,
      conditions = c(numerator, denominator)))
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
    ## It looks to me like noiseq flips the logFC compared to other methods.
    noiseq_result[["p"]] <- 1.0 - noiseq_result[["prob"]]
    if (isTRUE(invert)) {
      noiseq_result[["logFC"]] <- -1.0 * noiseq_result[["logFC"]]
    }
    noiseq_result[["adjp"]] <- p.adjust(noiseq_result[["p"]])
    result_list[[name]] <- noiseq_result
  }

  retlist <- list(
      "all_tables" = result_list,
      "batches" = batches,
      "batches_table" = batches_table,
      "conditions" = conditions,
      "conditions_table" = conditions_table,
      "contrast_list" = contrast_list,
      "contrasts" = apc,
      "contrasts_performed" = apc[["names"]],
      "input_data" = input,
      "method" = "noiseq",
      "model" = model_data,
      "model_string" = model_string)
  class(retlist) <- c("noiseq_pairwise", "list")
  return(retlist)
}

## EOF
