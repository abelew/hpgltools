## Use variancePartition's dream method with limma.

#' Set up a model matrix and set of contrasts for pairwise comparisons using
#' voom/limma.
#'
#' Creates the set of all possible contrasts and performs them using voom/limma.
#'
#' @param input Dataframe/vector or expt class containing count tables,
#'  normalization state, etc.
#' @param conditions Factor of conditions in the experiment.
#' @param batches Factor of batches in the experiment.
#' @param model_cond Include condition in the model?
#' @param model_batch Include batch in the model?  If this is a character
#'  instead of a logical, then it is passed to all_adjusers() to attempt to find
#'  model parameters which describe surrogate variables in the data.
#' @param model_intercept Perform a cell-means or intercept model? A little more
#'  difficult for me to understand.  I have tested and get the same answer
#'  either way.
#' @param alt_model Specify another model.
#' @param extra_contrasts Some extra contrasts to add to the list.
#'  This can be pretty neat, lets say one has conditions A,B,C,D,E
#'  and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like: "c_vs_b_ctrla = (C-B)-A, e_vs_d_ctrla = (E-D)-A,
#'  de_vs_cb = (E-D)-(C-B),"
#' @param annot_df Data frame for annotations.
#' @param libsize I've recently figured out that libsize is far more important
#'  than I previously realized.  Play with it here.
#' @param limma_method Choose one of limma's lm methods.
#' @param limma_robust Make the significance estimation robust?
#' @param voom_norm Use this method to normalize the voom inputs.
#' @param limma_trend Add trend lines to limma's voom plot?
#' @param force Force data which may not be appropriate for limma into it?
#' @param keepers Perform an explicit set of contrasts instead of all.
#' @param ... Use the elipsis parameter to feed options to write_limma().
#' @return List including the following information:
#'  macb = the mashing together of condition/batch so you can look at it
#'  macb_model = The result of calling model.matrix(~0 + macb)
#'  macb_fit = The result of calling lmFit(data, macb_model)
#'  voom_result = The result from voom()
#'  voom_design = The design from voom (redundant from voom_result, but convenient)
#'  macb_table = A table of the number of times each condition/batch pairing happens
#'  cond_table = A table of the number of times each condition appears (the
#'   denominator for the identities)
#'  batch_table = How many times each batch appears
#'  identities = The list of strings defining each condition by itself
#'  all_pairwise = The list of strings defining all the pairwise contrasts
#'  contrast_string = The string making up the makeContrasts() call
#'  pairwise_fits = The result from calling contrasts.fit()
#'  pairwise_comparisons = The result from eBayes()
#'  limma_result = The result from calling write_limma()
#' @seealso [limma] [Biobase] [deseq_pairwise()] [edger_pairwise()] [basic_pairwise()]
#'  DOI:10.1101/2023.03.17.533005
#' @examples
#' \dontrun{
#'  pretend <- dream_pairwise(expt)
#' }
#' @export
dream_pairwise <- function(input = NULL, conditions = NULL,
                           batches = NULL, model_cond = TRUE,
                           model_batch = TRUE, model_intercept = FALSE,
                           alt_model = NULL, extra_contrasts = NULL,
                           annot_df = NULL, libsize = NULL,
                           limma_method = "ls", limma_robust = FALSE, voom_norm = "quantile",
                           limma_trend = FALSE, force = FALSE, keepers = NULL, ...) {
  arglist <- list(...)
  ## This is used in the invocation of a voom() implementation for normalization.
  ## This is for the eBayes() call.

  message("Starting limma/varpart pairwise comparison.")
  san_input <- sanitize_expt(input)
  input_data <- choose_limma_dataset(san_input, force = force)
  design <- pData(san_input)
  if (is.null(conditions)) {
    conditions <- design[["condition"]]
  }
  if (is.null(batches)) {
    batches <- design[["batch"]]
  }

  ## The following small piece of logic is intended to handle situations where we use
  ## tximport for limma (kallisto/sailfish/salmon).
  if (is.null(san_input[["tximport"]])) {
    data <- input_data[["data"]]
  } else {
    data <- edgeR::DGEList(san_input[["tximport"]][["scaled"]][["counts"]])
    data <- edgeR::calcNormFactors(data)
  }

  if (is.null(libsize)) {
    message("libsize was not specified, this parameter has profound effects on limma's result.")
    if (!is.null(san_input[["best_libsize"]])) {
      message("Using the libsize from expt$best_libsize.")
      libsize <- san_input[["best_libsize"]]
    } else if (!is.null(input[["libsize"]])) {
      message("Using the libsize from expt$libsize.")
      libsize <- san_input[["libsize"]]
    } else if (!is.null(
      san_input[["normalized"]][["intermediate_counts"]][["normalization"]][["libsize"]])) {
      libsize <- colSums(data, na.rm = TRUE)
    } else {
      message("Using the libsize from expt$normalized$intermediate_counts$normalization$libsize")
      libsize <- san_input[["normalized"]][["intermediate_counts"]][["normalization"]][["libsize"]]
    }
  } else {
    message("libsize was specified.  This parameter has profound effects on limma's result.")
  }

  if (is.null(libsize)) {
    libsize <- colSums(data, na.rm = TRUE)
  }
  condition_table <- table(conditions)
  batch_table <- table(batches)
  conditions <- as.factor(conditions)
  batches <- as.factor(batches)

  message("Limma step 1/6: choosing model.")
  model <- choose_model(input = san_input,
                        conditions = conditions,
                        batches = batches,
                        model_batch = model_batch,
                        model_cond = model_cond,
                        model_intercept = model_intercept,
                        alt_model = alt_model,
                        ...)
  ##model <- choose_model(input, conditions, batches,
  ##                      model_batch = model_batch,
  ##                      model_cond = model_cond,
  ##                      model_intercept = model_intercept,
  ##                      alt_model = alt_model)
  chosen_model <- model[["chosen_model"]]
  model_string <- model[["chosen_string"]]
  my_formula <- as.formula(model_string)

  voom_plot <- NULL
  message("Attempting voomWithDreamWeights.")
  fun_voom <- variancePartition::voomWithDreamWeights(
    counts = data, formula = model_string,
    data = design, plot = TRUE)
  voom_plot <- grDevices::recordPlot()

  one_replicate <- FALSE
  fun_design <- pData(san_input)
  if (is.null(fun_voom)) {
    ## Apparently voom returns null where there is only 1 replicate.
    message("voom returned null, I am not sure what will happen.")
    one_replicate <- TRUE
    fun_voom <- data
  }

  ## Do the lmFit() using this model
  pairwise_fits <- NULL
  identity_fits <- NULL
  message("Limma/varpart step 3/6: running dream.")
  contrasts <- make_pairwise_contrasts(model = chosen_model, conditions = conditions,
                                       extra_contrasts = extra_contrasts, keepers = keepers)
  all_pairwise_contrasts <- contrasts[["all_pairwise"]]
  all_pairwise_contrasts <- gsub(x = all_pairwise_contrasts, pattern = ",$", replacement = "")
  contrast_vector <- c()
  for (f in all_pairwise_contrasts) {
    name_val <- strsplit(f, "=")[[1]]
    vec_name <- name_val[1]
    vec <- name_val[2]
    num_den <- strsplit(vec, "-")[[1]]
    num <- paste0("condition", num_den[1])
    den <- paste0("condition", num_den[2])
    expression <- paste0(num, " - ", den)
    names(expression) <- vec_name
    contrast_vector <- c(contrast_vector, expression)
  }

  varpart_contrasts <- variancePartition::makeContrastsDream(my_formula, design,
                                                contrasts = contrast_vector)
  fitted_data <- variancePartition::dream(
    exprObj = fun_voom, formula = model_string, data = design, L = varpart_contrasts)

  all_tables <- NULL
  all_pairwise_comparisons <- limma::eBayes(fitted_data,
                                            robust = limma_robust,
                                            trend = limma_trend)
  message("Limma step 6/6: Writing limma outputs.")
  pairwise_results <- make_limma_tables(fit = all_pairwise_comparisons, adjust = "BH",
                                        n = 0, coef = NULL, annot_df = NULL)
  contrasts_performed <- names(pairwise_results)

  retlist <- list(
    "all_pairwise" = all_pairwise,
    "all_tables" = pairwise_results,
    "batches" = batches,
    "batches_table" = batch_table,
    "conditions" = conditions,
    "conditions_table" = condition_table,
    "contrast_string" = contrast_vector,
    "contrasts_performed" = contrasts_performed,
    "dispersion_plot" = voom_plot,
    "fit" = fitted_data,
    "input_data" = input,
    "method" = "limma",
    "model" = model,
    "model_string" = model_string,
    "pairwise_comparisons" = all_pairwise_comparisons,
    "single_table" = all_tables,
    "voom_design" = fun_design,
    "voom_result" = fun_voom)
  class(retlist) <- c("dream_pairwise", "list")
  if (!is.null(arglist[["limma_excel"]])) {
    retlist[["dream_excel"]] <- write_limma(retlist, excel = arglist[["limma_excel"]])
  }
  return(retlist)
}

#' Writes out the results of a limma search using toptable().
#'
#' However, this will do a couple of things to make one's life easier:
#' 1.  Make a list of the output, one element for each comparison of the contrast matrix
#' 2.  Write out the toptable() output in separate .csv files and/or sheets in excel
#' 3.  Since I have been using qvalues a lot for other stuff, add a column for them.
#'
#' @param fit Result from lmFit()/eBayes()
#' @param adjust Pvalue adjustment chosen.
#' @param n Number of entries to report, 0 says do them all.
#' @param coef Which coefficients/contrasts to report, NULL says do them all.
#' @param annot_df Optional data frame including annotation information to
#'  include with the tables.
#' @param intercept Intercept model?
#' @return List of data frames comprising the toptable output for each
#'  coefficient, I also added a qvalue entry to these toptable() outputs.
#' @seealso [limma] [write_xlsx()]
#' @examples
#' \dontrun{
#'  finished_comparison = eBayes(limma_output)
#'  table = make_limma_tables(finished_comparison, adjust = "fdr")
#' }
make_limma_tables <- function(fit = NULL, adjust = "BH", n = 0, coef = NULL,
                              annot_df = NULL, intercept = FALSE) {
  ## Figure out the number of genes if not provided
  if (n == 0) {
    n <- nrow(fit[["coefficients"]])
  }

  ## If specific contrast(s) is/are not requested, get them all.
  if (is.null(coef)) {
    if (isTRUE(intercept)) {
      coef <- colnames(fit[["coefficients"]])
      coef <- coef[2:length(coef)]
    } else {
      coef <- colnames(fit[["contrasts"]])
    }
  } else {
    coef <- as.character(coef)
  }
  return_identities <- list()
  return_data <- list()
  end <- length(coef)
  data_tables <- list()
  if (isTRUE(intercept)) {

    ## If we do have an intercept model, then we get the data
    ## in a slightly different fashion.
    for (c in seq_len(ncol(fit[["coefficients"]]))) {
      data_table <-  limma::topTable(fit, adjust.method = adjust,
                                     n = n, coef = c, sort.by = "logFC")

      for (column in seq_len(ncol(data_table))) {
        data_table[[column]] <- signif(x = as.numeric(data_table[[column]]), digits = 4)
      }
      if (!is.null(annot_df)) {
        data_table <- merge(data_table, annot_df, by.x = "row.names", by.y = "row.names")
      }

      if (c == 1) {
        return_identities[[1]] <- data_table
      } else {
        comparison <- colnames(fit[["coefficients"]])[c]
        return_data[[comparison]] <- data_table
      }
    }

  } else {
    ## If we do not have an intercept (~ 0 + ...)
    ## Then extract the coefficients and identities separately.
    for (c in seq_len(end)) {
      comparison <- coef[c]
      message("Limma step 6/6: ", c, "/", end, ": Creating table: ",
              comparison, ".  Adjust = ", adjust)
      data_tables[[c]] <- limma::topTable(fit, adjust.method = adjust,
                                          n = n, coef = comparison, sort.by = "logFC")
      names(data_tables)[c] <- comparison
    }

    ## Take a moment to prettily format the numbers in the tables
    ## and fill in the identity table.
    for (d in seq_along(data_tables)) {
      comparison <- coef[d]
      table <- data_tables[[d]]
      for (column in seq_len(ncol(table))) {
        table[[column]] <- signif(x = as.numeric(table[[column]]), digits = 4)
      }
      if (!is.null(annot_df)) {
        table <- merge(table, annot_df, by.x = "row.names", by.y = "row.names")
      }
      if (grepl(pattern = "_vs_", x = comparison)) {
        return_data[[comparison]] <- table
      } else {
        return_identities[[comparison]] <- table
      }
    }

  } ## End checking for an intercept/nointercept model.

  retlist <- list(
    "identities" = return_identities,
    "contrasts" = return_data)
  return(retlist)
}

#' Writes out the results of a limma search using write_de_table()
#'
#' Looking to provide a single interface for writing tables from limma and friends.
#'
#' @param data Output from limma_pairwise()
#' @param ... Options for writing the xlsx file.
#' @seealso [write_de_table()]
#' @examples
#' \dontrun{
#'  finished_comparison = limma_pairwise(expressionset)
#'  data_list = write_limma(finished_comparison)
#' }
#' @export
write_limma <- function(data, ...) {
  result <- write_de_table(data, type = "limma", ...)
  return(result)
}

## EOF
