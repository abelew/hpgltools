## Use variancePartition's dream method with limma.

#' Set up a model matrix and set of contrasts for pairwise comparisons using
#' voom/limma via their modified functions from variancePartition.
#'
#' Creates the set of all possible contrasts and performs them using voom/limma.
#'
#' @param input Dataframe/vector or expt class containing count tables,
#'  normalization state, etc.
#' @param model_fstring Formula string describing the statistical
#'  model of interest.
#' @param null_fstring Formula string describing the null model.
#' @param model_svs Matrix of surrogates or character describing how
#'  to get them.
#' @param extra_contrasts Some extra contrasts to add to the list.
#'  This can be pretty neat, lets say one has conditions A,B,C,D,E
#'  and wants to do (C/B)/A and (E/D)/A or (E/D)/(C/B) then use this
#'  with a string like: "c_vs_b_ctrla = (C-B)-A, e_vs_d_ctrla = (E-D)-A,
#'  de_vs_cb = (E-D)-(C-B),"
#' @param annot_df Data frame for annotations.
#' @param libsize I've recently figured out that libsize is far more important
#'  than I previously realized.  Play with it here.
#' @param filter Filter the data before seeking SVs?
#' @param num_surrogates Number of SVs or way to guesstimate them.
#' @param limma_method Choose one of limma's lm methods.
#' @param limma_robust Make the significance estimation robust?
#' @param voom_norm Use this method to normalize the voom inputs.
#' @param limma_trend Add trend lines to limma's voom plot?
#' @param force Force data which may not be appropriate for limma into it?
#' @param keepers Perform an explicit set of contrasts instead of all.
#' @param keep_underscore Sanitize underscores from the model matrix?
#' @param adjust Use this p-value adjustment.
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
dream_pairwise <- function(input = NULL, model_fstring = "~ 0 + condition + batch",
                           null_fstring = "~", model_svs = NULL,
                           extra_contrasts = NULL, annot_df = NULL,
                           libsize = NULL, filter = TRUE,
                           num_surrogates = "be",
                           limma_method = "ls", limma_robust = FALSE, voom_norm = "none",
                           limma_trend = FALSE, force = FALSE, keepers = NULL,
                           keep_underscore = TRUE, adjust = "BH", ...) {
  arglist <- list(...)
  ## This is used in the invocation of a voom() implementation for normalization.
  ## This is for the eBayes() call.

  mesg("Starting dream pairwise comparisons.")
  current <- state(input)
  filteredp <- current[["filter"]]
  if ((is.null(filteredp) || filteredp == "raw") && isTRUE(filter)) {
    input <- sm(normalize(input, filter = filter))
  }

  input <- sanitize_se(input, keep_underscore = keep_underscore)
  input_data <- choose_binom_dataset(input, force = force)
  count_mtrx <- input_data[["data"]]
  fctrs <- get_formula_factors(model_fstring)
  condition_column <- fctrs[["factors"]][1]
  design <- pData(input)
  conditions <- droplevels(as.factor(design[[condition_column]]))
  batches <- droplevels(as.factor(design[["batch"]]))
  condition_table <- table(conditions)
  batch_table <- table(batches)
  condition_levels <- levels(conditions)

  ## The following small piece of logic is intended to handle situations where we use
  ## tximport for limma (kallisto/sailfish/salmon).
  if (is.null(input[["tximport"]])) {
    ## Adding an explicit as.data.frame() because otherwise this gets cast as an EList
    ## and fails the function 'varianceParition::filterObj' or whatever
    ## matrices are always class("matrix", "array") which when queried
    ## comes back as an EList, and so variancePartition expects to find the
    ## 'E' slot because it assumes it is therefore a DGEList; but no it is a matrix
    ## So, if I recast this as a dataframe, it is no longer multi-classed
    ## and gets passed through properly, which is a stupid solution to a stupid problem.
    data <- as.data.frame(input_data[["data"]])
  } else {
    data <- edgeR::DGEList(input[["tximport"]][["scaled"]][["counts"]])
    data <- edgeR::calcNormFactors(data)
  }

  if (is.null(libsize)) {
    libsize <- libsize(input)
  }

  mesg("Dream/limma step 1/6: choosing model.")
  ## for the moment, if someone choose an alt model, force it through.
  appended_fstring <- model_fstring
  if ("character" %in% class(model_svs)) {
    model_params <- adjuster_svs(input, model_fstring = model_fstring,
                                 null_fstring = null_fstring,
                                 model_svs = model_svs,
                                 num_surrogates = num_surrogates,
                                 filter = filter,
                                 ...)
    estimate_type <- model_svs
    model_svs <- model_params[["model_adjust"]]
    null_model <- model_params[["null_model"]]
    appended_fstring <- model_params[["appended_fstring"]]
    design <- pData(model_params[["modified_input"]])
  }
  model_mtrx <- model.matrix(as.formula(appended_fstring), data = design)
  fctrs <- get_formula_factors(model_fstring)
  ## Note, if we want to work like DESEq2, this should not be first, but last.
  contrast_factor <- fctrs[["contrast"]]
  simple_fstring <- glue("~ 0 + {contrast_factor}")
  model_formula <- as.formula(model_fstring)
  simple_model <- model.matrix(as.formula(simple_fstring), data = design)
  voom_plot <- NULL
  mesg("Dream/limma 2/6: Attempting voomWithDreamWeights.")
  tmp_file <- tmpmd5file(pattern = "voom_dream", fileext = ".png")
  this_plot <- png(filename = tmp_file)
  controlled <- dev.control("enable")
  voom_result <- variancePartition::voomWithDreamWeights(
    counts = data, formula = model_fstring,
    data = design, plot = TRUE)
  voom_plot <- grDevices::recordPlot()
  plotted <- dev.off()
  removed <- file.remove(tmp_file)
  one_replicate <- FALSE
  if (is.null(voom_result)) {
    ## Apparently voom returns null where there is only 1 replicate.
    message("voom returned null, I am not sure what will happen.")
    one_replicate <- TRUE
    voom_result <- data
  }

  ## Do the lmFit() using this model
  pairwise_fits <- NULL
  identity_fits <- NULL
  mesg("Dream/limma step 3/6: making limma and dream contrasts.")
  contrasts <- make_pairwise_contrasts(
    model = model_mtrx, conditions = conditions, contrast_factor = contrast_factor,
    extra_contrasts = extra_contrasts, keepers = keepers, keep_underscore = keep_underscore,
    do_identities = FALSE)
  all_pairwise_contrasts <- contrasts[["all_pairwise_contrasts"]]
  contrast_vector <- c()
  for (n in seq_along(contrasts[["all_pairwise"]])) {
    dream_contrast <- gsub(x = contrasts[["all_pairwise"]][n],
                           pattern = "\\,[[:space:]]*$", replacement = "")
    contrast_vector <- c(contrast_vector, dream_contrast)
  }
  varpart_contrasts <- variancePartition::makeContrastsDream(
    formula = as.formula(appended_fstring), data = design, contrasts = contrast_vector)
  mesg("Dream/limma step 4/6: Running dream.")
  fitted_data <- variancePartition::dream(
    exprObj = voom_result, formula = model_fstring, data = design, L = varpart_contrasts)
  mesg("Dream/limma step 4.2/6: Making identity contrasts.")
  identity_contrasts <- sm(make_pairwise_contrasts(
    model = simple_model, conditions = conditions,
    contrast_factor = contrast_factor,
    do_identities = TRUE, do_pairwise = FALSE,
    keep_underscore = keep_underscore))
  identities <- identity_contrasts[["all_pairwise_contrasts"]]
  mesg("Dream/limma step 4.5/6: Running dream for identities.")
  ## The following throws a warning: contrasts with only a single non-zero term
  ## are already evaluated by default.  This is true, but I have not yet
  ## figured out how to access them unless I explicitly ask; thus I am going
  ## to just do a suppressWarnings() until I get back to this and figure it out.
  identity_fits <- suppressWarnings(variancePartition::dream(
    exprObj = voom_result, formula = simple_fstring, data = design, L = identities))
  ##identity_fits <- limma::contrasts.fit(fit = fitted_data, contrasts = identities)
  mesg("Dream/limma step 5/6: Running eBayes.")
  if (isTRUE(one_replicate)) {
    all_pairwise_comparisons <- fitted_data[["coefficients"]]
    all_identity_comparisons <- identity_fits[["coefficients"]]
  }
  ## One might reasonably ask, wtf for the next few lines:
  ## Here is a snippet of the dream documentation:
  ## Since dream uses an estimated degrees of freedom value for each
  ## hypothsis test, the degrees of freedom is different for each gene
  ## here. Therefore, the t-statistics are not directly comparable since
  ## they have different degrees of freedom. In order to be able to
  ## compare test statistics, we report z.std which is the p-value
  ## transformed into a signed z-score. This can be used for downstream
  ## analysis.Note that if a random effect is not specified, dream()
  ## automatically uses lmFit(), but the user must run eBayes()
  ## afterward.
  all_tables <- NULL
  all_pairwise_comparisons <- variancePartition::eBayes(fitted_data,
                                                        robust = limma_robust,
                                                        trend = limma_trend)
  all_identity_comparisons <- limma::eBayes(identity_fits,
                                            robust = limma_robust,
                                            trend = limma_trend)
  mesg("Dream/limma step 6/6: Creating tables.")
  ## Make a list of the output, one element for each comparison of the contrast matrix
  pairwise_results <- make_varpart_tables(fit = all_pairwise_comparisons,
                                          adjust = adjust, n = 0, coef = NULL,
                                          annot_df = NULL)
  varpart_tables <- pairwise_results[["contrasts"]]
  ## Stop limma from complaining about already evaluated tables.
  identity_results <- suppressWarnings(
    make_limma_tables(fit = all_identity_comparisons, adjust = "BH",
                      n = 0, coef = NULL, annot_df = NULL))
  varpart_identities <- identity_results[["identities"]]

  contrasts_performed <- names(varpart_tables)
  retlist <- list(
    "all_pairwise" = all_pairwise,
    "all_tables" = varpart_tables,
    "batches" = batches,
    "batch_table" = batch_table,
    "conditions" = conditions,
    "condition_table" = condition_table,
    "contrast_string" = contrast_vector,
    "contrasts_performed" = contrasts_performed,
    "design" = design,
    "dispersion_plot" = voom_plot,
    "fit" = fitted_data,
    "identities" = identities,
    "identity_tables" = varpart_identities,
    "identity_comparisons" = all_identity_comparisons,
    "input_data" = input,
    "method" = "varpart",
    "model" = model_mtrx,
    "model_string" = model_fstring,
    "pairwise_comparisons" = all_pairwise_comparisons,
    "single_table" = all_tables,
    "voom_result" = voom_result)
  class(retlist) <- c("dream_pairwise", "list")
  if (!is.null(arglist[["limma_excel"]])) {
    retlist[["dream_excel"]] <- write_limma(retlist, excel = arglist[["limma_excel"]])
  }
  return(retlist)
}

#' Print a summary of the result from dream_pairwise().
#'
#' @param x List from dream_pairwise().
#' @param ... Other args for the generic.
#' @export
print.dream_pairwise <- function(x, ...) {
  summary_string <- glue("The results from the hybrid variancePartition/limma pairwise analysis.")
  message(summary_string)
  return(invisible(x))
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
make_varpart_tables <- function(fit = NULL, adjust = "BH", n = 0, coef = NULL,
                                annot_df = NULL, intercept = FALSE) {
  ## Figure out the number of genes if not provided
  if (n == 0 || is.null(n)) {
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
      data_table <- variancePartition::topTable(fit, adjust.method = adjust,
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
      comp_name <- strsplit(x = comparison, split = " = ")[[1]][1]
      mesg("Varpart/limma step 6/6: ", c, "/", end, ": Creating table: ",
           comp_name, ".  Adjust = ", adjust)
      data_tables[[comp_name]] <- variancePartition::topTable(
        fit, adjust.method = adjust,
        n = n, coef = comparison, sort.by = "logFC")
    }
    ## Take a moment to prettily format the numbers in the tables
    ## and fill in the identity table.
    for (d in seq_along(data_tables)) {
      comparison <- coef[d]
      comp_name <- strsplit(x = comparison, split = " = ")[[1]][1]
      table <- data_tables[[d]]
      for (column in seq_len(ncol(table))) {
        table[[column]] <- signif(x = as.numeric(table[[column]]), digits = 4)
      }
      if (!is.null(annot_df)) {
        table <- merge(table, annot_df, by.x = "row.names", by.y = "row.names")
      }
      if (grepl(pattern = "_vs_", x = comparison)) {
        return_data[[comp_name]] <- table
      } else {
        return_identities[[comp_name]] <- table
      }
    }
  } ## End checking for an intercept/nointercept model.

  retlist <- list(
    "identities" = return_identities,
    "contrasts" = return_data)
  return(retlist)
}

#' Nearly a copy of write_limma().
#'
#' This will add a couple of columns to the output table which are
#' specific to variancePartition's dream
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
write_dream <- function(data, ...) {
  result <- write_de_table(data, type = "dream", ...)
  return(result)
}

## EOF
