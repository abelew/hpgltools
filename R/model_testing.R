## model_testing.r: Some functions to get ready to pass models to various
## DE/test/etc methods.  These functions seek to catch some corner cases when
## playing with the various model types when using DESeq/etc.


#' Try out a few experimental models and return a likely working
#' option.
#'
#' This function is simultaneously too complicated and too stupid.
#'
#' The _pairwise family of functions all demand an experimental model.  This
#' tries to choose a consistent and useful model for all for them.  This does
#' not try to do multi-factor, interacting, nor dependent variable models, if
#' you want those do them yourself and pass them off as alt_model.
#'
#' Invoked by the _pairwise() functions.
#'
#' @param input Input data used to make the model.
#' @param conditions Factor of conditions in the putative model.
#' @param batches Factor of batches in the putative model.
#' @param model_batch Try to include batch in the model?
#' @param model_cond Try to include condition in the model? (Yes!)
#' @param model_intercept Use an intercept model instead of cell-means?
#' @param alt_model Use your own model.
#' @param alt_string String describing an alternate model.
#' @param intercept Choose an intercept for the model as opposed to 0.
#' @param reverse Reverse condition/batch in the model?  This shouldn't/doesn't
#'  matter but I wanted to test.
#' @param contr List of contrasts.arg possibilities.
#' @param surrogates Number of or method used to choose the number of surrogate
#'  variables.
#' @param verbose Print some information about what is happening?
#' @param ... Further options are passed to arglist.
#' @return List including a model matrix and strings describing cell-means and
#'  intercept models.
#' @seealso [stats::model.matrix()]
#' @examples
#' \dontrun{
#'  a_model <- choose_model(expt, model_batch = TRUE, model_intercept = FALSE)
#'  a_model$chosen_model
#'  ## ~ 0 + condition + batch
#' }
old_choose_model <- function(input, conditions = NULL, batches = NULL, model_batch = TRUE,
                         model_cond = TRUE, model_intercept = FALSE,
                         alt_model = NULL, alt_string = NULL, null_model = NULL,
                         intercept = 0, reverse = FALSE, contr = NULL,
                         surrogates = "be", verbose = TRUE, keep_underscore = FALSE, ...) {
  arglist <- list(...)
  design <- NULL
  if (class(input)[1] != "matrix" && class(input)[1] != "data.frame") {
    design <- pData(input)
  }
  if (is.null(design)) {
    conditions <- as.factor(conditions)
    batches <- as.factor(batches)
    design <- data.frame("condition" = conditions,
                         "batch" = batches,
                         stringsAsFactors = TRUE)
    message("Design was null, but is now:")
    print(design)
  }
  ## Make a model matrix which will have one entry for
  ## each of the condition/batches
  ## It would be much smarter to generate the models in the following if() {} blocks
  ## But I have it in my head to eventually compare results using different models.

  ## colname_exclude: used to simplify the model column names
  ## so that we need not do stupid things like set up contrasts with names like:
  ## 'visitnumbervisit3 - visitnumbervisit1' but may instead just do 'visit3 - visit1'
  ## This defaults to condition, but in the case of an alt model will be the first
  ## factor mentioned in the model

  ## The previous iteration of this had an explicit contrasts.arg set, like this:
  ## contrasts.arg = list(condition = "contr.treatment"))
  ## Which looked like this for a full invocation:
  ## condbatch_int_model <- try(stats::model.matrix(~ 0 + condition + batch,
  ##                                   contrasts.arg = list(condition = "contr.treatment",
  ##                                                      batch = "contr.treatment")),
  ## The contrasts.arg has been removed because it seems to result in the same model.

  ## Note, I do not have contrasts.arg here,
  ## I should make it the first element of the int_string
  model_colname_exclude <- "condition"
  contrast_element <- "condition"
  clist <- list("condition" = "contr.treatment")
  blist <- list("batch" = "contr.treatment")
  cblist <- list("condition" = "contr.treatment", "batch" = "contr.treatment")
  if (!is.null(contr)) {
    if (!is.null(contr[["condition"]]) && !is.null(contr[["batch"]])) {
      cblist <- list("condition" = contr[["condition"]], "batch" = contr[["batch"]])
    } else if (!is.null(contr[["condition"]])) {
      clist <- list("condition" = contr[["condition"]])
      cblist[["condition"]] <- contr[["condition"]]
    } else if (!is.null(contr[["batch"]])) {
      blist <- list("batch" = contr[["batch"]])
      cblist[["batch"]] <- contr[["batch"]]
    }
  }

  cond_noint_string <- "~ 0 + condition"
  cond_noint_model <- try(stats::model.matrix(
    object = as.formula(cond_noint_string),
    contrasts.arg = clist,
    data = design))

  batch_noint_string <- "~ 0 + batch"
  batch_noint_model <- try(stats::model.matrix(
    object = as.formula(batch_noint_string),
    contrasts.arg = blist,
    data = design))

  condbatch_noint_string <- "~ 0 + condition + batch"
  condbatch_noint_model <- try(stats::model.matrix(
    object = as.formula(condbatch_noint_string),
    contrasts.arg = cblist,
    data = design))

  batchcond_noint_string <- "~ 0 + batch + condition"
  batchcond_noint_model <- try(stats::model.matrix(
    object = as.formula(batchcond_noint_string),
    contrasts.arg = cblist,
    data = design))

  cond_int_string <- "~ condition"
  cond_int_model <- try(stats::model.matrix(
    object = as.formula(cond_int_string),
    contrasts.arg = clist,
    data = design))

  batch_int_string <- "~ batch"
  batch_int_model <- try(stats::model.matrix(
    object = as.formula(batch_int_string),
    contrasts.arg = blist,
    data = design))

  condbatch_int_string <- "~ condition + batch"
  condbatch_int_model <- try(stats::model.matrix(
    object = as.formula(condbatch_int_string),
    contrasts.arg = cblist,
    data = design), silent = TRUE)

  batchcond_int_string <- "~ batch + condition"
  batchcond_int_model <- try(stats::model.matrix(
    object = as.formula(batchcond_int_string),
    contrasts.arg = cblist,
    data = design), silent = TRUE)

  ## Default to the conditional model
  int_model <- cond_int_model
  noint_model <- cond_noint_model
  int_string <- cond_int_string
  noint_string <- cond_noint_string
  including <- "condition"
  if (!is.null(alt_model)) {
    message("We have an alt model.")
    chosen_model <- stats::model.matrix(object = as.formula(alt_model),
                                        data = design)
    if (!is.null(contr)) {
      chosen_model <- stats::model.matrix(object = as.formula(alt_model),
                                          contrasts.arg = contr,
                                          data = design)
    }
    int_model <- chosen_model
    noint_model <- chosen_model
    int_string <- alt_model
    noint_string <- alt_model
    including <- "alt"
  } else if (is.null(model_batch)) {
    message("Model_batch is null.")
    int_model <- cond_int_model
    noint_model <- cond_noint_model
    int_string <- cond_int_string
    noint_string <- cond_noint_string
    including <- "condition"
  } else if (isTRUE(model_cond) && isTRUE(model_batch)) {
    message("model_cond and model_batch are true.")
    if (class(condbatch_int_model)[1] == "try-error") {
      if (isTRUE(verbose)) {
        message("The condition+batch model failed. ",
                "Does your experimental design support both condition and batch? ",
                "Using only a conditional model.")
      }
      int_model <- cond_int_model
      noint_model <- cond_noint_model
      int_string <- cond_int_string
      noint_string <- cond_noint_string
      including <- "condition"
    } else if (isTRUE(reverse)) {
      message("Reversing the condition+batch model to batch+condition.")
      int_model <- batchcond_int_model
      noint_model <- batchcond_noint_model
      int_string <- batchcond_int_string
      noint_string <- batchcond_noint_string
      including <- "batch+condition"
    } else {
      message("Using condition+batch.")
      int_model <- condbatch_int_model
      noint_model <- condbatch_noint_model
      int_string <- condbatch_int_string
      noint_string <- condbatch_noint_string
      including <- "condition+batch"
    }
  }
  message("Finished initial if, int_string: ", int_string)

  if (class(model_batch)[1] == "character") {
    ## Then calculate the estimates using all_adjusters
    if (isTRUE(verbose)) {
      mesg("Extracting surrogate estimates from ", model_batch,
              " and adding them to the model.")
    }
    model_batch_info <- all_adjusters(input, estimate_type = model_batch,
                                      surrogates = surrogates,
                                      alt_model = alt_model)
    ## Changing model_batch from 'sva' to the resulting matrix.
    ## Hopefully this will simplify things later for me.
    model_batch <- model_batch_info[["model_adjust"]]
    int_append_string <- paste0(int_string, " + model_batch")

    contrast_element <- gsub(x = noint_string, pattern = "~\\s*\\d*\\s*\\+\\s*(\\w+).*$", replacement = "\\1")
    model_colname_exclude <- contrast_element
    message("Setting the contrasts.arg to: ", contrast_element, ".")
    clist <- list()
    clist[[contrast_element]] <- "contr.treatment"

    int_model <- stats::model.matrix(as.formula(int_append_string),
                                     contrasts.arg = clist,
                                     data = design)
    noint_append_string <- paste0(noint_string, " + model_batch")
    noint_model <- stats::model.matrix(as.formula(noint_append_string),
                                       contrasts.arg = clist,
                                       data = design)
    sv_names <- glue("SV{1:ncol(model_batch)}")
    sv_string <- ""
    for (sv in sv_names) {
      sv_string <- glue("{sv_string} + {sv}")
    }
    noint_string <- glue("{noint_string}{sv_string}")
    int_string <- glue("{int_string}{sv_string}")
    rownames(model_batch) <- rownames(int_model)
    including <- glue("{contrast_element}{sv_string}")
    message("Finishing alt_model with: ", int_string, " and ", noint_string)
  } else if (class(model_batch)[1] == "numeric" || class(model_batch)[1] == "matrix") {
    if (isTRUE(verbose)) {
      mesg("Including batch estimates from sva/ruv/pca in the model.")
    }
    int_append_string <- paste0(int_string, " + model_batch")
    ## Note, I do not have contrasts.arg here,
    ## I should make it the first element of the int_string
    contrast_element <- gsub(x = noint_string,
                             pattern = "~\\s*\\d*\\s*\\+\\s*(\\w+).*$",
                             replacement = "\\1")
    model_colname_exclude <- contrast_element
    message("Setting the contrasts.arg to: ", contrast_element, ".")
    clist <- list()
    clist[[contrast_element]] <- "contr.treatment"
    int_model <- stats::model.matrix(as.formula(int_append_string),
                                     contrasts.arg = clist,
                                     data = design)
    noint_append_string <- paste0(noint_string, " + model_batch")
    noint_model <- stats::model.matrix(as.formula(noint_append_string),
                                       contrasts.arg = clist,
                                       data = design)
    sv_names <- glue("SV{1:ncol(model_batch)}")
    sv_string <- ""
    for (sv in sv_names) {
      sv_string <- glue("{sv_string} + {sv}")
    }
    int_string <- glue("{int_string}{sv_string}")
    noint_string <- glue("{noint_string}{sv_string}")
    rownames(model_batch) <- rownames(int_model)
    including <- glue("{contrast_element}{sv_string}")
  } else if (isTRUE(model_cond)) {
    int_model <- cond_int_model
    noint_model <- cond_noint_model
    int_string <- cond_int_string
    noint_string <- cond_noint_string
    including <- "condition"
  } else if (isTRUE(model_batch)) {
    int_model <- batch_int_model
    noint_model <- batch_noint_model
    int_string <- batch_int_string
    noint_string <- batch_noint_string
    including <- "batch"
  } else {
    ## Default to the conditional model
    int_model <- cond_int_model
    noint_model <- cond_noint_model
    int_string <- cond_int_string
    noint_string <- cond_noint_string
    including <- "condition"
  }

  tmpnames <- colnames(int_model)
  tmpnames <- gsub(pattern = "model_batch", replacement = "", x = tmpnames)
  if (isTRUE(keep_underscore)) {
    tmpnames <- gsub(pattern = "data[^_[:^punct:]]", replacement = "", x = tmpnames, perl = TRUE)
  } else {
    tmpnames <- gsub(pattern = "data[[:punct:]]", replacement = "", x = tmpnames)
  }
  tmpnames <- gsub(pattern = "-", replacement = "", x = tmpnames)
  tmpnames <- gsub(pattern = "+", replacement = "", x = tmpnames)
  ## The next lines ensure that conditions/batches which are all numeric will
  ## not cause weird errors for contrasts. Ergo, if a condition is something
  ## like '111', now it will be 'c111' Similarly, a batch '01' will be 'b01'
  numeric_pattern <- glue("^{model_colname_exclude}(\\\\d+)$")
  tmpnames <- gsub(pattern = numeric_pattern, replacement = "c\\1", x = tmpnames)
  tmpnames <- gsub(pattern = "^batch(\\d+)$", replacement = "b\\1", x = tmpnames)
  tmpnames <- gsub(pattern = model_colname_exclude, replacement = "", x = tmpnames)
  tmpnames <- gsub(pattern = "batch", replacement = "", x = tmpnames)
  colnames(int_model) <- tmpnames

  tmpnames <- colnames(noint_model)
  tmpnames <- gsub(pattern = "model_batch", replacement = "", x = tmpnames)
  if (isTRUE(keep_underscore)) {
    tmpnames <- gsub(pattern = "data[^_[:^punct:]]", replacement = "", x = tmpnames, perl = TRUE)
  } else {
    tmpnames <- gsub(pattern = "data[[:punct:]]", replacement = "", x = tmpnames)
  }
  tmpnames <- gsub(pattern = "-", replacement = "", x = tmpnames)
  tmpnames <- gsub(pattern = "+", replacement = "", x = tmpnames)
  ## The next lines ensure that conditions/batches which are all numeric will
  ## not cause weird errors for contrasts. Ergo, if a condition is something
  ## like '111', now it will be 'c111' Similarly, a batch '01' will be 'b01'
  tmpnames <- gsub(pattern = numeric_pattern, replacement = "c\\1", x = tmpnames)
  tmpnames <- gsub(pattern = "batch^(\\d+)$", replacement = "b\\1", x = tmpnames)
  tmpnames <- gsub(pattern = model_colname_exclude, replacement = "", x = tmpnames)
  tmpnames <- gsub(pattern = "batch", replacement = "", x = tmpnames)
  colnames(noint_model) <- tmpnames

  chosen_model <- NULL
  chosen_string <- NULL
  if (isTRUE(model_intercept)) {
    if (isTRUE(verbose)) {
      mesg("Choosing the intercept containing model: ", int_string, ".")
    }
    chosen_model <- int_model
    chosen_string <- int_string
  } else {
    if (isTRUE(verbose)) {
      mesg("Choosing the non-intercept containing model: ", noint_string, ".")
    }
    chosen_model <- noint_model
    chosen_string <- noint_string
  }

  retlist <- list(
    "int_model" = int_model,
    "noint_model" = noint_model,
    "int_string" = int_string,
    "noint_string" = noint_string,
    "chosen_model" = chosen_model,
    "chosen_string" = chosen_string,
    "model_batch" = model_batch,
    "including" = including)
  return(retlist)
}

choose_alt_model <- function(input, alt_model, surrogates, keep_underscore, contr,
                             ...) {
  arglist <- list(...)
  message("Working inside choose_alt_model")
  design <- pData(input)
  model_factors <- get_formula_factors(alt_model)
  ## Remove this string from the colnames of the model:
  contrast_factor <- model_factors[["factors"]][1]
  contr_treatment_list <- list()
  contr_treatment_list[[contrast_factor]] = "contr.treatment"
  model_formula <- as.formula(alt_model)
  data_model <- stats::model.matrix(object = model_formula,
                                    contrasts.arg = contr_treatment_list,
                                    data = design)
  null_formula_string <- make_null_formula(alt_model)

  if (!is.null(surrogate_estimate)) {
    if (class(surrogate_estimate)[1] == "character") {
      ## Then calculate the estimates using all_adjusters
      if (isTRUE(verbose)) {
        mesg("Extracting surrogate estimates from ", surrogate_estimate,
             " and adding them to the model.")
      }
      surrogate_info <- all_adjusters(input, design = design, estimate_type = surrogate_estimate,
                                      batch1 = "batch", batch2 = NULL,
                                      surrogates = surrogates, low_to_zero = FALSE,
                                      na_to_zero = TRUE, adjust_method = NULL)

      ## Changing model_batch from 'sva' to the resulting matrix.
      ## Hopefully this will simplify things later for me.
      surrogate_estimate <- surrogate_info[["model_adjust"]]
      rownames(surrogate_estimate) <- rownames(data_model)
      ## Now we have a matrix of surrogate estimates which may be appended to the model.
    }

    sv_names <- glue("SV{seq_len(ncol(surrogate_estimate))}")
    colnames(surrogate_estimate) <- sv_names
    sv_string <- ""
    for (sv in sv_names) {
      sv_string <- glue("{sv_string} + {sv}")
    }
    including <- glue("{including}{sv_string}")
    message("Adding surrogates to the model via: ", sv_string, ".")
    model_string <- glue("{model_string} {sv_string}")
    data_model <- cbind(data_model, surrogate_estimate)
  }
  data_model <- sanitize_model(data_model)

  retlist <- list(
    "chosen_model" = data_model,
    "chosen_formula" = model_formula,
    "chosen_string" = model_string,
    "including" = including)
  return(retlist)
}

choose_model <- function(input, conditions = NULL, batches = NULL, model_batch = TRUE,
                         model_cond = TRUE, model_intercept = FALSE,
                         alt_model = NULL, null_model = NULL,
                         intercept = 0, contr = NULL, surrogate_estimate = NULL,
                         surrogates = "be", keep_underscore = FALSE, ...) {
  arglist <- list(...)
  ## Make a model matrix which will have one entry for
  ## each of the condition/batches
  if (isTRUE(model_intercept)) {
    mesg("Creating a default intercept model starting with ~ condition,")
  } else {
    mesg("Creating a default cell-means model starting with ~ 0 + condition.")
  }

  design <- pData(input)
  model_colname_exclude <- "condition"
  contrast_element <- "condition"
  contr_treatment_list <- list()
  if (is.null(contr)) {
    if (isTRUE(model_cond)) {
      contr_treatment_list[["condition"]] = "contr.treatment"
    }
    if (isTRUE(model_batch)) {
      contr_treatment_list[["batch"]] = "contr.treatment"
    }
  } else {
    if (!is.null(contr[["condition"]])) {
      contr_treatment_list[["condition"]] = contr[["condition"]]
    }
    if (!is.null(contr[["batch"]])) {
      contr_treatment_list[["batch"]] = contr[["batch"]]
    }
  }

  including <- ""
  model_string <- "~"
  ## This is another place we can simplify.
  if (isFALSE(model_intercept)) {
    model_string <- glue("{model_string} {intercept} +")
  }
  if (isTRUE(model_cond)) {
    model_string <- glue("{model_string} condition +")
    including <- glue("{including}condition ")
  }
  if (isTRUE(model_batch)) {
    model_string <- glue("{model_string} batch +")
    including <- glue("{including}batch ")
  }
  model_string <- gsub(x = model_string, pattern = "\\+$",
                       replacement = "")
  including <- gsub(x = including, pattern = " $", replacement = "")
  model_formula <- as.formula(model_string)
  data_model <- stats::model.matrix(object = model_formula,
                                    contrasts.arg = contr_treatment_list,
                                    data = design)

  if (!is.null(surrogate_estimate)) {
    if (class(surrogate_estimate)[1] == "character") {
      ## Then calculate the estimates using all_adjusters
      if (isTRUE(verbose)) {
        mesg("Extracting surrogate estimates from ", surrogate_estimate,
             " and adding them to the model.")
      }
      surrogate_info <- all_adjusters(input, design = design, estimate_type = surrogate_estimate,
                                      batch1 = "batch", batch2 = NULL,
                                      surrogates = surrogates, low_to_zero = FALSE,
                                      na_to_zero = TRUE, adjust_method = NULL)

      ## Changing model_batch from 'sva' to the resulting matrix.
      ## Hopefully this will simplify things later for me.
      surrogate_estimate <- surrogate_info[["model_adjust"]]
      rownames(surrogate_estimate) <- rownames(data_model)
      ## Now we have a matrix of surrogate estimates which may be appended to the model.
    }

    sv_names <- glue("SV{seq_len(ncol(surrogate_estimate))}")
    colnames(surrogate_estimate) <- sv_names
    sv_string <- ""
    for (sv in sv_names) {
      sv_string <- glue("{sv_string} + {sv}")
    }
    including <- glue("{including}{sv_string}")
    message("Adding surrogates to the model via: ", sv_string, ".")
    model_string <- glue("{model_string} {sv_string}")
    data_model <- cbind(data_model, surrogate_estimate)
    appended_design <- cbind(design, surrogate_estimate)
  } ## End testing for a surrogate estimate
  data_model <- sanitize_model(data_model)

  ## 202407: Include a new key, 'model_surrogates' to make it easier for
  ## things like DESeq to append surrogates to the pData.
  ## or perhaps just add the svs to the input and return it?
  retlist <- list(
    "chosen_model" = data_model,
    "chosen_formula" = model_formula,
    "chosen_string" = model_string,
    "model_surrogates" = surrogate_estimate,
    "input" = input,
    "including" = including)
  return(retlist)
}

#' A single place to extract count tables from a set of surrogate variables.
#'
#' Given an initial set of counts and a series of surrogates, what would the
#' resulting count table look like? Hopefully this function answers that
#' question.
#'
#' @param data Original count table, may be an expt/expressionset or df/matrix.
#' @param adjust Surrogates with which to adjust the data.
#' @param design Experimental design if it is not included in the expressionset.
#' @param method Which methodology to follow, ideally these agree but that seems untrue.
#' @param cond_column design column containing the condition data.
#' @param batch_column design column with the batch data, used for
#'  subtractive methods.
#' @param matrix_scale Was the input for the surrogate estimator on a log or linear scale?
#' @param return_scale Does one want the output linear or log?
#' @param ... Arguments passed to downstream functions.
#' @return A data frame of adjusted counts.
#' @seealso [sva] [RUVSeq] [crossprod()] [tcrossprod()] [solve()]
#' @export
counts_from_surrogates_model <- function(mtrx, adjust = NULL, method = "ruv", design = NULL,
                                         full_model = NULL, cond_model = NULL, null_model = NULL,
                                         matrix_scale = "linear", return_scale = "linear") {
  adjust_mtrx <- as.matrix(adjust)
  switchret <- switch(
    method,
    "cbcb_add" = {
      ## we also have conditional_model which may make this easier
      new_counts <- cbcb_batch(mtrx, full_model,
                               conditional_model = cond_model,
                               batch_model = null_model,
                               method = "add",
                               matrix_scale = matrix_scale,
                               return_scale = return_scale,
                               ...)
    },
    "cbcb_subtract" = {
      ##new_counts <- cbcb_batch(mtrx, my_design, method = "subtract",
      ##                         matrix_scale = matrix_scale, return_scale = return_scale,
      ##                         ...)
      new_counts <- cbcb_batch(mtrx, new_model,
                               conditional_model = cond_model,
                               batch_model = batch_model,
                               method = "subtract",
                               matrix_scale = matrix_scale, return_scale = return_scale,
                               ...)
    },
    "ruv" = {
      ## Here is the original code, as a reminder: W is the matrix of surrogates.
      ## W <- svdWa$u[, (first:k), drop = FALSE]
      ## alpha <- solve(t(W) %*% W) %*% t(W) %*% Y
      ## correctedY <- Y - W %*% alpha
      ## if(!isLog & all(.isWholeNumber(x))) {
      ##   if(round) {
      ##     correctedY <- round(exp(correctedY) - epsilon)
      ##     correctedY[correctedY<0] <- 0
      ##   } else {
      ##     correctedY <- exp(correctedY) - epsilon
      ##   }
      ## }
      ## colnames(W) <- paste("W", seq(1, ncol(W)), sep = "_")
      ## return(list(W = W, normalizedCounts = t(correctedY)))

      ##alpha <- try(solve(t(adjust_mtrx) %*% adjust_mtrx))
      ## Y <- t(mtrx)
      ## W <- ruv_model

      log_mtrx <- NULL
      if (matrix_scale != "log2") {
        ## Note that we are actually going on scale e.
        log_mtrx <- log(mtrx + 1)
      }

      ## original_alpha <- solve(t(adjust_mtrx) %*% adjust_mtrx) %*% t(adjust_mtrx) %*% t(data_mtrx)
      alpha <- try(solve(crossprod(adjust_mtrx)))
      if (class(alpha)[1] == "try-error") {
        warning("Data modification by the model failed. Leaving counts untouched.")
        return(mtrx)
      }
      beta <- tcrossprod(t(adjust_mtrx), log_mtrx)
      gamma <- alpha %*% beta
      delta <- t(log_mtrx) - (adjust_mtrx %*% gamma)
      new_counts <- t(delta)

      if (return_scale == "log2") {
        new_counts <- new_counts / log(2)
      } else if (return_scale == "linear") {
        new_counts <- exp(new_counts) - 1
      } else {
        stop("I do not understand the scale: ", return_scale, ".")
      }
    },
    "solve_crossproducts" = {
      ## I think that if I did the math out, this is the same as ruv above.

      ##data_modifier <- try(solve(t(new_model) %*% new_model) %*% t(new_model))
      ## In the previous code, this was: 'Hat <- solve(t(X) %*% X) %*% t(X)'
      ## Now it is in two separate lines, first the solve operation:

      if (matrix_scale != "log2") {
        ## Note that we are actually going on scale e.
        l2mtrx <- log2(mtrx + 1)
      }

      data_solve <- try(solve(t(new_model) %*% new_model), silent = TRUE)
      if (class(data_solve)[1] == "try-error") {
        warning("Data modification by the model failed. Leaving counts untouched.")
        return(mtrx)
      }
      ## If the solve operation passes, then the '%*% t(X)' is allowed to happen.
      data_modifier <- data_solve %*% t(new_model)
      transformation <- (data_modifier %*% t(l2mtrx))
      conds <- ncol(conditional_model)
      new_counts <- l2mtrx - t(as.matrix(new_model[, -c(1:conds)]) %*%
                                 transformation[-c(1:conds), ])
      if (return_scale == "linear") {
        new_counts <- (2 ^ new_counts) - 1
      } else if (return_scale != "log2") {
        stop("I do not understand the return scale: ", return_scale, ".")
      }
    },
    {
      stop("I do not understand method: ", method, ".")
    }
  ) ## End the switch

  ## If the matrix state and return state are not the same, fix it.
  ## It appears to me that the logic of this is wrong, but I am not yet certain why.
  return(new_counts)
}





#' Simplify a model formula string to a simple vector of factors.
#'
#' I would like to be able to make my various pairwise methods
#' generalizable, therefore I need to no longer rely on an expected
#' pair of factors: 'condition' and 'batch', instead I want to observe
#' the first factor and anything which follows.  Thus, I will string
#' split the formula string by punctuation and yank out the names of
#' interest.
#'
#' @param formula_string Formula describing the formula of interest.
#' @return Ordered vector of factors in the formula ignoring
#' interactions etc.
#' @export
get_formula_factors <- function(formula_string = NULL) {
  if (is.null(formula_string)) {
    warning("No formula string was provided, this is an error, but for now this will assume condition+batch.")
    formula_string = "~ 0 + condition + batch"
  }
  message("Getting factors from: ", formula_string, ".")
  split_result <- strsplit(x = formula_string, split = "[[:punct:]]", perl = TRUE)
  factor_vector <- split_result[[1]]
  retlist <- list(
    "type" = "cellmeans",
    "interaction" = FALSE,
    "mixed" = FALSE)
  interaction_test <- grepl(x = formula_string, pattern = "\\*|\\:")
  if (isTRUE(interaction_test)) {
    retlist[["interaction"]] <- TRUE
    interactor_split <- strsplit(x = formula_string, split = "\\*|\\:")[[1]]
    interactors <- c("first", "second")
    first_interactor <- strsplit(x = interactor_split[1], split = "[[:punct:]|[:space:]]")[[1]]
    interactors[1] <- first_interactor[length(first_interactor)]
    second_interactor <- strsplit(x = interactor_split[2], split = "[[:punct:]|[:space:]]")[[1]]
    interactors[2] <- second_interactor[1]
    retlist[["interactors"]] <- interactors
  }
  mixed_test <- grepl(x = formula_string, pattern = "\\|")
  if (isTRUE(mixed_test)) {
    retlist[["mixed"]] <- TRUE
    mixed_split <- strsplit(x = formula_string, split = "\\|")[[1]]
    mixers <- c("first", "second")
    first_mix <- strsplit(x = mixed_split[1], split = "[[:punct:]|[:space:]]")[[1]]
    mixers[1] <- first_mix[length(first_mix)]
    second_mix <- strsplit(x = mixed_split[2], split = "[[:punct:]|[:space:]]")[[1]]
    mixers[2] <- second_mix[1]
    retlist[["mixers"]] <- mixers
  }

  fct_vector <- c()
  for (i in seq_len(length(factor_vector))) {
    fct <- gsub(x = factor_vector[i], pattern = "[[:space:]]", replacement = "")
    ## Check for an cell-means model intercept: ~ 0 + f1 + f2
    ##                                            ^
    if (grepl(x = fct, pattern = "^[[:digit:]]+$")) {
      retlist[["type"]] <- "cellmeans"
      retlist[["cellmeans_intercept"]] <- fct
      next
    }

    if (grepl(x = fct, pattern = "^[[:alnum:]]+$")) {
      fct_vector <- c(fct_vector, fct)
    }
  }
  retlist[["factors"]] <- unique(fct_vector)
  return(retlist)
}

make_null_formula <- function(formula_string, starting = 2) {
  formula_info <- get_formula_factors(formula_string)
  factors <- formula_info[["factors"]]
  message("Extracting a new formula string starting with ", factors[starting], ".")

  start <- formula_info[["cellmeans_intercept"]]
  null_string <- "~"
  if (formula_info[["type"]] == "cellmeans") {
    null_string <- glue("{null_string} {start}")
  }
  if (starting <= length(factors)) {
    for (f in starting:length(factors)) {
      fct <- factors[f]
      message("Adding fct: ", fct)
      null_string <- glue("{null_string} + {fct}")
    }
  }
  return(null_string)
}

#' Make sure a given experimental factor and design will play together.
#'
#' Have you ever wanted to set up a differential expression analysis and after
#' minutes of the computer churning away it errors out with some weird error
#' about rank?  Then this is the function for you!
#'
#' @param design Dataframe describing the design of the experiment.
#' @param goal Experimental factor you actually want to learn about.
#' @param factors Experimental factors you rather wish would just go away.
#' @param ... I might decide to add more options from other functions.
#' @return List of booleans telling if the factors + goal will work.
#' @seealso [model.matrix()] [qr()]
#' @export
model_test <- function(design, goal = "condition", factors = NULL, ...) {
  arglist <- list(...)
  ## For testing, use some existing matrices/data
  message("There are ", length(levels(as.factor(design[, goal]))),
          " levels in the goal: ", goal, ".")
  ret_list <- list()
  if (is.null(factors)) {
    message("Testing an experimental design with only ", goal, ".")
    matrix_all_formula <- as.formula(glue("~ 0 + {goal}"))
    matrix_test <- model.matrix(matrix_all_formula, data = design)
    num_columns <- ncol(matrix_test)
    matrix_decomp <- qr(matrix_test)
    message("The model of ", goal, "has ", num_columns,
            " levels and rank ", matrix_decomp[["rank"]], ".")
    if (matrix_decomp[["rank"]] < num_columns) {
      message("This will not work, a different factor should be used.")
      ret_list[[goal]] <- FALSE
    } else {
      ret_list[[goal]] <- TRUE
    }
    for (factor in colnames(design)) {
      if (factor == goal) {
        next
      }
      message("Testing an experimental design with ", goal, " and ", factor, ".")
      matrix_goal <- design[, goal]
      matrix_factor <- design[, factor]
      if (length(levels(as.factor(matrix_factor))) == 1) {
        message("Factor ", factor, " has only 1 level, skipping it.")
        next
      }
      matrix_all_formula <- as.formula(glue("~ 0 + {goal} + {factor}"))
      matrix_test <- model.matrix(matrix_all_formula, data = design)
      num_columns <- ncol(matrix_test)
      matrix_decomp <- qr(matrix_test)
      message("The model of ", goal, " and ", factor, " has ", num_columns,
              " and rank ", matrix_decomp[["rank"]])
      if (matrix_decomp[["rank"]] < num_columns) {
        message("This will not work, a different factor should be used.")
        ret_list[[factor]] <- FALSE
      } else {
        ret_list[[factor]] <- TRUE
      }
    } ## End for loop
  } else {
    for (factor in factors) {
      matrix_goal <- design[, goal]
      matrix_factor <- design[, factor]
      matrix_all_formula <- as.formula(glue("~ 0 + {goal} + {factor}"))
      matrix_test <- model.matrix(matrix_all_formula, data = design)
      num_columns <- ncol(matrix_test)
      matrix_decomp <- qr(matrix_test)
      message("The model of ", goal, " and ", factor, " has ", num_columns,
              " and rank ", matrix_decomp[["rank"]])
      if (matrix_decomp[["rank"]] < num_columns) {
        message("This will not work, a different factor should be used.")
        ret_list[[factor]] <- FALSE
      } else {
        ret_list[[factor]] <- TRUE
      }
    }
  }
  return(ret_list)
}

sanitize_model <- function(model, keep_underscores = FALSE,
                           exclude_strings = c("condition", "batch")) {
    startnames <- colnames(model)
    ## Get rid of punctuation in the model column names
    ## optionally leave underscores alone.
    if (isTRUE(keep_underscore)) {
      startnames <- gsub(pattern = "[^_[:^punct:]]", replacement = "", x = startnames, perl = TRUE)
    } else {
      startnames <- gsub(pattern = "[[:punct:]]", replacement = "", x = startnames)
    }
    ## The next lines ensure that conditions/batches which are all numeric will
    ## not cause weird errors for contrasts. Ergo, if a condition is something
    ## like '111', now it will be 'c111' Similarly, a batch '01' will be 'b01'

    ## I think the following two lines are mistaken, I wanted to change the
    ## elements of the model, not the colnames when they are just numeric.
    ##startnames <- gsub(pattern = "^condition(\\d+)$", replacement = "c\\1", x = startnames)
    ##startnames <- gsub(pattern = "^batch(\\d+)$", replacement = "b\\1", x = startnames)
    count <- 1
    for (exclude in exclude_strings) {
      numeric_regex <- paste0("^", exclude, "(\\d+)$")
      startchar <- substr(exclude, 1, 1)
      replace_string <- paste0(startchar, "\\1")
      startnames <- gsub(pattern = numeric_regex, replacement = replace_string, x = startnames)
      regex_string <- paste0("^", exclude)
      startnames <- gsub(pattern = regex_string, replacement = "", x = startnames)
      count <- count + 1
    }
    return(model)
  }

## EOF
