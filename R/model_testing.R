## model_testing.r: Some functions to get ready to pass models to various
## DE/test/etc methods.  These functions seek to catch some corner cases when
## playing with the various model types when using DESeq/etc.

#' Gather the models and perform forest plots to look at various regression analyses.
#'
#' @param meta Experimental design
#' @param query Factor to query against
#' @param multivariable If set to FALSE, this will iterate over every factor individually.
#' @param intercept Set an intercept? (unused)
#' @param factors Set of factors to query
#' @param excel output xlsx file.
#' @export
extract_linear_regression <- function(meta, query = "condition", multivariable = TRUE,
                                      scale = TRUE, intercept = FALSE, factors = NULL,
                                      excel = NULL) {
  if (isFALSE(multivariable)) {
    serial <- iterate_linear_regression(design, query = query, factors = factors, family = family,
                                        conf = conf, excel = excel)
    return(serial)
  }
  initial_fstring <- glue("{query} ~ .")
  if (isTRUE(scale)) {
    initial_fstring <- glue("scale({query}) ~ .")
  }
  if (!is.null(factors)) {
    initial_fstring <- glue("{query} ~ ")
    if (isTRUE(scale)) {
      initial_fstring <- glue("scale({query}) ~ ")
    }
    for (fct in factors) {
      message("Adding: ", fct)
      if (isTRUE(scale)) {
        initial_fstring <- glue("{initial_fstring} scale({fct}) +")
      } else {
        initial_fstring <- glue("{initial_fstring} {fct} +")
      }
    }
    initial_fstring <- gsub(x = initial_fstring, pattern = " \\+$", replacement = "")
  }
  mesg("Testing regression coefficients with model string: ", initial_fstring, ".")

  ## Sanitize the meta and warn the user
  for (f in colnames(meta)) {
    column_class <- class(meta[[f]])[1]
    if (column_class != "numeric") {
      warning("Coercing the column ", f, " to numeric in order to perform lm.")
      meta[[f]] <- as.numeric(meta[[f]])
    }
  }

  initial_lm <- lm(as.formula(initial_fstring), data = meta)
  ## In the container, this seems to fail with 'no applicable method
  ## for tidy applied to object of class summary.lm
  initial_summary <- summary(initial_lm) %>%
    generics::tidy(conf.int = TRUE)
  ## FIXME: Figure this out, presumably I am feeding lm a model which is not full rank in some way?
  stepwise_result <- try(step(initial_lm), silent = TRUE)
  forest_df <- initial_summary[2:nrow(initial_summary), ]
  colnames(forest_df) <- c("term", "estimate", "std_error", "z", "pr_z", "conf_low", "conf_high")
  forest <- plot_forest_from_regression(forest_df, iterate = FALSE,
                                        type = "linear", intercept = intercept, query = query)
  written <- NULL
  if (!is.null(excel)) {
    xlsx <- init_xlsx(excel)
    wb <- xlsx[["wb"]]
    excel_basename <- xlsx[["basename"]]
    written <- write_xlsx(data = initial_summary, wb = wb, sheet = "lm_summary")
    new_column <- written[["end_col"]] + 2
    try_result <- xlsx_insert_png(
      a_plot = forest, wb = wb, start_col = new_column, sheet = "lm_summary")
    written <- write_xlsx(data = meta, wb = wb, sheet = "input")
    excel_ret <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  }
  retlist <- list(
    "initial_lm" = initial_lm,
    "summary" = initial_summary,
    "stepwise_result" = stepwise_result,
    "forest" = forest,
    "excel" = written)
  return(retlist)
}

#' Invoke what I think is an appropriate logistical regression model.
#'
#' The current implementation uses lm() and assumes everything is a
#' linear model, this will attempt to invoke an appropriate logistic
#' model via glm() and provide similar/identical tables/plots.
#'
#' A reference to myself regarding families:
#'   gaussian: identity, log, inverse
#'   binomial: logit, probit, cauchit
#'   Gamma: inverse, identity, log ??  ooo Gamma.log etc
#'   quasi: logit, probit, cloglog, identity, inverse, log 1/mu^2,
#'   sqrt
#'
#'  For the purposes of a 'normal' logistic regression, I think
#' 'binomial' is sufficient.
#'
#' @param design Experimental design, I need to change this it is not a
#'  matrix.
#' @param query Response variable.
#' @param multivariable When not true, this will iterate over every factor individually.
#' @param factors set of factors to query against the query.
#' @param family The family passed to glm.
#' @param conf Confidence interval chosen for plotting.
#' @param excel Output xlsx file to which we print the f values etc.
#' @export
extract_logistic_regression <- function(design, query = "condition", multivariable = TRUE,
                                        factors = NULL, family = "binomial", conf = 0.95,
                                        scale = TRUE, excel = NULL, intercept = FALSE) {
  if (isFALSE(multivariable)) {
    serial <- iterate_logistic_regression(design, query = query, factors = factors, family = family,
                                          conf = conf, excel = excel)
    return(serial)
  }

  ## TODO: Make this smart enough to accept a formula string in addition to a
  ## vector of factors of interest.
  initial_fstring <- glue("{query} ~ .")
  if (isTRUE(scale) && class(design[[query]])[1] == "numeric") {
    initial_fstring <- glue("scale({query}) ~ .")
  }
  if (!is.null(factors)) {
    initial_fstring <- glue("{query} ~ ")
    if (isTRUE(scale) && class(design[[query]])[1] == "numeric") {
      initial_fstring <- glue("scale({query}) ~ ")
    }
    for (fct in factors) {
      message("Adding: ", fct)
      if (isTRUE(scale) && class(design[[fct]])[1] == "numeric") {
        initial_fstring <- glue("{initial_fstring} scale({fct}) +")
      } else {
        initial_fstring <- glue("{initial_fstring} {fct} +")
      }
    }
    initial_fstring <- gsub(x = initial_fstring, pattern = " \\+$", replacement = "")
  }

  mesg("Testing regression coefficients with model string: ", initial_fstring, ".")
  initial_glm <- glm(as.formula(initial_fstring), data = design,
                     family = family)
  initial_summary <- summary(initial_glm)
  summary_df <- initial_summary[["coefficients"]]
  initial_conf <- stats::confint(initial_glm, level = conf)
  colnames(initial_conf) <- c("conf.low", "conf.high")
  full_summary <- as.data.frame(cbind(summary_df, initial_conf))
  full_summary[["term"]] <- rownames(full_summary)
  percent <- conf * 100
  complete_idx <- complete.cases(full_summary)
  if (sum(!complete_idx) > 0) {
    plot_df <- full_summary[complete_idx, ]
    dropped <- full_summary[!complete_idx, ]
    warning("Dropped the row(s): ", dropped, " from plotting.")
  } else {
    plot_df <- full_summary
  }
  colnames(plot_df) <- c("estimate", "std_error", "z", "pr_z", "conf_low", "conf_high", "term")
  forest <- plot_forest_from_regression(plot_df, percent = percent, intercept = intercept,
                                        family = family, iterate = FALSE, query = query)
  written <- NULL
  if (!is.null(excel)) {
    xlsx <- init_xlsx(excel)
    wb <- xlsx[["wb"]]
    excel_basename <- xlsx[["basename"]]
    written <- write_xlsx(data = full_summary, wb = wb, sheet = "logistic_summary")
    new_column <- written[["end_col"]] + 2
    try_result <- xlsx_insert_png(
      a_plot = forest, wb = wb, start_col = new_column, sheet = "logistic_summary")
    written <- write_xlsx(data = design, wb = wb, sheet = "input")
    excel_ret <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  }
  retlist <- list(
    "initial_glm" = initial_glm,
    "summary" = full_summary,
    "forest" = forest,
    "excel" = written)
  return(retlist)
}

#' Use a vector of factors and design to count up the levels.
get_degrees <- function(design, fctrs) {
  degrees <- 0
  for (f in fctrs) {
    if (is.null(design[[f]])) {
      message("Factor ", f, " is not in the experimental design.")
      next
    }
    fctr <- as.factor(design[[f]])
    range <- table(fctr)
    adder <- length(levels(fctr))
    mesg("Factor: ", f, " has ", adder, " elements.")
    mesg("They range from ", min(range), " to ", max(range), " elements/level.")
    degrees <- degrees + adder
  }
  retlist <- list(
    "sum" = degrees,
    "dof" = degrees - length(fctrs))
  return(retlist)
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
    first_mix <- strsplit(x = mixed_split[1], split = "[[:punct:]|[:space:]]+")[[1]]
    mixers[1] <- first_mix[length(first_mix)]
    mixed_split[2] <- gsub(x = mixed_split[2], pattern = "^[[:space:]]+", replacement = "")
    second_mix <- strsplit(x = mixed_split[2], split = "[[:punct:]|[:space:]]")[[1]]
    mixers[2] <- second_mix[1]
    retlist[["mixers"]] <- mixers
  }

  fct_vector <- c()
  for (i in seq_len(length(factor_vector))) {
    factor_vector[i] <- gsub(x = factor_vector[i], pattern = "[[:space:]]", replacement = "")
    if (grepl(x = factor_vector[i], pattern = "^[[:digit:]]+$")) {
      next
    }
    if (grepl(x = factor_vector[i], pattern = "^[[:alnum:]]+$")) {
      fct_vector <- c(fct_vector, factor_vector[i])
    }
  }
  ## Check for an cell-means model intercept: ~ 0 + f1 + f2
  ##                                            ^
  if (grepl(x = factor_vector[2], pattern = "^[[:digit:]]+$")) {
    retlist[["type"]] <- "cellmeans"
    retlist[["cellmeans_intercept"]] <- factor_vector[2]
  }

  retlist[["factors"]] <- unique(fct_vector)
  return(retlist)
}
setGeneric("get_formula_factors")

#' Perform a series of single regression analyses and tabulate/plot the results.
#'
#' @param design Experimental design.
#' @param query Factor of primary interest.
#' @param factors Set of factors to query against (if not set, then query will be
#'  the first design column, and these will be 2:end)
#' @param conf Choose the confidence interval.
#' @param excel Write the results to this file.
#' @export
iterate_linear_regression <- function(design, query = "condition", factors = NULL,
                                      conf = 0.95, excel = NULL) {
  if (is.null(factors)) {
    all_factors <- colnames(design)
    message("No factors of interest were provided,
using all and assuming the first column (", all_factors[1], ") is the query.")
    query <- all_factors[1]
    factors <- all_factors[2:length(all_factors)]
  }

  summary_df <- data.frame(row.names = "intercept")
  summary_df[["estimate"]] <- 0
  summary_df[["std_error"]] <- 0
  summary_df[["z"]] <- 0
  summary_df[["pr_z"]] <- 0
  summary_df[["conf_low"]] <- 0
  summary_df[["conf_high"]] <- 0
  summary_df[["term"]] <- "intercept"
  summary_colnames <- colnames(summary_df)

  initial_fstring <- glue("{query} ~ ")
  for (fct in factors) {
    message("Adding: ", fct)
    test_fstring <- glue("{initial_fstring} {fct}")
    mesg("Testing regression coefficients with model string: ", test_fstring, ".")
    num_levels <- length(levels(as.factor(design[[fct]])))
    if (num_levels < 2) {
      message("Cannot use ", fct, " it has only 1 level.")
      next
    }
    initial_lm <- lm(as.formula(test_fstring), data = design)
    initial_summary <- summary(initial_lm)
    ## Note that I am explicitly removing the first row, which is the intercept.
    ## In the case where a factor has two elements, this will result in a vector
    ## In the case where a factor has > 2 elements, this results in a df of n-1 rows...
    summary_names <- rownames(initial_summary[["coefficients"]])[-1]
    test_summary_df <- initial_summary[["coefficients"]][-1, ]
    ## Repeat intercept removal
    initial_conf <- as.data.frame(stats::confint(initial_lm, level = conf))
    colnames(initial_conf) <- c("conf_low", "conf_high")
    test_confidence <- initial_conf[-1, ]
    if (class(test_summary_df)[1] == "numeric") {
      test_summary_df <- t.data.frame(as.data.frame(test_summary_df))
      full_summary <- cbind.data.frame(as.data.frame(test_summary_df), test_confidence)
      full_summary[["term"]] <- summary_names
      rownames(full_summary) <- summary_names
      colnames(full_summary) <- summary_colnames
    } else {
      full_summary <- cbind.data.frame(test_summary_df, test_confidence)
      full_summary[["term"]] <- rownames(full_summary)
      colnames(full_summary) <- summary_colnames
    }
    summary_df <- rbind.data.frame(summary_df, full_summary)
  }
  ## Drop the intercept, which I only put in as a place keeper
  summary_df <- summary_df[-1, ]

  complete_idx <- complete.cases(summary_df)
  if (sum(!complete_idx) > 0) {
    plot_df <- summary_df[complete_idx, ]
    dropped <- summary_df[!complete_idx, ]
    warning("Dropped the row(s): ", dropped, " from plotting.")
  } else {
    plot_df <- summary_df
  }
  percent <- conf * 100
  forest <- plot_forest_from_regression(plot_df, percent = percent,
                                        type = "linear", iterate = TRUE, query = query)
  written <- NULL
  if (!is.null(excel)) {
    xlsx <- init_xlsx(excel)
    wb <- xlsx[["wb"]]
    excel_basename <- xlsx[["basename"]]
    written <- write_xlsx(data = full_summary, wb = wb, sheet = "logistic_summary")
    new_column <- written[["end_col"]] + 2
    try_result <- xlsx_insert_png(
      a_plot = forest, wb = wb, start_col = new_column, sheet = "logistic_summary")
    written <- write_xlsx(data = design, wb = wb, sheet = "input")
    excel_ret <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  }
  retlist <- list(
    "initial_lm" = initial_lm,
    "summary" = full_summary,
    "forest" = forest,
    "excel" = written)
  return(retlist)
}

#' Perform a series of single regression analyses and tabulate/plot the results.
#'
#' @param design Experimental design.
#' @param query Factor of primary interest.
#' @param factors Set of factors to query against (if not set, then query will be
#'  the first design column, and these will be 2:end)
#' @param conf Choose the confidence interval.
#' @param excel Write the results to this file.
#' @export
iterate_logistic_regression <- function(design, query = "condition",
                                        factors = NULL, family = "binomial", conf = 0.95,
                                        excel = NULL) {
  if (is.null(factors)) {
    all_factors <- colnames(design)
    message("No factors of interest were provided,
using all and assuming the first column (", all_factors[1], ") is the query.")
    query <- all_factors[1]
    factors <- all_factors[2:length(all_factors)]
  }

  summary_df <- data.frame(row.names = "intercept")
  summary_df[["estimate"]] <- 0
  summary_df[["std_error"]] <- 0
  summary_df[["z"]] <- 0
  summary_df[["pr_z"]] <- 0
  summary_df[["conf_low"]] <- 0
  summary_df[["conf_high"]] <- 0
  summary_df[["term"]] <- "intercept"
  summary_colnames <- colnames(summary_df)

  initial_fstring <- glue("{query} ~ ")
  for (fct in factors) {
    message("Adding: ", fct)
    test_fstring <- glue("{initial_fstring} {fct}")
    mesg("Testing regression coefficients with model string: ", test_fstring, ".")
    num_levels <- length(levels(as.factor(design[[fct]])))
    if (num_levels < 2) {
      message("Cannot use ", fct, " it has only 1 level.")
      next
    }
    initial_glm <- glm(as.formula(test_fstring), data = design,
                       family = family)
    initial_summary <- summary(initial_glm)
    ## Note that I am explicitly removing the first row, which is the intercept.
    ## In the case where a factor has two elements, this will result in a vector
    ## In the case where a factor has > 2 elements, this results in a df of n-1 rows...
    summary_names <- rownames(initial_summary[["coefficients"]])[-1]
    test_summary_df <- initial_summary[["coefficients"]][-1, ]
    ## Repeat intercept removal
    initial_conf <- as.data.frame(stats::confint(initial_glm, level = conf))
    colnames(initial_conf) <- c("conf_low", "conf_high")
    test_confidence <- initial_conf[-1, ]
    if (class(test_summary_df)[1] == "numeric") {
      test_summary_df <- t.data.frame(as.data.frame(test_summary_df))
      full_summary <- cbind.data.frame(as.data.frame(test_summary_df), test_confidence)
      full_summary[["term"]] <- summary_names
      rownames(full_summary) <- summary_names
      colnames(full_summary) <- summary_colnames
    } else {
      full_summary <- cbind.data.frame(test_summary_df, test_confidence)
      full_summary[["term"]] <- rownames(full_summary)
      colnames(full_summary) <- summary_colnames
    }
    summary_df <- rbind.data.frame(summary_df, full_summary)
  }
  ## Drop the intercept, which I only put in as a place keeper
  summary_df <- summary_df[-1, ]

  complete_idx <- complete.cases(summary_df)
  if (sum(!complete_idx) > 0) {
    plot_df <- summary_df[complete_idx, ]
    dropped <- summary_df[!complete_idx, ]
    warning("Dropped the row(s): ", dropped, " from plotting.")
  } else {
    plot_df <- summary_df
  }
  percent <- conf * 100
  forest <- plot_forest_from_regression(plot_df, percent = percent, iterate = TRUE, query = query)
  written <- NULL
  if (!is.null(excel)) {
    xlsx <- init_xlsx(excel)
    wb <- xlsx[["wb"]]
    excel_basename <- xlsx[["basename"]]
    written <- write_xlsx(data = summary_df, wb = wb, sheet = "logistic_summary")
    new_column <- written[["end_col"]] + 2
    try_result <- xlsx_insert_png(
      a_plot = forest, wb = wb, start_col = new_column, sheet = "logistic_summary")
    written <- write_xlsx(data = design, wb = wb, sheet = "input")
    excel_ret <- try(openxlsx::saveWorkbook(wb, excel, overwrite = TRUE))
  }
  retlist <- list(
    "initial_glm" = initial_glm,
    "summary" = summary_df,
    "forest" = forest,
    "excel" = written)
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

#' Ensure an experimental model is safe to use
#'
#' @param model input model
#' @param keep_underscore I previously dropped all punctuation including underscores.
#' @param exclude_strings simplify the strings for these factors a little.
sanitize_model <- function(model, keep_underscore = FALSE,
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
test_model_rank <- function(design, goal = "condition", factors = NULL, ...) {
  arglist <- list(...)
  ## For testing, use some existing matrices/data
  message("There are ", length(levels(as.factor(design[, goal]))),
          " levels in the goal: ", goal, ".")
  retlist <- list()
  if (is.null(factors)) {
    message("Testing an experimental design with only ", goal, ".")
    matrix_all_formula <- as.formula(glue("~ {goal}"))
    matrix_test <- model.matrix(matrix_all_formula, data = design)
    num_columns <- ncol(matrix_test)
    matrix_decomp <- qr(matrix_test)
    message("The model of ", goal, " has ", num_columns,
            " levels and rank ", matrix_decomp[["rank"]], ".")
    if (matrix_decomp[["rank"]] < num_columns) {
      message("This will not work, a different factor should be used.")
      retlist[[goal]] <- FALSE
    } else {
      retlist[[goal]] <- TRUE
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
      matrix_all_formula <- as.formula(glue("~ {goal} + {factor}"))
      matrix_test <- model.matrix(matrix_all_formula, data = design)
      num_columns <- ncol(matrix_test)
      matrix_decomp <- qr(matrix_test)
      message("The model of ", goal, " and ", factor, " has ", num_columns,
              " and rank ", matrix_decomp[["rank"]])
      if (matrix_decomp[["rank"]] < num_columns) {
        message("This will not work, a different factor should be used.")
        retlist[[factor]] <- FALSE
      } else {
        retlist[[factor]] <- TRUE
      }
    } ## End for loop
  } else {
    for (factor in factors) {
      matrix_goal <- design[, goal]
      matrix_factor <- design[, factor]
      matrix_all_formula <- as.formula(glue("~ {goal} + {factor}"))
      matrix_test <- model.matrix(matrix_all_formula, data = design)
      num_columns <- ncol(matrix_test)
      matrix_decomp <- qr(matrix_test)
      message("The model of ", goal, " and ", factor, " has ", num_columns,
              " and rank ", matrix_decomp[["rank"]])
      if (matrix_decomp[["rank"]] < num_columns) {
        message("This will not work, a different factor should be used.")
        retlist[[factor]] <- FALSE
      } else {
        retlist[[factor]] <- TRUE
      }
    }
  }
  return(retlist)
}

test_expt_model_rank <- function(expt, fstring = "~ donor") {
  retlist <- list()
  design <- pData(expt)
  fctrs <- get_formula_factors(fstring)
  degrees <- get_degrees(design, fctrs[["factors"]])
  test_formula <- as.formula(fstring)
  matrix_test <- model.matrix(test_formula, data = design)
  num_columns <- ncol(matrix_test)
  matrix_decomp <- qr(matrix_test)
  message("The model of ", fstring, " has ", num_columns,
          " column
and rank ", matrix_decomp[["rank"]])
  if (matrix_decomp[["rank"]] < num_columns) {
    message("This will not work, a different factor should be used.")
    retlist[["full_rank"]] <- FALSE
  } else {
    retlist[["full_rank"]] <- TRUE
  }
  retlist[["matrix"]] <- matrix_test
  retlist[["matrix_decomp"]] <- matrix_decomp
  retlist[["matrix_num_column"]] <- num_columns
  return(retlist)
}

#' Given the result from one of the regression testers, plot it!
#'
#' @param df The primary dataframe from one of the sister regression functions above.
#' @param percent Confidence interval chosen.
#' @param type Either linear or logistic.
#' @param iterate Was this a series of single-variable regressions, or all in one?
#' @param family Only currently used for logistic.
#' @param intercept Include the intercept in the plot?
#' @export
plot_forest_from_regression <- function(plot_df, percent = 95, type = "logistic",
                                        iterate = TRUE, family = "binomial",
                                        intercept = FALSE, base_size = 18, title_size = 22,
                                        axis_size = 20, query = "condition") {
  ## On may not wish to see the intercept in the plot, if not, remove it here...
  if (!isTRUE(intercept)) {
    mesg("Removing intercept row(s).")
    intercept_idx <- grepl(x = rownames(plot_df), pattern = "ntercept")
    plot_df <- plot_df[!intercept_idx, ]
  }
  iterate_string <- "Iterative Regression Models"
  if (!isTRUE(iterate)) {
    iterate_string <- "Combined Regression Model"
  }
  title <- glue("Logistic ({family}) {iterate_string}\n Estimating Effects on {query}")
  if (type != "logistic") {
    title <- glue("Linear {iterate_string}\n Estimating Effects on {query}")
  }
  ylabel <- glue("Coefficient {percent}% confidence interval")
  forest <- ggplot(plot_df, aes(x = term, y = estimate, ymin = conf_low, ymax = conf_high)) +
    ggplot2::geom_pointrange(color = "black", size = 0.5) +
    ## scale_y_continuous(limits = c(min_val, max_val)) +
    ggplot2::geom_hline(yintercept = 0, color = "steelblue") +
    ggplot2::coord_flip() +
    ggplot2::xlab("") +
    ggplot2::ylab(ylabel) +
    ggplot2::labs(title = title) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, face = "bold"),
      axis.text.x = ggplot2::element_text(size = axis_size),
      axis.text.y = ggplot2::element_text(size = axis_size))
  return(forest)
}
## EOF
