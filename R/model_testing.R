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
                                      intercept = TRUE, factors = NULL, excel = NULL) {
  if (isFALSE(multivariable)) {
    serial <- iterate_linear_regression(design, query = query, factors = factors, family = family,
                                        conf = conf, excel = excel)
    return(serial)
  }
  initial_fstring <- glue("{query} ~ .")
  if (!is.null(factors)) {
    initial_fstring <- glue("{query} ~ ")
    for (fct in factors) {
      message("Adding: ", fct)
      initial_fstring <- glue("{initial_fstring} {fct} +")
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
  initial_summary <- summary(initial_lm) %>%
    generics::tidy(conf.int = TRUE)
  ## FIXME: Figure this out, presumably I am feeding lm a model which is not full rank in some way?
  stepwise_result <- try(step(initial_lm), silent = TRUE)
  forest_df <- initial_summary[2:nrow(initial_summary), ]
  colnames(forest_df) <- c("estimate", "std_error", "z", "pr_z", "conf_low", "conf_high", "term")
  forest <- plot_forest_from_regression(plot_df)
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
    "initial_summary" = initial_summary,
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
                                        excel = NULL) {
  if (isFALSE(multivariable)) {
    serial <- iterate_logistic_regression(design, query = query, factors = factors, family = family,
                                          conf = conf, excel = excel)
    return(serial)
  }
  initial_fstring <- glue("{query} ~ .")
  if (!is.null(factors)) {
    initial_fstring <- glue("{query} ~ ")
    for (fct in factors) {
      message("Adding: ", fct)
      initial_fstring <- glue("{initial_fstring} {fct} +")
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
  stepwise_result <- step(initial_glm)
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
  forest <- plot_forest_from_regression(plot_df, percent = percent, family = family)
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
    "stepwise_result" = stepwise_result,
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
  forest <- plot_forest_from_regression(plot_df, percent = percent)

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
    "stepwise_result" = stepwise_result,
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
  forest <- plot_forest_from_regression(plot_df, percent = percent)

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
    "stepwise_result" = stepwise_result,
    "forest" = forest,
    "excel" = written)
  return(retlist)
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
      ret_list[[goal]] <- 0
    } else {
      ret_list[[goal]] <- 1
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
        ret_list[[factor]] <- 0
      } else {
        ret_list[[factor]] <- 1
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
        ret_list[[factor]] <- 0
      } else {
        ret_list[[factor]] <- 1
      }
    }
  }
  return(ret_list)
}

#' Given the result from one of the regression testers, plot it!
#'
#' @param df The primary dataframe from one of the sister regression functions above.
#' @param percent Confidence interval chosen.
#' @param type Either linear or logistic.
#' @param family Only currently used for logistic.
plot_forest_from_regression <- function(df, percent = 95,
                                        type = "logistic", family = "binomial") {
  title <- glue("Logistic ({family}) Regression Models Estimating \n Effects on {query}")
  if (type != "logistic") {
    title <- glue("Linear Regression Models Estimating \n Effects on {query}")
  }
  ylabel <- glue("Coefficient {percent}% confidence interval")
  forest <- ggplot(plot_df, aes(x = term, y = estimate, ymin = conf_low, ymax = conf_high)) +
    geom_pointrange(color = "black", size = 0.5) +
    ## scale_y_continuous(limits = c(min_val, max_val)) +
    geom_hline(yintercept = 0, color = "steelblue") +
    coord_flip() +
    xlab("") +
    ylab(ylabel) +
    labs(title = title) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12))
  return(forest)
}
## EOF
