## model_testing.r: Some functions to get ready to pass models to various
## DE/test/etc methods.  These functions seek to catch some corner cases when
## playing with the various model types when using DESeq/etc.

extract_stepwise_regression <- function(mtrx, query = "condition",
                                        factors = NULL,
                                        excel = "excel/tc_regression_table.xlsx") {
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
  initial_lm <- lm(as.formula(initial_fstring), data = mtrx)
  initial_summary <- summary(initial_lm) %>% tidy(conf.int = TRUE)
  stepwise_result <- step(initial_lm)
  written <- write_xlsx(data = initial_summary, excel = excel)
  forest_df <- initial_summary[2:nrow(initial_summary), ]
  forest <- ggplot(forest_df, aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
    geom_pointrange(color = "black", size = 0.5) +
    geom_hline(yintercept = 0, color = "steelblue") +
    coord_flip() +
    xlab("") +
    ylab("Coefficient (95% confidence interval)") +
    labs(title = glue("Stepwise Linear Regression Models Estimating \n Effects on {query}")) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12))
  retlist <- list(
    "initial_lm" = initial_lm,
    "initial_summary" = initial_summary,
    "stepwise_result" = stepwise_result,
    "forest" = forest)
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

## EOF
