#' Annotation sanitizers
#'
#' @export
setGeneric("sanitize_annotations", function(input, ...) {
  message("This function is intended to sanitize annotations for a data structure.")
  message("It was passed an object of type ", class(input),
          " and does not know what to do with it.")
  return(NULL)
  standardGeneric("sanitize_annotations")
})

#' Given an expressionset, sanitize the gene information data.
#'
#' @param expt Input expressionset.
#' @param columns Set of columns to sanitize, otherwise all of them.
#' @param na_value Fill in NA with this.
#' @param lower sanitize capitalization.
#' @param punct Remove punctuation?
#' @param factorize Convert columns to factors?  When set to 'heuristic'
#'  this tries out as.factor and sees if the number of levels is silly.
#' @param max_levels The definition of 'silly' above.
#' @param spaces Allow spaces in the data?
#' @param numbers Sanitize number formats (e.g. 1.000.000,0 vs. 1,000,000.0)
#' @param numeric Set columns to numeric when possible?
#' @export
setMethod(
  "sanitize_annotations", signature = signature(input = "expt"),
  definition =  function(input, columns = NULL, na_value = "notapplicable",
                         lower = TRUE, punct = TRUE, factorize = "heuristic",
                         max_levels = NULL, spaces = FALSE, numbers = NULL,
                         numeric = FALSE) {
    expt <- input
    meta <- fData(expt)
    sanitized <- sanitize_metadata(meta, columns = columns, na_value = na_value,
                                   lower = lower, punct = punct, factorize = factorize,
                                   max_levels = max_levels, spaces = spaces,
                                   numbers = numbers, numeric = numeric)
    fData(expt) <- sanitized
    return(expt)
  })

#' @export
setMethod(
  "sanitize_annotations", signature = signature(input = "SummarizedExperiment"),
  definition =  function(input, columns = NULL, na_value = "notapplicable",
                         lower = TRUE, punct = TRUE, factorize = "heuristic",
                         max_levels = NULL, spaces = FALSE, numbers = NULL,
                         numeric = FALSE) {
    se <- input
    meta <- rowData(se)
    sanitized <- sanitize_metadata(meta, columns = columns, na_value = na_value,
                                   lower = lower, punct = punct, factorize = factorize,
                                   max_levels = max_levels, spaces = spaces,
                                   numbers = numbers, numeric = numeric)
    rowData(se) <- sanitized
    return(se)
  })

#' Metadata sanitizers
#'
#' @export
sanitize_metadata <- function(input, ...) {
  message("This function is intended to sanitize metadata for a data structure.")
  message("It was passed an object of type ", class(input),
          " and does not know what to do with it.")
  return(NULL)
  standardGeneric("sanitize_metadata")
}
setGeneric("sanitize_metadata")

#' Given an expressionset, sanitize pData columns of interest.
#'
#' I wrote this function after spending a couple of hours confused
#' because one cell in my metadata said 'cure ' instead of 'cure' and
#' I could not figure out why chaos reigned in my analyses.  There is
#' a sister to this somewhere else which checks that the expected
#' levels of a metadata factor are consistent; this is because in
#' another analysis we essentially had a cell which said 'cyre' and a
#' similar data explosion occurred.
#'
#' @param meta Input metadata
#' @param columns Set of columns to check, if left NULL, all columns
#'  will be molested.
#' @param na_value Fill NA values with a string.
#' @param lower Set everything to lowercase?
#' @param punct Remove punctuation?
#' @param factorize Set some columns to factors?  If set to a vector
#'  of length >=1, then set all of the provided columns to factors.
#'  When set to 'heuristic', set any columns with <= max_levels
#'  different elements to factors.
#' @param max_levels When heuristically setting factors, use this as
#'  the heuristic, when NULL it is the number of samples / 6
#' @param spaces Remove any spaces in this column?
#' @param numbers Sanitize numbers by adding a prefix character to them?
#' @param numeric Recast the values as numeric when possible?
#' @export
setMethod(
  "sanitize_metadata", signature = signature(input = "data.frame"),
  definition = function(input, columns = NULL, na_value = "notapplicable",
                        lower = TRUE, punct = TRUE, factorize = "heuristic",
                        max_levels = NULL, spaces = FALSE, numbers = NULL,
                        numeric = FALSE) {
    meta <- input
    if (is.null(max_levels)) {
      max_levels <- nrow(meta) / 6.0
    }
    if (is.null(columns)) {
      columns <- colnames(meta)
    }
    for (col in seq_along(columns)) {
      todo <- columns[col]
      mesg("Sanitizing metadata column: ", todo, ".")
      if (! todo %in% colnames(meta)) {
        mesg("The column ", todo, " is missing, skipping it (also warning this).")
        warning("The column ", todo, " is missing, skipping it.")
        next
      }
      ## First get rid of trailing/leading spaces, those anger me and are crazy hard to find
      meta[[todo]] <- gsub(pattern = "^[[:space:]]", replacement = "", x = meta[[todo]])
      meta[[todo]] <- gsub(pattern = "[[:space:]]$", replacement = "", x = meta[[todo]])
      ## Set the column to lowercase, I have recently had a rash of mixed case sample sheet data.
      if (isTRUE(numeric)) {
        if (!is.null(na_value)) {
          na_idx <- is.na(meta[[todo]])
          meta[na_idx, todo] <- na_value
          mesg("Setting numeric NAs to ", na_value, ".")
          ## I assume it is ok to suppress warnings here given my explicit setting of NA values.
          meta[[todo]] <- suppressWarnings(as.numeric(meta[[todo]]))
        }
      } else {
        if (isTRUE(lower)) {
          mesg("Setting everything to lowercase.")
          meta[[todo]] <- tolower(meta[[todo]])
        }
        ## I think punctuation needs to go
        if (isTRUE(punct)) {
          mesg("Removing punctuation.")
          meta[[todo]] <- gsub(pattern = "[^_[:^punct:]]",
                               replacement = "", x = meta[[todo]],
                               perl = TRUE)
        }
        if (!is.null(numbers)) {
          if (isTRUE(numbers)) {
            ## Use the first letter of the column name.
            numbers <- gsub(x = todo, pattern = "^(\\w{1}).*$", replacement = "\\1")
          }
          mesg("Adding a prefix to bare numbers.")
          meta[[todo]] <- gsub(pattern = "^([[:digit:]]+)$",
                               replacement = glue("{numbers}\\1"), x = meta[[todo]])
        }

        if (!is.null(na_value)) {
          na_idx <- is.na(meta[[todo]])
          meta[na_idx, todo] <- na_value
          mesg("Setting NAs to character/factor value ", na_value, ".")
        }
        ## Handle spaces after the previous changes.
        if (isTRUE(spaces)) {
          mesg("Removing all spaces.")
          meta[[todo]] <- gsub(pattern = "[[:space:]]", replacement = "", x = meta[[todo]])
        }
        if (!is.null(factorize) &&
              (length(factorize) == 1 && factorize[1] == "heuristic")) {
          nlevels <- length(levels(as.factor(meta[[todo]])))
          if (nlevels <= max_levels) {
            mesg("Setting column ", todo, " to a factor.")
            meta[[todo]] <- as.factor(meta[[todo]])
          }
        }
      } ## End checking if we are sanitizing numeric or other data.
    } ## End iterating over the columns of interest

    return(meta)
  })

#' Metadata sanitizers for an expt
#' @export
setMethod(
  "sanitize_metadata", signature = signature(input = "expt"),
  definition = function(input, columns = NULL, na_string = "notapplicable",
                        lower = TRUE, punct = TRUE, factorize = "heuristic",
                        max_levels = NULL, spaces = FALSE, numbers = NULL) {
    expt <- input
    old_meta <- pData(expt)
    new_meta <- sanitize_metadata(old_meta, columns = columns, na_string = na_string,
                                  lower = lower, punct = punct,
                                  factorize = factorize, max_levels = max_levels,
                                  spaces = spaces, numbers = numbers)
    pData(expt) <- new_meta
    return(expt)
  })

#' Metadata sanitizers for an expt
#' @export
setMethod(
  "sanitize_metadata", signature = signature(input = "SummarizedExperiment"),
  definition = function(input, columns = NULL, na_value = "notapplicable",
                        lower = TRUE, punct = TRUE, factorize = "heuristic",
                        max_levels = NULL, spaces = FALSE, numbers = NULL,
                        numeric = FALSE) {
    se <- input
    old_meta <- as.data.frame(colData(se))
    new_meta <- sanitize_metadata(old_meta, columns = columns, na_value = na_value,
                                  lower = lower, punct = punct,
                                  factorize = factorize, max_levels = max_levels,
                                  spaces = spaces, numbers = numbers, numeric = numeric)
    colData(se) <- new_meta
    return(se)
  })

#' Metadata sanitizers for an expressionset
#' @export
setMethod(
  "sanitize_metadata", signature = signature(input = "ExpressionSet"),
  definition = function(input, columns = NULL, na_string = "notapplicable",
                        lower = TRUE, punct = TRUE, spaces = FALSE,
                        numbers = NULL) {
    exp <- input
    old_meta <- pData(meta)
    new_meta <- sanitize_metadata(old_meta, columns = columns, na_string = na_string,
                                  lower = lower, punct = punct, spaces = spaces,
                                  numbers = numbers)
    pData(exp) <- new_meta
    return(exp)
  })
