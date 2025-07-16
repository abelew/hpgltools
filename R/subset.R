#' Subset an expt
#'
#' Working on improving my understanding of how R sets up functions by type.
#' @param expt expt to subset
#' @param i Set of genes to keep
#' @param j Set of samples to keep.
#' @param ... Parameters to pass to subset_genes/subset_expt.
#' @param value New subset value
#' @export
`[.expt` <- function(expt, i, j, ..., value) {
  if (!missing(i)) {
    message("Subsetting on features.")
    expt <- subset_genes(expt, ids = i, method = "keep", ...)
  }
  if (!missing(j)) {
    message("Subsetting on samples.")
    expt <- subset_expt(expt, ids = j, ...)
  }
  if (missing(i) && missing(j)) {
    message("No subset was provided, returning the original.")
  }
  return(expt)
}

#' Simplifying subset on metadata when i is a character.
#'
#' @param x an expt
#' @param i Column to extract
#' @param j not sure
#' @param ... extra arguments.
#' @export
setMethod(
  "[[", signature = c(x = "expt", i = "character"),
  definition = function(x, i, j, ...) {
    pData(x)[[i]]
  })

#' Simplifying subset on metadata.
#'
#' @param x an expt
#' @param i Column to extract
#' @param j not sure
#' @param ... extra args.
#' @export
setMethod(
  "[[", signature = c(x = "expt", i = "ANY"),
  definition = function(x, i, j, ...) {
    pData(x)[[i]]
  })

#' Simplifying subset on metadata when i is anything and j is missing.
#'
#' @param x an expt
#' @param i Column to extract
#' @param j not sure
#' @param ... extra arguments.
#' @export
setMethod(
  "[[", signature = c("expt", "ANY", "missing"),
  definition = function(x, i, j, ...) {
    pData(x)[[i]]
  })

#' Simplifying subset on metadata when i is anything and j is missing.
#'
#' @param x an expt
#' @param i Column to extract
#' @param j not sure
#' @param ... extra arguments.
#' @param value new value
#' @export
setReplaceMethod(
  "[[", c("expt", "ANY", "missing"),
  function(x, i, j, ..., value) {
    pData(x)[[i, ...]] <- value
    return(x)
  })
