## All of these functions have been deprecated for about a year, I am going to start
## culling them.

## Ensure imports are loaded from 01_hpgltools.R
#' @include 01_hpgltools.R
NULL

#' combine summarized experiments when using deprecated function names.
#'
#' @param ... Arguments passed to combine_se()
#' @export
combine_expts <- function(...) {
  combine_se(...)
}

#' Create a summarized experiment when using deprecated function names.
#'
#' @param ... Arguments passed to create_se()
#' @export
create_expt <- function(...) {
  create_se(...)
}

#' Filter a dataset by gene when using deprecated function names.
#'
#' @param ... Arguments passed to semantic_filter()
#' @export
semantic_expt_filter <- function(...) {
  semantic_filter(...)
}

#' Subset a summarizedExperiment when using deprecated function names.
#'
#' @param ... Arguments passed to subset_se()
#' @export
subset_expt <- function(...) {
  subset_se(...)
}

## EOF
