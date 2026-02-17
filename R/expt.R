## All of these functions have been deprecated for about a year, I am going to start
## culling them.

## Ensure imports are loaded from 01_hpgltools.R
#' @include 01_hpgltools.R
NULL

#' @export
combine_expts <- function(...) {
  combine_se(...)
}

#' @export
create_expt <- function(...) {
  create_se(...)
}

#' @export
make_pombe_expt <- function(...) {
  make_pombe_se(...)
}

#' @export
semantic_expt_filter <- function(...) {
  semantic_filter(...)
}

#' @export
subset_expt <- function(...) {
  subset_se(...)
}

#' @export
write_expt <- function(...) {
  write_se(...)
}

## EOF
