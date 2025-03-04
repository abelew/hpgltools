% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize_shared.R
\name{do_batch}
\alias{do_batch}
\title{Actually runs the batch method, this more than anything shows that hpgl_norm is too complicated.}
\usage{
do_batch(
  count_table,
  method = "raw",
  design,
  batch1 = "batch",
  current_state = NULL,
  current_design = NULL,
  expt_state = NULL,
  surrogate_method = "be",
  surrogates = NULL,
  low_to_zero = FALSE,
  cpus = 4,
  batch2 = NULL,
  noscale = TRUE,
  adjust_method = "ruv",
  batch_step = 4
)
}
\arguments{
\item{count_table}{The counts in their current state.}

\item{method}{Batch/SV method to employ.}

\item{design}{Experimental design, requiring columns named 'condition' and 'batch'.}

\item{current_state}{State of the data before messing with it.}

\item{adjust_method}{Method to use to modify the counts after finding the surrogates.}

\item{batch_step}{Choose when to perform this in the set of normalization tasks.}

\item{...}{Extra arguments passed to sva and friends.}
}
\description{
Actually runs the batch method, this more than anything shows that hpgl_norm is too complicated.
}
