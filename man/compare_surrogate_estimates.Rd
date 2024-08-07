% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normalize_batch.R
\name{compare_surrogate_estimates}
\alias{compare_surrogate_estimates}
\title{Perform a comparison of the surrogate estimators demonstrated by Jeff Leek.}
\usage{
compare_surrogate_estimates(
  expt,
  extra_factors = NULL,
  filter_it = TRUE,
  filter_type = TRUE,
  do_catplots = FALSE,
  surrogates = "be",
  ...
)
}
\arguments{
\item{expt}{Experiment containing a design and other information.}

\item{extra_factors}{Character list of extra factors which may be included in
the final plot of the data.}

\item{filter_it}{Most of the time these surrogate methods get mad if there
are 0s in the data.  Filter it?}

\item{filter_type}{Type of filter to use when filtering the input data.}

\item{do_catplots}{Include the catplots?  They don't make a lot of sense yet,
so probably no.}

\item{surrogates}{Use 'be' or 'leek' surrogate estimates, or choose a
number.}

\item{...}{Extra arguments when filtering.}
}
\value{
List of the results.
}
\description{
This is entirely derivative, but seeks to provide similar estimates for one's
own actual data and catch corner cases not taken into account in that
document (for example if the estimators don't converge on a surrogate
variable). This will attempt each of the surrogate estimators described by
Leek: pca, sva supervised, sva unsupervised, ruv supervised, ruv residuals,
ruv empirical. Upon completion it will perform the same limma expression
analysis and plot the ranked t statistics as well as a correlation plot
making use of the extracted estimators against condition/batch/whatever
else. Finally, it does the same ranking plot against a linear fitting Leek
performed and returns the whole pile of information as a list.
}
\seealso{
[normalize_expt()] [plot_pca()] [all_adjuster()] [corrplot] [ffpe]
}
