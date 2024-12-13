% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expt.R
\name{set_expt_batches}
\alias{set_expt_batches}
\title{Change the batches of an expt.}
\usage{
set_expt_batches(expt, fact, ids = NULL, ...)
}
\arguments{
\item{expt}{Expt to modify.}

\item{fact}{Batches to replace using this factor.}

\item{ids}{Specific samples to change.}

\item{...}{Extra options are like spinach.}
}
\value{
The original expt with some new metadata.
}
\description{
When exploring differential analyses, it might be useful to play with the
conditions/batches of the experiment.  Use this to make that easier.
}
\examples{
\dontrun{
 expt = set_expt_batches(big_expt, factor = c(some,stuff,here))
}
}
\seealso{
[create_expt()] [set_expt_conditions()] [Biobase]
}
