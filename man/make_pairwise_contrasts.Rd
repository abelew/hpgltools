% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_shared.R
\name{make_pairwise_contrasts}
\alias{make_pairwise_contrasts}
\title{Run makeContrasts() with all pairwise comparisons.}
\usage{
make_pairwise_contrasts(
  model,
  conditions,
  do_identities = FALSE,
  do_extras = TRUE,
  do_pairwise = TRUE,
  keepers = NULL,
  extra_contrasts = NULL,
  keep_underscore = FALSE,
  ...
)
}
\arguments{
\item{model}{Describe the conditions/batches/etc in the experiment.}

\item{conditions}{Factor of conditions in the experiment.}

\item{do_identities}{Include all the identity strings? Limma can
use this information while edgeR can not.}

\item{do_extras}{Include extra contrasts?  This seems redundant with extra_contrasts
below, but there is a reason for it.}

\item{do_pairwise}{Include all pairwise strings? This shouldn't
need to be set to FALSE, but just in case.}

\item{keepers}{Only extract this subset of all possible pairwise contrasts.}

\item{extra_contrasts}{Optional string of extra contrasts to include.}

\item{...}{Extra arguments passed here are caught by arglist.}
}
\value{
List including the following information:
\enumerate{
 \item all_pairwise_contrasts = the result from makeContrasts(...)
 \item identities = the string identifying each condition alone
 \item all_pairwise = the string identifying each pairwise comparison alone
 \item contrast_string = the string passed to R to call makeContrasts(...)
 \item names = the names given to the identities/contrasts
}
}
\description{
In order to have uniformly consistent pairwise contrasts, I decided
to avoid potential human erors(sic) by having a function generate
all contrasts.
}
\details{
Invoked by the _pairwise() functions.
}
\examples{
\dontrun{
 pretend <- make_pairwise_contrasts(model, conditions)
}
}
\seealso{
[limma::makeContrasts()]
}
