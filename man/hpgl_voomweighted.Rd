% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_limma.R
\name{hpgl_voomweighted}
\alias{hpgl_voomweighted}
\title{A minor change to limma's voom with quality weights to attempt to address some corner cases.}
\usage{
hpgl_voomweighted(
  data,
  fun_model,
  libsize = NULL,
  normalize.method = "none",
  plot = TRUE,
  span = 0.5,
  var.design = NULL,
  method = "genebygene",
  maxiter = 50,
  tol = 1e-10,
  trace = FALSE,
  replace.weights = TRUE,
  col = NULL,
  ...
)
}
\arguments{
\item{data}{Some data!}

\item{fun_model}{A model for voom() and arrayWeights()}

\item{libsize}{Library sizes passed to voom().}

\item{normalize.method}{Passed to voom()}

\item{plot}{Do the plot of mean variance?}

\item{span}{yes}

\item{var.design}{maybe}

\item{method}{kitty!}

\item{maxiter}{50 is good}

\item{tol}{I have no tolerance.}

\item{trace}{no trace for you.}

\item{replace.weights}{Replace the weights?}

\item{col}{yay columns!}

\item{...}{more arguments!}
}
\value{
a voom return
}
\description{
This copies the logic employed in hpgl_voom().  I suspect one should not use it.
}
\examples{
\dontrun{
## No seriously, dont run this, I think it is wiser to use the functions
## provided by limma. But this provides a place to test stuff out.
 voom_result <- hpgl_voomweighted(dataset, model)
}
}
\seealso{
[limma::voom()]
}
