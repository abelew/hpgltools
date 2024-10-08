% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_plots.R
\name{plot_de_pvals}
\alias{plot_de_pvals}
\title{Given a DE table with p-values, plot them.}
\usage{
plot_de_pvals(
  combined_data,
  type = "limma",
  p_type = "both",
  columns = NULL,
  ...
)
}
\arguments{
\item{combined_data}{Table to extract the values from.}

\item{type}{If provided, extract the {type}_p and {type}_adjp columns.}

\item{p_type}{Which type of pvalue to show (adjusted, raw, or all)?}

\item{columns}{Otherwise, extract whatever columns are provided.}

\item{...}{Arguments passed through to the histogram plotter}
}
\value{
Multihistogram of the result.
}
\description{
Plot a multi-histogram containing (adjusted)p-values.
}
\details{
The assumption of this plot is that the adjustment will significantly
decrease the representation of genes in the 'highly significant' range of
p-values.  However, it is hoped that it will not utterly remove them.
}
\seealso{
[plot_histogram()]
}
