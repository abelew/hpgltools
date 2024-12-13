% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tnseq.R
\name{tnseq_multi_saturation}
\alias{tnseq_multi_saturation}
\title{Plot the saturation of multiple libraries simultaneously.}
\usage{
tnseq_multi_saturation(
  meta,
  meta_column,
  ylimit = 100,
  column = "Reads",
  adjust = 1,
  ggstatsplot = FALSE
)
}
\arguments{
\item{meta}{Experimental metadata}

\item{meta_column}{Metadata column containing the filenames to query.}

\item{ylimit}{Maximum y axis}

\item{column}{Data file column to use for density calculation.}

\item{adjust}{Density adjustment.}

\item{ggstatsplot}{Include pretty ggstatsplot plot?}
}
\value{
a plot and table of the saturation for all samples.
}
\description{
Plot the saturation of multiple libraries simultaneously.
}
