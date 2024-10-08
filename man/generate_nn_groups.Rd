% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify.R
\name{generate_nn_groups}
\alias{generate_nn_groups}
\title{Given an n-dimensional matrix, try some KNN-esque clustering on it.}
\usage{
generate_nn_groups(
  mtrx,
  resolution = 1,
  k = 10,
  type = "snn",
  full = TRUE,
  merge_to = NULL,
  ...
)
}
\arguments{
\item{mtrx}{Matrix to cluster, usually 2d from a point plot.}

\item{resolution}{Used after cluster generation for making neighbor
groups.}

\item{k}{Used during cluster generation.}

\item{type}{Define the type of clustering to perform, currently
only KNN/SNN}

\item{full}{Get the full set of metrics from bluster.}

\item{merge_to}{Use the neighborhood collapse function to set a
hard ceiling on the number of clusters in the final result.}

\item{...}{Extra args for bluster.}
}
\value{
List containing the resulting groups and some information about them.
}
\description{
I want some functions to help me understand clustering.  This is a
first pass at that goal.
}
