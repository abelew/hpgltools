% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_plots.R
\name{get_plot_columns}
\alias{get_plot_columns}
\title{A small rat's nest of if statements intended to figure out what columns
are wanted to plot a MA/Volcano from any one of a diverse set of possible
input types.}
\usage{
get_plot_columns(data, type, p_type = "adj", adjp = TRUE)
}
\arguments{
\item{data}{Data structure in which to hunt columns/data.}

\item{type}{Type of method used to make the data.}

\item{p_type}{Use adjusted p-values?}

\item{adjp}{I think this is reundant.}
}
\description{
I split this function away from the main body of extract_de_plots()
so that I can come back to it and strip it down to something a bit
more legible.  Eventually I want to dispatch this logic off to
separate functions depending on the class of the input.
}
\details{
This function should die in a fire.
}
