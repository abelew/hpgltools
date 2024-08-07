% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_proteomics.R
\name{plot_cleaved}
\alias{plot_cleaved}
\title{Plot the average mass and expected intensity of a set of sequences given an
enzyme.}
\usage{
plot_cleaved(pep_sequences, enzyme = "trypsin", start = 600, end = 1500)
}
\arguments{
\item{pep_sequences}{Set of protein sequences.}

\item{enzyme}{One of the allowed enzymes for cleaver.}

\item{start}{Limit the set of fragments from this point}

\item{end}{to this point.}
}
\value{
List containing the distribution of weights and the associated plot.
}
\description{
This uses the cleaver package to generate a plot of expected intensities
vs. weight for a list of protein sequences.
}
