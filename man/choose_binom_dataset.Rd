% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_shared.R
\name{choose_binom_dataset}
\alias{choose_binom_dataset}
\title{A sanity check that a given set of data is suitable for methods which assume
a negative binomial distribution of input.}
\usage{
choose_binom_dataset(input, force = FALSE, ...)
}
\arguments{
\item{input}{Expressionset containing expt object.}

\item{force}{Ignore every warning and just use this data.}

\item{...}{Extra arguments passed to arglist.}
}
\value{
dataset suitable for limma analysis
}
\description{
Take an expt and poke at it to ensure that it will not result in troubled
results.
}
\details{
Invoked by deseq_pairwise() and edger_pairwise().
}
\seealso{
[DESeq2] [edgeR] [choose_basic_dataset()] [choose_limma_dataset()]
}
