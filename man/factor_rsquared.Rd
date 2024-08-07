% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimension_reduction.R
\name{factor_rsquared}
\alias{factor_rsquared}
\title{Collect the r^2 values from a linear model fitting between a singular
value decomposition and factor.}
\usage{
factor_rsquared(datum, fact, type = "factor")
}
\arguments{
\item{datum}{Result from corpcor::fast.svd.}

\item{fact}{Experimental factor from the original data.}

\item{type}{Make this categorical or continuous with factor/continuous.}
}
\value{
The r^2 values of the linear model as a percentage.
}
\description{
Collect the r^2 values from a linear model fitting between a singular
value decomposition and factor.
}
\seealso{
[corpcor] [stats::lm()]
}
