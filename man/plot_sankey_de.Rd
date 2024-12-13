% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_plots.R
\name{plot_sankey_de}
\alias{plot_sankey_de}
\title{Make a sankey plot showing how the number of genes deemed significant is constrained.}
\usage{
plot_sankey_de(
  de_table,
  lfc = 1,
  p = 0.05,
  lfc_column = "deseq_logfc",
  p_column = "deseq_adjp"
)
}
\arguments{
\item{de_table}{The result from combine_de_tables()}

\item{lfc}{FC constraint.}

\item{p}{P-value constraint.}

\item{lfc_column}{Dataframe column from which to acquire the FC
values.}

\item{p_column}{Dataframe column from which to acquire the
p-values.}
}
\value{
A fun sankey plot!
}
\description{
Ideally, this should show how adding various Fc/p-value constraints
on the definition of 'significant' decreases the number of genes
one is likely to look at.
}
