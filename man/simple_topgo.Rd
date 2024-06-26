% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ontology_topgo.R
\name{simple_topgo}
\alias{simple_topgo}
\title{Perform a simplified topgo analysis.}
\usage{
simple_topgo(
  sig_genes,
  goid_map = "id2go.map",
  go_db = NULL,
  pvals = NULL,
  limitby = "fisher",
  limit = 0.1,
  signodes = 100,
  sigforall = TRUE,
  numchar = 300,
  selector = "topDiffGenes",
  pval_column = "deseq_adjp",
  overwrite = FALSE,
  densities = FALSE,
  pval_plots = TRUE,
  parallel = TRUE,
  excel = NULL,
  ...
)
}
\arguments{
\item{sig_genes}{Data frame of differentially expressed genes, containing IDs
any other columns.}

\item{goid_map}{File containing mappings of genes to goids in the format
expected by topgo.}

\item{go_db}{Data frame of the goids which may be used to make the goid_map.}

\item{pvals}{Set of pvalues in the DE data which may be used to improve the
topgo results.}

\item{limitby}{Test to index the results by.}

\item{limit}{Ontology pvalue to use as the lower limit.}

\item{signodes}{I don't remember right now.}

\item{sigforall}{Provide the significance for all nodes?}

\item{numchar}{Character limit for the table of results.}

\item{selector}{Function name for choosing genes to include.}

\item{pval_column}{Column from which to acquire scores.}

\item{overwrite}{Yeah I do not remember this one either.}

\item{densities}{Densities, yeah, the densities...}

\item{pval_plots}{Include pvalue plots of the results a la clusterprofiler?}

\item{excel}{Print the results to an excel file?}

\item{...}{Other options which I do not remember right now!}
}
\value{
Big list including the various outputs from topgo
}
\description{
This will attempt to make it easier to run topgo on a set of genes.
The way I organized these data structures is completely stupid. I
want to convert the data from topgo to clusterprofiler for ease of
plotting, but because of the terrible way I organized everything
that is likely to be difficult.
}
\seealso{
[topGO]
}
