% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_xlsx.R
\name{extract_abundant_genes}
\alias{extract_abundant_genes}
\title{Extract the sets of genes which are significantly more abundant than the rest.}
\usage{
extract_abundant_genes(
  pairwise,
  according_to = "deseq",
  n = 100,
  z = NULL,
  unique = FALSE,
  excel = "excel/abundant_genes.xlsx",
  ...
)
}
\arguments{
\item{pairwise}{Output from _pairwise()().}

\item{according_to}{What tool(s) define 'most?'  One may use deseq, edger,
limma, basic, all.}

\item{n}{How many genes to pull?}

\item{z}{Instead take the distribution of abundances and pull those past the
given z score.}

\item{unique}{One might want the subset of unique genes in the top-n which
are unique in the set of available conditions.  This will attempt to
provide that.}

\item{excel}{Excel file to write.}

\item{...}{Arguments passed into arglist.}
}
\value{
The set of most/least abundant genes by contrast/tool.
}
\description{
Given the output of something_pairwise(), pull out the genes for each contrast
which are the most/least abundant.  This is in contrast to extract_significant_genes().
That function seeks out the most changed, statistically significant genes.
}
\seealso{
\pkg{openxlsx}
}
