% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variants.R
\name{snps_vs_genes}
\alias{snps_vs_genes}
\title{Make a summary of the observed snps by gene ID.}
\usage{
snps_vs_genes(
  expt,
  snp_result,
  start_col = "start",
  end_col = "end",
  snp_name_col = "seqnames",
  observed_in = NULL,
  expt_name_col = "chromosome",
  ignore_strand = TRUE
)
}
\arguments{
\item{expt}{The original expressionset.}

\item{snp_result}{The result from get_snp_sets().}

\item{start_col}{Which column provides the start of each gene?}

\item{end_col}{and the end column of each gene?}

\item{snp_name_col}{Name of the column in the metadata with the sequence names.}

\item{observed_in}{Minimum proportion of samples required before this is deemed real.}

\item{expt_name_col}{Name of the metadata column with the chromosome names.}

\item{ignore_strand}{Ignore strand information when returning?}
}
\value{
List with some information by gene.
}
\description{
Instead of cross referencing variant positions against experimental
condition, one might be interested in seeing what variants are observed per
gene.  This function attempts to answer that question.
}
\examples{
 \dontrun{
 expt <- create_expt(metadata, gene_information)
 snp_expt <- count_expt_snps(expt)
 snp_result <- get_snp_sets(snp_expt)
 gene_intersections <- snps_vs_genes(expt, snp_result)
}
}
\seealso{
[GenomicRanges::makeGRangesFromDataFrame()] [IRanges::subsetByOverlaps()]
 [IRanges::mergeByOverlaps()] [IRanges::countOverlaps()]
}
