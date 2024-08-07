% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gene_sets.R
\name{make_gsc_from_abundant}
\alias{make_gsc_from_abundant}
\title{Given a pairwise result, make a gene set collection.}
\usage{
make_gsc_from_abundant(
  pairwise,
  according_to = "deseq",
  annotation_name = "org.Hs.eg.db",
  researcher_name = "elsayed",
  study_name = "macrophage",
  category_name = "infection",
  phenotype_name = NULL,
  pair_names = "high",
  current_id = "ENSEMBL",
  required_id = "ENTREZID",
  ...
)
}
\arguments{
\item{pairwise}{A pairwise result, or combined de result, or extracted genes.}

\item{according_to}{When getting significant genes, use this method.}

\item{annotation_name}{Name of the annotation dataset.}

\item{researcher_name}{Prefix of the name for the generated set(s).}

\item{study_name}{Second element in the name of the generated set(s).}

\item{category_name}{Third element in the name of the generated set(s).}

\item{phenotype_name}{Optional phenotype data for the generated set(s).}

\item{pair_names}{The suffix of the generated set(s).}

\item{current_id}{What type of ID is the data currently using?}

\item{required_id}{What type of ID should the use?}

\item{...}{Extra arguments for extract_abundant_genes().}
}
\value{
List containing 3 GSCs, one containing both the highs/lows called
 'colored', one of the highs, and one of the lows.
}
\description{
If I want to play with gsva and friends, then I need GeneSetCollections!
Much like make_gsc_from_significant(), this function extract the genes deemed
'abundant' and generates gene sets accordingly.
}
\seealso{
[extract_abundant_genes()] [make_gsc_from_ids()] [GSEABase]
}
