% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation_biomart.R
\name{load_biomart_orthologs}
\alias{load_biomart_orthologs}
\title{Use biomart to get orthologs between supported species.}
\usage{
load_biomart_orthologs(
  gene_ids = NULL,
  first_species = "hsapiens",
  second_species = "mmusculus",
  host = NULL,
  trymart = "ENSEMBL_MART_ENSEMBL",
  archive = TRUE,
  default_hosts = c("useast.ensembl.org", "uswest.ensembl.org", "www.ensembl.org",
    "asia.ensembl.org"),
  year = NULL,
  month = NULL,
  trydataset = NULL,
  attributes = "ensembl_gene_id"
)
}
\arguments{
\item{gene_ids}{List of gene IDs to translate.}

\item{first_species}{Linnean species name for one species.}

\item{second_species}{Linnean species name for the second species.}

\item{host}{Ensembl server to query.}

\item{trymart}{Assumed mart name to use.}

\item{archive}{Use an archive server?}

\item{default_hosts}{Set of default hosts to query.}

\item{year}{When using an archive server, use this year (otherwise it will choose last year).}

\item{month}{When using an archive server, use this month (otherwise, this month).}

\item{trydataset}{Choose a dataset to query.}

\item{attributes}{Key to query}
}
\value{
list of 4 elements:  The first is the set of all ids, as getLDS seems
 to always send them all; the second is the subset corresponding to the
 actual ids of interest, and the 3rd/4th are other, optional ids from other datasets.
}
\description{
Biomart's function getLDS is incredibly powerful, but it makes me think very
polite people are going to start knocking on my door, and it fails weirdly
pretty much always. This function attempts to alleviate some of that frustration.
}
\details{
Tested in test_40ann_biomart.R
As with my other biomart functions, this one grew out of frustrations when
attempting to work with the incredibly unforgiving biomart service.  It does
not attempt to guarantee a useful biomart connection, but will hopefully
point out potentially correct marts and attributes to use for a successful
query.  I can say with confidence that it works well between mice and
humans.
}
\examples{
 mouse_yeast_orthologs <- load_biomart_orthologs(gene_ids = NULL, first_species = "mmusculus",
                                                 second_species = "scerevisiae")
 head(mouse_yeast_orthologs$all_linked_genes)
}
\seealso{
[biomaRt::getLDS()]
}
