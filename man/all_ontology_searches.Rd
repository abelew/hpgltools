% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ontology_shared.R
\name{all_ontology_searches}
\alias{all_ontology_searches}
\title{Perform ontology searches given the results of a differential expression analysis.}
\usage{
all_ontology_searches(
  de_out,
  gene_lengths = NULL,
  goids = NULL,
  n = NULL,
  z = NULL,
  lfc = NULL,
  p = NULL,
  overwrite = FALSE,
  species = "unsupported",
  orgdb = "org.Dm.eg.db",
  goid_map = "reference/go/id2go.map",
  gff_file = NULL,
  gff_type = "gene",
  do_goseq = TRUE,
  do_cluster = TRUE,
  do_topgo = TRUE,
  do_gostats = TRUE,
  do_gprofiler = TRUE,
  do_trees = FALSE,
  ...
)
}
\arguments{
\item{de_out}{List of topTables comprising limma/deseq/edger outputs.}

\item{gene_lengths}{Data frame of gene lengths for goseq.}

\item{goids}{Data frame of goids and genes.}

\item{n}{Number of genes at the top/bottom of the fold-changes to define 'significant.'}

\item{z}{Number of standard deviations from the mean fold-change used to
define 'significant.'}

\item{lfc}{Log fold-change used to define 'significant'.}

\item{p}{Maximum pvalue to define 'significant.'}

\item{overwrite}{Overwrite existing excel results file?}

\item{species}{Supported organism used by the tools.}

\item{orgdb}{Provide an organismDbi/Orgdb to hold the various annotation
data, in response to the shift of clusterprofiler and friends towards using them.}

\item{goid_map}{Mapping file used by topGO, if it does not exist then
goids_df creates it.}

\item{gff_file}{gff file containing the annotations used by gff2genetable
from clusterprofiler.}

\item{gff_type}{Column to use from the gff file for the universe of genes.}

\item{do_goseq}{Perform simple_goseq()?}

\item{do_cluster}{Perform simple_clusterprofiler()?}

\item{do_topgo}{Perform simple_topgo()?}

\item{do_gostats}{Perform simple_gostats()?}

\item{do_gprofiler}{Perform simple_gprofiler()?}

\item{do_trees}{make topGO trees from the data?}

\item{...}{Arguments to pass through in arglist.}
}
\value{
a list of up/down ontology results from
 goseq/clusterprofiler/topgo/gostats, and associated trees.
}
\description{
This takes a set of differential expression results, extracts a subset of
up/down expressed genes; passes them to goseq, clusterProfiler, topGO,
GOstats, and gProfiler; collects the outputs; and returns them in a
(hopefully) consistent fashion.  It attempts to handle the differing required
annotation/GOid inputs required for each tool and/or provide supported
species in ways which the various tools expect.
}
\examples{
\dontrun{
 many_comparisons = limma_pairwise(expt = an_expt)
 tables = many_comparisons$limma
 this_takes_forever = limma_ontology(tables, gene_lengths = lengthdb,
                                     goids = goids_df, z = 1.5, gff_file='length_db.gff')
}
}
\seealso{
[goseq] [clusterProfiler] [topGO] [goStats] [gProfiler] [GO.db]
}
