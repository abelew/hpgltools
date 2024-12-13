% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/de_xlsx.R
\name{combine_mapped_table}
\alias{combine_mapped_table}
\title{Combine data taken from map_keepers() into a single large table.}
\usage{
combine_mapped_table(
  entry,
  includes,
  adjp = TRUE,
  padj_type = "fdr",
  annot_df = NULL,
  excludes = NULL,
  lfc_cutoff = 1,
  p_cutoff = 0.05,
  format_sig = 4,
  sheet_count = 0,
  keep_underscore = TRUE,
  wanted_genes = NULL,
  scale_p = FALSE
)
}
\arguments{
\item{entry}{Single entry from map_keepers() which provides
orientation information about the table from all_pairwise(), along
with the actual data.}

\item{includes}{List of methods to include.}

\item{adjp}{Used adjusted pvalues when defining 'significant.?}

\item{padj_type}{Perform this type of pvalue adjustment.}

\item{annot_df}{Include these annotations in the result tables.}

\item{excludes}{When provided as a list, remove any rows with values in the column defined by the list names, otherwise exclude rownames.}

\item{lfc_cutoff}{Use this value for a log2FC significance cutoff.}

\item{p_cutoff}{Use this value for a(n adjusted) pvalue
significance cutoff.}

\item{format_sig}{Use this many significant digits for some of the
unwieldy numbers.}

\item{sheet_count}{Start with these sheet number and increment for excel.}
}
\description{
This is part of an ongoing attempt to simplify and clean up the
combine_de_tables() function.  I am hoping that map_keepers and
this will be able to take over all the logic currently held in the
various extract_keepers_xxx() functions.
}
