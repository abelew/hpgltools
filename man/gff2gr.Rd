% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotation_gff.R
\name{gff2gr}
\alias{gff2gr}
\title{Rewrite a gff file as a granges with full seqinfo if possible.}
\usage{
gff2gr(gff, type = NULL, type_column = "type")
}
\arguments{
\item{gff}{Input gff file.}

\item{type}{Feature type to extract.}

\item{type_column}{Tag from the gff file to use when extracting the type.}
}
\description{
Rewrite a gff file as a granges with full seqinfo if possible.
}
