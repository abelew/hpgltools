% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alt_splicing.R
\name{plot_rmats}
\alias{plot_rmats}
\title{Given some psi and tpm data from rMATS, make a pretty plot!}
\usage{
plot_rmats(
  se = NULL,
  a5ss = NULL,
  a3ss = NULL,
  mxe = NULL,
  ri = NULL,
  sig_threshold = 0.05,
  dpsi_threshold = 0.7,
  label_type = NULL,
  alpha = 0.7
)
}
\arguments{
\item{se}{Table of skipped exon data from rmats.}

\item{a5ss}{Table of alternate 5p exons.}

\item{a3ss}{Table of alternate 3p exons.}

\item{mxe}{Table of alternate exons.}

\item{ri}{Table of retained introns.}

\item{sig_threshold}{Use this significance threshold.}

\item{dpsi_threshold}{Use a delta threshold.}

\item{label_type}{Choose a type of event to label.}

\item{alpha}{How see-through should the points be in the plot?}
}
\value{
List containing the plot and some of the requisite data.
}
\description{
This should take either a dataframe or filename for the psi data from rMATS.
This was mostly copy/pasted from plot_suppa().
}
\examples{
 \dontrun{
 rmats_plot <- plot_rmats(se_table, a5_table, a3_table)
}
}
\seealso{
[plot_supps()]
}
