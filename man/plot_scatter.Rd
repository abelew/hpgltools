% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_point.R
\name{plot_scatter}
\alias{plot_scatter}
\title{Make a pretty scatter plot between two sets of numbers.}
\usage{
plot_scatter(
  df,
  color = "black",
  xlab = NULL,
  xcol = NULL,
  ycol = NULL,
  ylab = NULL,
  alpha = 0.6,
  size = 2
)
}
\arguments{
\item{df}{Dataframe likely containing two columns.}

\item{color}{Color of the dots on the graph.}

\item{xlab}{Alternate x-axis label.}

\item{ylab}{Alternate x-axis label.}

\item{alpha}{Define how see-through the dots are.}

\item{size}{Size of the dots on the graph.}
}
\value{
Ggplot2 scatter plot.
}
\description{
This function tries to supplement a normal scatterplot with some information
describing the relationship between the columns of data plotted.
}
\examples{
\dontrun{
 plot_scatter(lotsofnumbers_intwo_columns)
}
}
\seealso{
[plot_linear_scatter()] [all_pairwise()]
}
