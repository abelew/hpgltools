% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_hist.R
\name{plot_multihistogram}
\alias{plot_multihistogram}
\title{Make a pretty histogram of multiple datasets.}
\usage{
plot_multihistogram(
  data,
  log = FALSE,
  binwidth = NULL,
  bins = NULL,
  colors = NULL
)
}
\arguments{
\item{data}{Dataframe of lots of pretty numbers, this also accepts lists.}

\item{log}{Plot the data on the log scale?}

\item{binwidth}{Set a static bin width with an unknown # of bins?  If neither of these are
provided, then bins is set to 500, if both are provided, then bins wins.}

\item{bins}{Set a static # of bins of an unknown width?}

\item{colors}{Change the default colors of the densities?}
}
\value{
List of the ggplot histogram and some statistics describing the distributions.
}
\description{
If there are multiple data sets, it might be useful to plot them on a
histogram together and look at the t.test results between distributions.
}
\examples{
\dontrun{
 kittytime = plot_multihistogram(df)
}
}
\seealso{
[stats::pairwise.t.test()] [ggplot2]
}
