% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_bar.R
\name{plot_sample_bars}
\alias{plot_sample_bars}
\title{The actual library size plotter.}
\usage{
plot_sample_bars(
  sample_df,
  condition = NULL,
  colors = NULL,
  integerp = FALSE,
  order = NULL,
  text = TRUE,
  plot_title = NULL,
  yscale = NULL,
  ...
)
}
\arguments{
\item{sample_df}{Expt, dataframe, or expressionset of samples.}

\item{condition}{Vector of sample condition names.}

\item{colors}{Color scheme if the data is not an expt.}

\item{integerp}{Is this comprised of integer values?}

\item{order}{Explicitly set the order of samples in the plot?}

\item{text}{Add the numeric values inside the top of the bars of the plot?}

\item{plot_title}{Title for the plot.}

\item{yscale}{Whether or not to log10 the y-axis.}

\item{...}{Used to catch random arguments which are unused here.}
}
\description{
This makes a ggplot2 plot of library sizes.
}
