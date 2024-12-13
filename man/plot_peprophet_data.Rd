% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_proteomics.R
\name{plot_peprophet_data}
\alias{plot_peprophet_data}
\title{Plot some data from the result of extract_peprophet_data()}
\usage{
plot_peprophet_data(
  table,
  xaxis = "precursor_neutral_mass",
  xscale = NULL,
  yaxis = "num_matched_ions",
  yscale = NULL,
  size_column = "prophet_probability",
  ...
)
}
\arguments{
\item{table}{Big honking data table from extract_peprophet_data()}

\item{xaxis}{Column to plot on the x-axis}

\item{xscale}{Change the scale of the x-axis?}

\item{yaxis}{guess!}

\item{yscale}{Change the scale of the y-axis?}

\item{size_column}{Use a column for scaling the sizes of dots in the plot?}

\item{...}{extra options which may be used for plotting.}
}
\value{
a plot!
}
\description{
extract_peprophet_data() provides a ridiculously large data table of a comet
result after processing by RefreshParser and xinteract/peptideProphet.
This table has some 37-ish columns and I am not entirely certain which ones
are useful as diagnostics of the data.  I chose a few and made options to
pull some/most of the rest.  Lets play!
}
