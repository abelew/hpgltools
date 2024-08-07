% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_misc.R
\name{pp}
\alias{pp}
\title{Plot a picture, with hopefully useful options for most(any) format.}
\usage{
pp(file, image = NULL, width = 9, height = 9, res = 180, ...)
}
\arguments{
\item{file}{Filename to write}

\item{image}{Optionally, add the image you wish to plot and this will both
print it to file and screen.}

\item{width}{How wide?}

\item{height}{How high?}

\item{res}{The chosen resolution.}

\item{...}{Arguments passed to the image plotters.}
}
\value{
a png/svg/eps/ps/pdf with height = width=9 inches and a high resolution
}
\description{
This calls svg/png/postscript/etc according to the filename provided.
}
\seealso{
[png()] [svg()] [postscript()] [cairo_ps()] [cairo_pdf()] [tiff()] [devEMF::emf()]
 [jpg()] [bmp()]
}
