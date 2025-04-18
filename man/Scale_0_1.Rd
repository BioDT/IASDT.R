% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spat_scale_0_1.R
\name{scale_0_1}
\alias{scale_0_1}
\title{Scale a \code{SpatRaster} object values to a range between 0 and 1}
\usage{
scale_0_1(raster)
}
\arguments{
\item{raster}{SpatRaster; The SpatRaster object to be scaled.}
}
\value{
SpatRaster; A SpatRaster object with all values scaled between 0 and
1.
}
\description{
Scale a \code{SpatRaster} object (from the \code{terra} R package in) values to a range
between 0 and 1.
}
\note{
This function takes a SpatRaster object as input, calculates its
minimum and maximum values, and scales all its values to a range between 0
and 1. This is useful for normalization purposes in spatial analysis and
modeling. The function relies on the \code{terra} package for spatial data
manipulation.
}
\examples{
library(terra)
r <- rast(ncols = 10, nrows = 10)
values(r) <- seq_len(terra::ncell(r))

scale_0_1(r)
}
\author{
Ahmed El-Gabbas
}
