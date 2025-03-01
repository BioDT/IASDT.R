% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Spat_setRastCRS.R
\name{setRastCRS}
\alias{setRastCRS}
\title{sets CRS for a SpatRaster}
\usage{
setRastCRS(R, CRS = "epsg:3035")
}
\arguments{
\item{R}{A SpatRaster object whose CRS needs to be set.}

\item{CRS}{Character. CRS value to be set, default is "epsg:3035".}
}
\value{
The SpatRaster object with the updated CRS.
}
\description{
This function sets the coordinate reference system (CRS) for a SpatRaster
object using the specified EPSG code. This is a wrapper function for
\code{terra::crs(R) <- CRS} but allowing to set the CRS in the pipe.
}
\author{
Ahmed El-Gabbas
}
