% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Spat_DownBoundary.R
\name{DownBoundary}
\alias{DownBoundary}
\title{Determine the boundaries of the requested GBIF data}
\usage{
DownBoundary(Left = NULL, Right = NULL, Bottom = NULL, Top = NULL)
}
\arguments{
\item{Left, Right, Bottom, Top}{Numeric, the left, right, bottom, and top
boundary of the area.}
}
\value{
A character string representing the WKT of the polygon that outlines
the specified boundaries.
}
\description{
This function constructs a Well-Known Text (WKT) string representing a
polygon that outlines the specified boundaries. It is used to define the area
of interest for downloading GBIF data through the \code{\link[rgbif:download_predicate_dsl]{rgbif::pred_within()}}
function.
}
\examples{
IASDT.R::DownBoundary(Left = 20, Right = 30, Bottom = 40, Top = 50)
}
\author{
Ahmed El-Gabbas
}
