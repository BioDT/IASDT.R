% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Spat_GridCross.R
\name{GridCross}
\alias{GridCross}
\title{Create a \code{multilinestring} sf object representing cross in the middle of each
grid cell}
\usage{
GridCross(DT = NULL)
}
\arguments{
\item{DT}{An \code{sf} object (tibble) representing grid cells. The function
expects this object to have a geometry column with polygon geometries. If
\code{NULL}, the function will stop with an error message.}
}
\value{
An \code{sf} object with \code{multilinestring} geometries representing crosses
in the middle of each input grid cell. The returned object has the same CRS
(Coordinate Reference System) as the input \code{DT}
}
\description{
Create a \code{multilinestring} (cross in the middle of the grid) sf object from
each grid cell
}
\note{
The function requires the \code{sf}, \code{dplyr}, \code{purrr}, \code{tibble}, and \code{tidyr}
packages to be installed and loaded.
}
\examples{
IASDT.R::LoadPackages(List = c("dplyr", "sf", "raster", "ggplot2"))

Grid <- raster::raster(
    nrows = 10, ncols = 10, xmn = 0, xmx = 10,
    ymn = 0, ymx = 10, crs = 4326) \%>\%
  setNames("Grid") \%>\%
  raster::setValues(1) \%>\%
  raster::rasterToPolygons() \%>\%
  sf::st_as_sf()

ggplot2::ggplot() +
  ggplot2::geom_sf(Grid, mapping = ggplot2::aes(), color = "black",
                   linewidth = 0.5, fill = "transparent") +
  ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
  ggplot2::scale_y_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
  ggplot2::theme_minimal()

Grid_X <- GridCross(Grid)

ggplot2::ggplot() +
  ggplot2::geom_sf(Grid, mapping = ggplot2::aes(), color = "black",
                   linewidth = 0.5, fill = "transparent") +
  ggplot2::geom_sf(Grid_X, mapping = ggplot2::aes(), color = "red",
                   linewidth = 0.5, inherit.aes = TRUE) +
  ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
  ggplot2::scale_y_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
  ggplot2::theme_minimal()
}
\seealso{
\code{\link{GridDiagOff}} for creating diagonal lines in grid cells.
}
\author{
Ahmed El-Gabbas
}
