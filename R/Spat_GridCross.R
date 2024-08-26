## |------------------------------------------------------------------------| #
# GridCross ------
## |------------------------------------------------------------------------| #

#' Create a `multilinestring` sf object representing cross in the middle of each
#' grid cell
#'
#' Create a `multilinestring` (cross in the middle of the grid) sf object from
#' each grid cell
#' @name GridCross
#' @param DT An `sf` object (tibble) representing grid cells. The function
#'   expects this object to have a geometry column with polygon geometries. If
#'   `NULL`, the function will stop with an error message.
#' @return An `sf` object with `multilinestring` geometries representing crosses
#'   in the middle of each input grid cell. The returned object has the same CRS
#'   (Coordinate Reference System) as the input `DT`
#' @author Ahmed El-Gabbas
#' @export
#' @seealso \code{\link{GridDiagOff}} for creating diagonal lines in grid cells.
#' @note The function requires the `sf`, `dplyr`, `purrr`, `tibble`, and `tidyr`
#'   packages to be installed and loaded.
#' @examples
#' IASDT.R::LoadPackages(List = c("dplyr", "sf", "raster", "ggplot2"))
#'
#' Grid <- raster::raster(
#'     nrows = 10, ncols = 10, xmn = 0, xmx = 10,
#'     ymn = 0, ymx = 10, crs = 4326) %>%
#'   setNames("Grid") %>%
#'   raster::setValues(1) %>%
#'   raster::rasterToPolygons() %>%
#'   sf::st_as_sf()
#'
#' ggplot2::ggplot() +
#'   ggplot2::geom_sf(Grid, mapping = ggplot2::aes(), color = "black",
#'                    linewidth = 0.5, fill = "transparent") +
#'   ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::scale_y_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::theme_minimal()
#'
#' Grid_X <- GridCross(Grid)
#'
#' ggplot2::ggplot() +
#'   ggplot2::geom_sf(Grid, mapping = ggplot2::aes(), color = "black",
#'                    linewidth = 0.5, fill = "transparent") +
#'   ggplot2::geom_sf(Grid_X, mapping = ggplot2::aes(), color = "red",
#'                    linewidth = 0.5, inherit.aes = TRUE) +
#'   ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::scale_y_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::theme_minimal()

GridCross <- function(DT = NULL) {

  if (is.null(DT)) {
    stop("Input DT  cannot be NULL", call. = FALSE)
  }

  InputCRS <- sf::st_crs(DT)

  DT %>%
    dplyr::pull("geometry") %>%
    purrr::map(
      .f = ~{
        SS2 <- .x %>%
          sf::st_coordinates() %>%
          tibble::as_tibble() %>%
          dplyr::select(X, Y)

        Horiz <- SS2 %>%
          dplyr::reframe(
            X = X,
            Y = (min(Y) + (max(Y) - min(Y)) / 2)) %>%
          dplyr::distinct() %>%
          as.matrix() %>%
          sf::st_linestring()
        Vert <- SS2 %>%
          dplyr::reframe(
            X = (min(X) + (max(X) - min(X)) / 2),
            Y = Y) %>%
          dplyr::distinct() %>%
          as.matrix() %>%
          sf::st_linestring()

        list(Horiz, Vert) %>%
          sf::st_multilinestring() %>%
          sf::st_geometry()  %>%
          sf::st_set_crs(InputCRS) %>%
          sf::st_sfc() %>%
          return()
      }) %>%
    tibble::tibble(geometry = .) %>%
    tidyr::unnest("geometry") %>%
    sf::st_as_sf() %>%
    return()
}
