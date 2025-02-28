## |------------------------------------------------------------------------| #
# GridDiagOff ------
## |------------------------------------------------------------------------| #

#' Create a `multilinestring` sf object for the diagonal and off-diagonal lines
#' for each grid cell
#'
#' This function takes an `sf` object representing a grid and creates a new sf
#' object where each grid cell is represented by a `multilinestring` geometry
#' consisting of its diagonal and off-diagonal lines.
#' @name GridDiagOff
#' @param DT An sf object (tibble) representing a spatial grid. The function
#'   expects this object to have a geometry column with polygons representing
#'   grid cells. If `NULL`, the function will stop with an error message.
#' @return An `sf` object where each row corresponds to a grid cell from the
#'   input, represented by a `multilinestring` geometry of its diagonal and
#'   off-diagonal lines. The returned object retains the coordinate reference
#'   system (CRS) of the input.
#' @author Ahmed El-Gabbas
#' @export
#' @note The function requires the `sf`, `dplyr`, `purrr`, `tibble`, and `tidyr`
#'   packages to be installed and loaded.
#' @seealso GridCross
#' @examples
#' IASDT.R::LoadPackages(List = c("dplyr", "sf", "raster", "ggplot2"))
#'
#' Grid <- raster::raster(nrows = 10, ncols = 10, xmn = 0, xmx = 10,
#'                        ymn = 0, ymx = 10, crs = 4326) %>%
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
#' Grid_X <- GridDiagOff(DT = Grid)
#'
#' ggplot2::ggplot() +
#'   ggplot2::geom_sf(Grid, mapping = ggplot2::aes(), color = "black",
#'                    linewidth = 0.5, fill = "transparent") +
#'   ggplot2::geom_sf(Grid_X, mapping = ggplot2::aes(), color = "red",
#'                    linewidth = 0.5, inherit.aes = TRUE) +
#'   ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::scale_y_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::theme_minimal()

GridDiagOff <- function(DT = NULL) {

  if (is.null(DT)) {
    stop("Input DT  cannot be NULL", call. = FALSE)
  }

  DT %>%
    dplyr::pull("geometry") %>%
    purrr::map(
      .f = ~{
        SS2 <- .x %>%
          sf::st_coordinates() %>%
          tibble::as_tibble() %>%
          dplyr::select(X, Y)

        OffDiag <- dplyr::bind_rows(
          dplyr::filter(SS2, X == min(X), Y == max(Y)),
          dplyr::filter(SS2, X == max(X), Y == min(Y))) %>%
          as.matrix() %>%
          sf::st_linestring()

        Diag <- dplyr::bind_rows(
          dplyr::filter(SS2, X == min(X), Y == min(Y)),
          dplyr::filter(SS2, X == max(X), Y == max(Y))) %>%
          as.matrix() %>%
          sf::st_linestring()

        Out <- list(OffDiag, Diag) %>%
          sf::st_multilinestring() %>%
          sf::st_geometry()  %>%
          sf::st_set_crs(sf::st_crs(DT)) %>%
          sf::st_sfc()

        Out

      }) %>%
    tibble::tibble(geometry = .) %>%
    tidyr::unnest("geometry") %>%
    sf::st_as_sf()
}
