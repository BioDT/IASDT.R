## |------------------------------------------------------------------------| #
# GridDiagOff ------
## |------------------------------------------------------------------------| #

#' Create a `multilinestring` (diagonal and off-diagonal lines) sf object from each grid cell
#'
#' Create a `multilinestring` (diagonal and off-diagonal lines) sf object from each grid cell
#'
#' @name GridDiagOff
#' @param DT input sf tibble
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' IASDT.R::LoadPackages(dplyr, sf, raster, ggplot2)
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
#' Grid_X <- GridDiagOff(Grid)
#'
#' ggplot2::ggplot() +
#'   ggplot2::geom_sf(Grid, mapping = ggplot2::aes(), color = "black",
#'                    linewidth = 0.5, fill = "transparent") +
#'   ggplot2::geom_sf(Grid_X, mapping = ggplot2::aes(), color = "red",
#'                    linewidth = 0.5, inherit.aes = TRUE) +
#'   ggplot2::scale_x_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::scale_y_continuous(expand = c(0, 0, 0, 0), limits = c(0, 10)) +
#'   ggplot2::theme_minimal()

GridDiagOff <- function(DT) {

  InputCRS <- sf::st_crs(DT)

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

        list(OffDiag, Diag) %>%
          sf::st_multilinestring() %>%
          sf::st_geometry()  %>%
          sf::st_set_crs(InputCRS) %>%
          sf::st_sfc() %>%
          return()

      }) %>%
    tibble::tibble(geometry = .) %>%
    tidyr::unnest("geometry") %>%
    sf::st_as_sf()
}
