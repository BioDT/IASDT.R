# |---------------------------------------------------| #
# PrepKnots ----
# |---------------------------------------------------| #

#' Prepare knot locations for GPP models.
#'
#' Prepare knot locations for GPP models.
#' 
#' @param Coords A matrix or data frame for the coordinates of the sampling units.
#' @param MinDist Numeric. Minimum distance in meters for the `knotDist` and `minKnotDist` arguments of the `Hmsc::constructKnots` function.
#' @param JitterDist Numeric. Distance in meter to be used in the `sf::st_jitter` function.
#' @name PrepKnots
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

PrepKnots <- function(Coords = DT_xy, MinDist, JitterDist = 100) {
  # coordinates of the sampling units as sf object
  SU_Sf <- Coords %>%
    tibble::as_tibble() %>%
    sf::st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 3035)

  # Prepare knots
  Knots <- Hmsc::constructKnots(
    Coords, knotDist = MinDist, minKnotDist = MinDist) %>%
    tibble::as_tibble() %>%
    # Hmsc::constructKnots may return duplicated points; discard them
    dplyr::distinct() %>%
    dplyr::mutate(
      Identical = purrr::map2_lgl(
        .x = Var1, .y = Var2,
        .f = ~{
          SU_Sf %>%
            dplyr::filter(x == .x, y == .y) %>%
            nrow() %>%
            magrittr::is_greater_than(0)
        }, .progress = TRUE))

  if (any(Knots$Identical)) {
    Knots <- Knots %>%
      sf::st_as_sf(coords = c("Var1", "Var2"), remove = FALSE, crs = 3035) %>%
      dplyr::mutate(
        geometry = sf::st_sfc(
          dplyr::if_else(
            Identical, sf::st_jitter(geometry, JitterDist), geometry))) %>%
      sf::st_coordinates() %>%
      as.data.frame() %>%
      stats::setNames(c("Var1", "Var2"))
  } else {
    Knots <- Knots %>%
      dplyr::select(Var1, Var2) %>%
      as.data.frame()
  }
  Hmsc::HmscRandomLevel(sData = Coords, sMethod = "GPP", sKnot = Knots)
}
