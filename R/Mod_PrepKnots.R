## |------------------------------------------------------------------------| #
# PrepKnots ----
## |------------------------------------------------------------------------| #

#' Prepare knot locations for Hmsc GPP models
#'
#' This function prepares the locations of knots for use in GPP models within
#' the HMSC framework. It ensures that knots are spaced at a minimum specified
#' distance and applies jitter to any identical coordinates to avoid overlap.
#' @param Coords A matrix or data frame containing the coordinates (x, y) of the
#'   sampling units.
#' @param MinDist A numeric value specifying the minimum distance (in the same
#'   units as Coords; here in meters) between knots. This distance is used for
#'   both `knotDist` and `minKnotDist` parameters of the [Hmsc::constructKnots]
#'   function.
#' @param JitterDist A numeric value specifying the distance (in the same units
#'   as Coords; here in meters) for jittering identical coordinates. Defaults to
#'   100.
#' @name PrepKnots
#' @author Ahmed El-Gabbas
#' @return An object suitable for specifying the random level in HMSC GPP
#'   models. This object contains the prepared knot locations.
#' @export

PrepKnots <- function(Coords = NULL, MinDist = NULL, JitterDist = 100) {

  if (is.null(Coords) || is.null(MinDist)) {
    stop("Both 'Coords' and 'MinDist' must be provided and cannot be NULL")
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  DT_xy <- Var1 <- Var2 <- Identical <- geometry <- NULL

  # coordinates of the sampling units as sf object
  SU_Sf <- tibble::as_tibble(Coords) %>%
    stats::setNames(c("x", "y")) %>%
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
        }))

  if (any(Knots$Identical)) {
    Knots <- sf::st_as_sf(
      x = Knots, coords = c("Var1", "Var2"), remove = FALSE, crs = 3035) %>%
      dplyr::mutate(
        geometry = sf::st_sfc(
          dplyr::if_else(
            Identical, sf::st_jitter(geometry, JitterDist), geometry))) %>%
      sf::st_coordinates() %>%
      as.data.frame() %>%
      stats::setNames(c("Var1", "Var2"))
  } else {
    Knots <- as.data.frame(dplyr::select(Knots, Var1, Var2))
  }
  return(Hmsc::HmscRandomLevel(sData = Coords, sMethod = "GPP", sKnot = Knots))
}
