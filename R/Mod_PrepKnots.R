## |------------------------------------------------------------------------| #
# PrepKnots ----
## |------------------------------------------------------------------------| #

#' Prepare knot locations for Hmsc GPP models
#'
#' This function prepares the locations of knots for use in GPP models within
#' the HMSC framework. It ensures that knots are spaced at a minimum specified
#' distance and applies jitter to any identical coordinates to avoid overlap.
#' @param Coords A numeric matrix or data frame containing the (x, y)
#'   coordinates of sampling units.
#' @param MinDist A numeric value specifying the minimum distance between knots
#'   in meters. This distance is used for both `knotDist` and `minKnotDist`
#'   parameters of the [Hmsc::constructKnots] function.
#' @param JitterDist A numeric value for the jitter distance applied to
#'   overlapping coordinates to avoid exact duplicates. Defaults to 100 meters.
#' @param MinLF,MaxLF integer. Minimum and maximum number of latent factors to
#'   be used. Both default to `NULL` which means that the number of latent
#'   factors will be estimated from the data. If either is provided, the
#'   respective values will be used as arguments to [Hmsc::setPriors].
#' @param Alphapw Prior for the alpha parameter. Defaults to a list with `Prior
#'   = NULL`, `Min = 20`, `Max = 1200`, and `Samples = 200`. If all list items
#'   are NULL, the default prior will be used. If `Prior` is a matrix, it will
#'   be used as the prior. If `Prior` is `NULL`, the prior will be generated
#'   using `Min`, `Max`, and `Samples`. `Min` and `Max` are the minimum and
#'   maximum values of the alpha parameter (in kilometer). `Samples` is the
#'   number of samples to be used in the prior.
#'
#' @name PrepKnots
#' @author Ahmed El-Gabbas
#' @return An object suitable for specifying the random level in HMSC GPP
#'   models. This object contains the prepared knot locations.
#' @export

PrepKnots <- function(
    Coords = NULL, MinDist = NULL, JitterDist = 100,
    MinLF = NULL, MaxLF = NULL,
    Alphapw = list(Prior = NULL, Min = 20, Max = 1200, Samples = 200)) {
  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Var1 <- Var2 <- Identical <- geometry <- NULL

  # # ..................................................................... ###

  AlphaPrior <- function(MinVal, MaxVal, NSamples) {
    seq(MinVal, MaxVal, length.out = (NSamples - 1)) %>%
      magrittr::multiply_by(1000) %>%
      round() %>%
      c(0, .) %>%
      cbind(c(0.5, rep(0.5 / (NSamples - 1), (NSamples - 1))))
  }

  # # ..................................................................... ###

  # Validate Alphapw input ----

  if (all(purrr::map_lgl(Alphapw, is.null))) {
    Prior <- NULL
  } else {
    if (is.null(Alphapw$Prior)) {
      Prior <- AlphaPrior(
        MinVal = Alphapw$Min, MaxVal = Alphapw$Max, NSamples = Alphapw$Samples)
    } else {
      Prior <- Alphapw$Prior

      # Check if Prior is a matrix
      if (!inherits(Prior, "matrix")) {
        stop("`Prior` must be a matrix", call. = FALSE)
      }

      # Check if Prior has 2 columns
      if (ncol(Prior) != 2) {
        stop("`Prior` should have exactly 2 columns.", call. = FALSE)
      }

      # Check if Prior has at least 100 rows
      if (nrow(Prior) < 100) {
        stop("`Prior` must have >= 100 rows.", call. = FALSE)
      }

      # Check if the first element of the second column is 0.5
      # if (Prior[1, 2] != 0.5) {
      #   stop("The first element of the second column is not 0.5.",
      #        call. = FALSE)
      # }

      # Check if all values are positive
      if (any(as.vector(Prior) < 0)) {
        stop("`Prior` can not contain negative numbers", call. = FALSE)
      }

      # Check if the sum of the second column is equal to 1
      if (sum(Prior[, 2]) != 1) {
        stop("The sum of the second column is not equal to 1.", call. = FALSE)
      }
    }
  }

  # # ..................................................................... ###

  # Validation for latent factors -----

  if (!is.null(MaxLF) && !is.numeric(MaxLF)) {
    stop("`MaxLF` must be NULL or a numeric value", call. = FALSE)
  }

  if (!is.null(MinLF) && !is.numeric(MinLF)) {
    stop("`MinLF` must be NULL or a numeric value", call. = FALSE)
  }

  # Ensure MinLF >= 1 and < MaxLF
  if (!is.null(MinLF) && MinLF < 1) {
    stop("`MinLF` must be >= 1", call. = FALSE)
  }

  if (!is.null(MinLF) && !is.null(MaxLF) && MinLF >= MaxLF) {
    stop("`MinLF` must be less than `MaxLF`", call. = FALSE)
  }

  # convert MaxLF and MinLF to integer
  if (!is.null(MaxLF)) {
    MaxLF <- as.integer(MaxLF)
  }

  if (!is.null(MinLF)) {
    MinLF <- as.integer(MinLF)
  }

  # # ..................................................................... ###

  # coordinates of the sampling units as sf object -----

  SU_Sf <- tibble::as_tibble(Coords) %>%
    stats::setNames(c("x", "y")) %>%
    sf::st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 3035)

  # # ..................................................................... ###

  # Prepare knots ----
  Knots <- Hmsc::constructKnots(
    sData = Coords, knotDist = MinDist, minKnotDist = MinDist
  ) %>%
    tibble::as_tibble() %>%
    # Hmsc::constructKnots may return duplicated points; discard them
    dplyr::distinct() %>%
    dplyr::mutate(
      Identical = purrr::map2_lgl(
        .x = Var1, .y = Var2,
        .f = ~ {
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

  # # ..................................................................... ###

  # Prepare random level ----

  rL <- Hmsc::HmscRandomLevel(sData = Coords, sMethod = "GPP", sKnot = Knots)

  if (is.null(MinLF) && !is.null(MaxLF)) {
    rL <- Hmsc::setPriors(rL, nfMax = MaxLF)
  }
  if (!is.null(MinLF) && is.null(MaxLF)) {
    rL <- Hmsc::setPriors(rL, nfMin = MinLF)
  }
  if (!is.null(MinLF) && !is.null(MaxLF)) {
    rL <- Hmsc::setPriors(rL, nfMin = MinLF, nfMax = MaxLF)
  }

  # # ..................................................................... ###

  # Set prior for alphapw ----
  if (!is.null(Prior)) {
    rL <- Hmsc::setPriors(rL, alphapw = Prior)
  }

  # # ..................................................................... ###

  # Return the random level object
  return(rL)
}
