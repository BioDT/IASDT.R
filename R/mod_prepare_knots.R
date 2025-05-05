## |------------------------------------------------------------------------| #
# prepare_knots ----
## |------------------------------------------------------------------------| #

#' Prepare knot locations for Hmsc GPP models
#'
#' Prepare the locations of knots for use in Gaussian Predictive Process (GPP)
#' models within the HMSC framework. It ensures that knots are spaced at a
#' minimum specified distance and applies jitter to any identical coordinates to
#' avoid overlap.
#' @param coordinates Numeric matrix or data frame containing the (x, y)
#'   coordinates of sampling units.
#' @param min_distance Numeric. Minimum distance between knots
#'   in meters. This distance is used for both `knotDist` and `minKnotDist`
#'   parameters of the [Hmsc::constructKnots] function.
#' @param jitter_distance Numeric. The jitter distance applied to
#'   overlapping coordinates to avoid exact duplicates. Defaults to 100 meters.
#' @param min_LF,max_LF Integer. Minimum and maximum number of latent factors to
#'   be used. Both default to `NULL` which means that the number of latent
#'   factors will be estimated from the data. If either is provided, the
#'   respective values will be used as arguments to [Hmsc::setPriors].
#' @param alphapw Prior for the alpha parameter. Defaults to a list with `Prior
#'   = NULL`, `Min = 20`, `Max = 1200`, and `Samples = 200`. If `alphapw` is
#'   `NULL` or a list with all `NULL` list items, the default prior will be
#'   used. If `Prior` is a matrix, it will be used as the prior. If `Prior =
#'   NULL`, the prior will be generated using `Min`, `Max`, and `Samples`. `Min`
#'   and `Max` are the minimum and maximum values of the alpha parameter (in
#'   kilometre). `Samples` is the number of samples to be used in the prior.
#'
#' @name prepare_knots
#' @author Ahmed El-Gabbas
#' @return An object suitable for specifying the random level in HMSC GPP
#'   models. This object contains the prepared knot locations.
#' @export

prepare_knots <- function(
    coordinates = NULL, min_distance = NULL, jitter_distance = 100,
    min_LF = NULL, max_LF = NULL,
    alphapw = list(Prior = NULL, Min = 20, Max = 1200, Samples = 200)) {
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

  # Validate alphapw input ----

  if (is.null(alphapw)) {
    Prior <- NULL
  } else {
    if (all(purrr::map_lgl(alphapw, is.null))) {
      Prior <- NULL
    } else {

      if (is.null(alphapw$Prior)) {

        if (is.null(alphapw$Min) || !is.numeric(alphapw$Min)) {
          IASDT.R::stop_ctx(
            "`alphapw$Min` must be numeric", alphapw_min = alphapw$Min)
        }

        if (is.null(alphapw$Samples) || !is.numeric(alphapw$Samples)) {
          IASDT.R::stop_ctx(
            "`alphapw$Samples` must be numeric",
            alphapw_samples = alphapw$Samples)
        }

        if (is.null(alphapw$Max) || !is.numeric(alphapw$Max)) {
          IASDT.R::stop_ctx(
            "`alphapw$Max` must be numeric", alphapw_max = alphapw$Max)
        }

        Prior <- AlphaPrior(
          MinVal = alphapw$Min, MaxVal = alphapw$Max,
          NSamples = alphapw$Samples)

      } else {

        Prior <- alphapw$Prior

        # Check if Prior is a matrix
        if (!inherits(Prior, "matrix")) {
          IASDT.R::stop_ctx(
            "`Prior` must be a matrix",
            Prior = Prior, class_prior = class(Prior))
        }

        # Check if Prior has 2 columns
        if (ncol(Prior) != 2) {
          IASDT.R::stop_ctx(
            "`Prior` should have exactly 2 columns.",
            Prior = Prior, ncol_Prior = ncol(Prior))
        }

        # Check if Prior has at least 100 rows
        if (nrow(Prior) < 100) {
          IASDT.R::stop_ctx(
            "`Prior` must have >= 100 rows.",
            Prior = Prior, nrow_Prior = nrow(Prior))
        }

        # Check if the first element of the second column is 0.5
        # if (Prior[1, 2] != 0.5) {
        #   IASDT.R::stop_ctx(
        #      "The first element of the second column is not 0.5.",
        #        prior_1_2 = Prior[1, 2])
        # }

        # Check if all values are positive
        if (any(as.vector(Prior) < 0)) {
          IASDT.R::stop_ctx(
            "`Prior` can not contain negative numbers",
            sum_negative = sum(as.vector(Prior) < 0))
        }

        # Check if the sum of the second column is equal to 1
        if (sum(Prior[, 2]) != 1) {
          IASDT.R::stop_ctx(
            "The sum of the second column is not equal to 1.",
            sum_prior_2 = sum(Prior[, 2]))
        }
      }
    }
  }

  # # ..................................................................... ###

  # Validation for latent factors -----

  if (!is.null(max_LF) && !is.numeric(max_LF)) {
    IASDT.R::stop_ctx(
      "`max_LF` must be NULL or a numeric value", max_LF = max_LF)
  }

  if (!is.null(min_LF) && !is.numeric(min_LF)) {
    IASDT.R::stop_ctx(
      "`min_LF` must be NULL or a numeric value", min_LF = min_LF)
  }

  # Ensure min_LF >= 1 and < max_LF
  if (!is.null(min_LF) && min_LF < 1) {
    IASDT.R::stop_ctx("`min_LF` must be >= 1", min_LF = min_LF)
  }

  if (!is.null(min_LF) && !is.null(max_LF) && min_LF >= max_LF) {
    IASDT.R::stop_ctx(
      "`min_LF` must be less than `max_LF`",
      min_LF = min_LF, max_LF = max_LF)
  }

  # convert max_LF and min_LF to integer
  if (!is.null(max_LF)) {
    max_LF <- as.integer(max_LF)
  }

  if (!is.null(min_LF)) {
    min_LF <- as.integer(min_LF)
  }

  # # ..................................................................... ###

  # coordinates of the sampling units as sf object -----

  SU_Sf <- tibble::as_tibble(coordinates) %>%
    stats::setNames(c("x", "y")) %>%
    sf::st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 3035)

  # # ..................................................................... ###

  # Prepare knots ----
  Knots <- Hmsc::constructKnots(
    sData = coordinates, knotDist = min_distance,
    minKnotDist = min_distance) %>%
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
            Identical, sf::st_jitter(geometry, jitter_distance), geometry))) %>%
      sf::st_coordinates() %>%
      as.data.frame() %>%
      stats::setNames(c("Var1", "Var2"))
  } else {
    Knots <- as.data.frame(dplyr::select(Knots, Var1, Var2))
  }

  # # ..................................................................... ###

  # Prepare random level ----

  rL <- Hmsc::HmscRandomLevel(
    sData = coordinates, sMethod = "GPP", sKnot = Knots)

  if (is.null(min_LF) && !is.null(max_LF)) {
    rL <- Hmsc::setPriors(rL, nfMax = max_LF)
  }
  if (!is.null(min_LF) && is.null(max_LF)) {
    rL <- Hmsc::setPriors(rL, nfMin = min_LF)
  }
  if (!is.null(min_LF) && !is.null(max_LF)) {
    rL <- Hmsc::setPriors(rL, nfMin = min_LF, nfMax = max_LF)
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
