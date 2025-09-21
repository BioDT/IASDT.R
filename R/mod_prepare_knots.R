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
#' @param min_distance Numeric. Minimum distance between knots in meters. This
#'   distance is used for both `knotDist` and `minKnotDist` parameters of the
#'   [Hmsc::constructKnots] function.
#' @param jitter_distance Numeric. The jitter distance applied to overlapping
#'   coordinates to avoid exact duplicates. Defaults to 100 meters.
#' @param min_lf,max_lf Integer. Minimum and maximum number of latent factors to
#'   be used. Both default to `NULL` which means that the number of latent
#'   factors will be estimated from the data. If either is provided, the
#'   respective values will be used as arguments to [Hmsc::setPriors].
#' @param alphapw Prior for the alpha parameter. Defaults to a list with `Prior
#'   = NULL`, `Min = 10`, `Max = 1500`, and `Samples = 101`. If `alphapw` is
#'   `NULL` or a list with all `NULL` list items, the default prior will be
#'   used. If `Prior` is a matrix, it will be used as the prior. If `Prior =
#'   NULL`, the prior will be generated using the minimum and maximum values of
#'   the alpha parameter (`min` and `max`, respectively; in kilometre) and  the
#'   number of samples (`Samples`). Defaults to a prior with 101 samples ranging
#'   from 10 to 1500 km, with the first value in the second column set to 0.5.
#'
#' @name prepare_knots
#' @author Ahmed El-Gabbas
#' @return An object of class `HmscRandomLevel`, suitable for specifying the
#'   random level in HMSC GPP models. This object contains the prepared knot
#'   locations as a data frame with columns `var_1` and `var_2` (numeric
#'   coordinates), after possible jittering and conversion to avoid overlap.
#' @export

prepare_knots <- function(
    coordinates = NULL, min_distance = NULL, jitter_distance = 100,
    min_lf = NULL, max_lf = NULL,
    alphapw = list(Prior = NULL, Min = 10, Max = 1500, Samples = 101)) {

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  var_1 <- var_2 <- identical <- geometry <- NULL

  # # ..................................................................... ###

  # alphapw Prior values
  alpha_prior <- function(min_val, max_val, n_samples) {
    values <- round(seq(min_val, max_val, length.out = (n_samples - 1)) * 1000)
    prs <- cbind(
      c(0, values),
      c(0.5, rep(0.5 / (n_samples - 1), (n_samples - 1))))
    dimnames(prs) <- NULL
    prs
  }

  if (ncol(coordinates) != 2L) {
    ecokit::stop_ctx(
      "coordinates must be a numeric matrix or data frame with 2 columns",
      ncol_coordinates = ncol(coordinates), include_backtrace = TRUE)
  }

  # # ..................................................................... ###

  # Validate alphapw input ----

  if (is.null(alphapw)) {
    prior <- NULL
  } else {
    if (all(purrr::map_lgl(alphapw, is.null))) {
      prior <- NULL
    } else {

      if (is.null(alphapw$Prior)) {

        if (is.null(alphapw$Min) || !is.numeric(alphapw$Min)) {
          ecokit::stop_ctx(
            "`alphapw$Min` must be numeric", alphapw_min = alphapw$Min,
            include_backtrace = TRUE)
        }

        if (is.null(alphapw$Samples) || !is.numeric(alphapw$Samples)) {
          ecokit::stop_ctx(
            "`alphapw$Samples` must be numeric",
            alphapw_samples = alphapw$Samples, include_backtrace = TRUE)
        }

        if (is.null(alphapw$Max) || !is.numeric(alphapw$Max)) {
          ecokit::stop_ctx(
            "`alphapw$Max` must be numeric", alphapw_max = alphapw$Max,
            include_backtrace = TRUE)
        }

        prior <- alpha_prior(
          min_val = alphapw$Min, max_val = alphapw$Max,
          n_samples = alphapw$Samples)

      } else {

        prior <- alphapw$Prior

        # Check if Prior is a matrix
        if (!inherits(prior, "matrix")) {
          ecokit::stop_ctx(
            "`prior` must be a matrix",
            prior = prior, class_prior = class(prior), include_backtrace = TRUE)
        }

        # Check if prior has 2 columns
        if (ncol(prior) != 2) {
          ecokit::stop_ctx(
            "`prior` should have exactly 2 columns.",
            prior = prior, ncol_prior = ncol(prior), include_backtrace = TRUE)
        }

        # Check if prior has at least 100 rows
        if (nrow(prior) < 100) {
          ecokit::stop_ctx(
            "`prior` must have >= 100 rows.",
            prior = prior, nrow_prior = nrow(prior), include_backtrace = TRUE)
        }

        # Check if the first element of the second column is 0.5
        # if (prior[1, 2] != 0.5) {
        #   ecokit::stop_ctx(
        #      "The first element of the second column is not 0.5.",
        #        prior_1_2 = prior[1, 2], include_backtrace = TRUE)
        # }

        # Check if all values are positive
        if (any(as.vector(prior) < 0)) {
          ecokit::stop_ctx(
            "`prior` can not contain negative numbers",
            sum_negative = sum(as.vector(prior) < 0), include_backtrace = TRUE)
        }

        # Check if the sum of the second column is equal to 1
        if (isFALSE(all.equal(sum(prior[, 2]), 1))) {
          ecokit::stop_ctx(
            paste0(
              "The sum of the second column is not equal to 1 ",
              "(within floating-point tolerance)."),
            sum_prior_2 = sum(prior[, 2]), include_backtrace = TRUE)
        }
      }
    }
  }

  # # ..................................................................... ###

  # Validation for latent factors -----

  if (!is.null(max_lf) && !is.numeric(max_lf)) {
    ecokit::stop_ctx(
      "`max_lf` must be NULL or a numeric value", max_lf = max_lf,
      include_backtrace = TRUE)
  }

  if (!is.null(min_lf) && !is.numeric(min_lf)) {
    ecokit::stop_ctx(
      "`min_lf` must be NULL or a numeric value", min_lf = min_lf,
      include_backtrace = TRUE)
  }

  # Ensure min_lf >= 1 and < max_lf
  if (!is.null(min_lf) && min_lf < 1) {
    ecokit::stop_ctx(
      "`min_lf` must be >= 1", min_lf = min_lf, include_backtrace = TRUE)
  }

  if (!is.null(min_lf) && !is.null(max_lf) && min_lf >= max_lf) {
    ecokit::stop_ctx(
      "`min_lf` must be strictly less than `max_lf`",
      min_lf = min_lf, max_lf = max_lf)
  }

  # convert max_lf and min_lf to integer
  if (!is.null(max_lf)) {
    max_lf <- as.integer(max_lf)
  }

  if (!is.null(min_lf)) {
    min_lf <- as.integer(min_lf)
  }

  # # ..................................................................... ###

  # coordinates of the sampling units as sf object -----

  su_sf <- tibble::as_tibble(coordinates) %>%
    stats::setNames(c("x", "y")) %>%
    sf::st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 3035L)

  # # ..................................................................... ###

  # Prepare knots ----
  knots <- Hmsc::constructKnots(
    sData = coordinates, knotDist = min_distance,
    minKnotDist = min_distance) %>%
    tibble::as_tibble() %>%
    dplyr::rename(var_1 = Var1, var_2 = Var2)
    # Hmsc::constructKnots may return duplicated points; discard them
    dplyr::distinct() %>%
    dplyr::mutate(
      identical = purrr::map2_lgl(
        .x = var_1, .y = var_2,
        .f = ~ {
          ident0 <- dplyr::filter(su_sf, x == .x, y == .y)
          nrow(ident0) > 0
        }))
  rm(su_sf, envir = environment())

  if (any(knots$identical)) {
    knots <- sf::st_as_sf(
      x = knots, coords = c("var_1", "var_2"), remove = FALSE, crs = 3035) %>%
      dplyr::mutate(
        geometry = sf::st_sfc(
          dplyr::if_else(
            identical, sf::st_jitter(geometry, jitter_distance), geometry))) %>%
      sf::st_coordinates() %>%
      as.data.frame() %>%
      stats::setNames(c("var_1", "var_2"))
  } else {
    knots <- as.data.frame(dplyr::select(knots, var_1, var_2))
  }

  # # ..................................................................... ###

  # Prepare random level ----

  rl <- Hmsc::HmscRandomLevel(
    sData = coordinates, sMethod = "GPP", sKnot = knots)

  if (is.null(min_lf) && !is.null(max_lf)) {
    rl <- Hmsc::setPriors(rl, nfMax = max_lf)
  }
  if (!is.null(min_lf) && is.null(max_lf)) {
    rl <- Hmsc::setPriors(rl, nfMin = min_lf)
  }
  if (!is.null(min_lf) && !is.null(max_lf)) {
    rl <- Hmsc::setPriors(rl, nfMin = min_lf, nfMax = max_lf)
  }

  # # ..................................................................... ###

  # Set prior for alphapw ----
  if (!is.null(prior)) {
    rl <- Hmsc::setPriors(rl, alphapw = prior)
  }

  # # ..................................................................... ###

  rl
}
