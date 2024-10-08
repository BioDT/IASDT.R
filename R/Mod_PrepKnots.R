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
#' @param MinLF,MaxLF integer. Minimum and maximum number of latent factors to
#'   be used. Both default to `NULL` which means that the number of latent
#'   factors will be estimated from the data. If either is provided, the
#'   respective values will be used as arguments to [Hmsc::setPriors].
#' @param SetAlphapw Logical; whether to set the prior for alphapw using
#'   [PriorAlpha]. Default: `TRUE`.
#' @param AlphapwPar List; parameters for generating prior distribution using
#'   `PriorAlpha` function. Default values are as follows:  `list(n_samples =
#'   300, prob_norm = 0.8, tr_mean = 200, tr_sd = 50, tr_min = 5, tr_max = 400,
#'   seed = NULL)`. See [PriorAlpha] for more details.
#' @name PrepKnots
#' @author Ahmed El-Gabbas
#' @return An object suitable for specifying the random level in HMSC GPP
#'   models. This object contains the prepared knot locations.
#' @export

PrepKnots <- function(
    Coords = NULL, MinDist = NULL, JitterDist = 100,
    MinLF = NULL, MaxLF = NULL, SetAlphapw = TRUE,
    AlphapwPar = list(
      n_samples = 300, prob_norm = 0.8, tr_mean = 200, tr_sd = 50, tr_min = 5,
      tr_max = 400, seed = NULL)) {

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Var1 <- Var2 <- Identical <- geometry <- NULL

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Validate AlphapwPar prior ----

  if (SetAlphapw) {
    if (AlphapwPar$n_samples <= 0) {
      stop("`AlphapwPar$n_samples` must be a positive number.", call. = FALSE)
    }
    if (AlphapwPar$prob_norm < 0 || AlphapwPar$prob_norm > 1) {
      stop("`prob_norm` must be between 0 and 1.", call. = FALSE)
    }
    if (AlphapwPar$tr_mean <= 0) {
      stop("`AlphapwPar$tr_mean` must be a positive number.", call. = FALSE)
    }
    if (AlphapwPar$tr_sd <= 0) {
      stop("`AlphapwPar$tr_sd` must be a positive number.", call. = FALSE)
    }
    if (AlphapwPar$tr_min <= 0) {
      stop("`AlphapwPar$tr_min` must be a positive number.", call. = FALSE)
    }
    if (AlphapwPar$tr_max <= AlphapwPar$tr_min) {
      stop(
        "`AlphapwPar$tr_max` must be greater than `AlphapwPar$tr_min`.",
        call. = FALSE)
    }
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Validation for latent factors -----

  if (!is.null(MaxLF) && !is.numeric(MaxLF)) {
    stop("`MaxLF` must be NULL or a numeric value", call. = FALSE)
  } else {
    if (!is.null(MaxLF)) {
      MaxLF <- as.integer(MaxLF)
    }
  }

  if (!is.null(MinLF) && !is.numeric(MinLF)) {
    stop("`MinLF` must be NULL or a numeric value", call. = FALSE)
  } else {
    if (!is.null(MinLF)) {
      MinLF <- as.integer(MinLF)
    }
  }

  # Ensure MinLF is greater than or equal to 1 and less than MaxLF if both are
  # set
  if (!is.null(MinLF) && MinLF < 1) {
    stop("`MinLF` must be >= 1", call. = FALSE)
  }
  if (!is.null(MinLF) && !is.null(MaxLF) && MinLF >= MaxLF) {
    stop("`MinLF` must be less than `MaxLF`", call. = FALSE)
  }
  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # coordinates of the sampling units as sf object -----

  # nolint start
  SU_Sf <- tibble::as_tibble(Coords) %>%
    stats::setNames(c("x", "y")) %>%
    sf::st_as_sf(coords = c("x", "y"), remove = FALSE, crs = 3035)
  # nolint end

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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
          dplyr::filter(SU_Sf, x == .x, y == .y) %>%
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

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

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

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Set prior for alphapw ----

  if (SetAlphapw) {
    rL <- Hmsc::setPriors(
      rL,
      alphapw = PriorAlpha(
        n_samples = AlphapwPar$n_samples, prob_norm = AlphapwPar$prob_norm,
        tr_mean = AlphapwPar$tr_mean, tr_sd = AlphapwPar$tr_sd,
        tr_min = AlphapwPar$tr_min, tr_max = AlphapwPar$tr_max,
        dist_max = max(rL$alphapw[, 1]) / 1000,
        seed = AlphapwPar$seed))
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Return the random level object

  return(rL)
}



## |------------------------------------------------------------------------| #
# PriorAlpha ----
## |------------------------------------------------------------------------| #

#' Generate a prior distribution using a mix of truncated normal and uniform
#' samples
#'
#' This function generates a prior alpha distribution using a mix of samples
#' from a truncated normal distribution and a uniform distribution.
#'
#' @param n_samples Integer; total number of samples to generate. Default: 300.
#' @param prob_norm Numeric; proportion of samples drawn from the truncated
#'   normal distribution. Default: 0.8.
#' @param tr_mean,tr_sd,tr_min,tr_max Numeric; parameters for the truncated
#'   normal distribution in the kilometer unit. See [truncnorm::rtruncnorm] for
#'   more details.
#' @param dist_max Numeric; maximum value for the uniform distribution samples.
#' @param seed Integer; an optional random seed for reproducibility.
#' @name PriorAlpha
#' @author Ahmed El-Gabbas
#' @export

PriorAlpha <- function(
    n_samples = 300, prob_norm = 0.8, tr_mean = 200, tr_sd = 50, tr_min = 5,
    tr_max = 400, dist_max = NULL, seed = NULL) {

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Validate inputs
  if (n_samples <= 0) {
    stop("`n_samples` must be a positive number.", call. = FALSE)
  }
  if (prob_norm < 0 || prob_norm > 1) {
    stop("`prob_norm` must be between 0 and 1.", call. = FALSE)
  }
  if (tr_mean <= 0) {
    stop("`tr_mean` must be a positive number.", call. = FALSE)
  }
  if (tr_sd <= 0) {
    stop("`tr_sd` must be a positive number.", call. = FALSE)
  }
  if (tr_min <= 0) {
    stop("`tr_min` must be a positive number.", call. = FALSE)
  }
  if (tr_max <= tr_min) {
    stop("`tr_max` must be greater than `tr_min`.", call. = FALSE)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # samples from a truncated normal distribution
  N_Normal <- ceiling(n_samples * prob_norm)
  SamplesNormal <- truncnorm::rtruncnorm(
    N_Normal, a = tr_min, b = tr_max, mean = tr_mean, sd = tr_sd)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


  # samples from a uniform distribution
  N_Uniform <- n_samples - N_Normal - 1
  SamplesUniform <- stats::runif(n = N_Uniform, min = tr_max, max = dist_max)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # combine both sets of samples
  Samples <- 1000 * sort(c(SamplesNormal, SamplesUniform))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # return the samples

  return(
    cbind(c(0, Samples), c(0.5, rep(0.5 / (n_samples - 1), (n_samples - 1)))))
}
