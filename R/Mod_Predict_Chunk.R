## |------------------------------------------------------------------------| #
# Predict_Chunk ----
## |------------------------------------------------------------------------| #

#' Generate predictions for a small chunk of the study area coordinates
#'
#' This function generates predictions for a small chunk of the study area,
#' returning an `sf` object. It is designed to be used internally within the
#' [Predict_Scenario] function and should not be called directly by the user.
#'
#' @param ChunkData A data frame (tibble) containing the data for the current
#'   chunk.
#' @param Model The trained `Hmsc` model object.
#' @param ModelVars Variables used in the model.
#' @return A simple features (`sf`) object containing predictions for each
#'   species as well as species richness. For each species and for species
#'   richness, three maps are produced representing the mean, standard deviation
#'   (sd), and coefficient of variation (cov) of the prediction.
#' @noRd
#' @name Predict_Chunk
#' @author Ahmed El-Gabbas

Predict_Chunk <- function(ChunkData = NULL, Model = NULL, ModelVars = NULL) {

  # # ..................................................................... ###

  if (is.null(ChunkData)) {
    stop("`ChunkData` must be provided and cannot be NULL", call. = FALSE)
  }
  if (is.null(Model)) {
    stop("`Model` must be provided and cannot be NULL", call. = FALSE)
  }
  if (is.null(ModelVars)) {
    stop("`ModelVars` must be provided and cannot be NULL", call. = FALSE)
  }


  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  SR_sd <- SR_mean <- NULL

  # # ..................................................................... ###

  # Coordinates of the current chunk
  XY <- tibble::tibble(x = ChunkData$x, y = ChunkData$y)

  # prepareGradient
  Gradient <- try(
    Hmsc::prepareGradient(
      hM = Model, XDataNew = as.data.frame(ChunkData[, ModelVars]),
      sDataNew = list(sample = as.data.frame(XY))),
    silent = TRUE)

  if (inherits(Gradient, "try-error")) {
    stop("Failed to prepare the gradient.", call. = FALSE)
  }

  # # ..................................................................... ###

  # Predicted probability of occurrence

  Preds <- try(
    stats::predict(
      object = Model, Gradient = Gradient,
      nParallel = 1, expected = TRUE, predictEtaMean = TRUE),
    silent = TRUE)

  if (inherits(Preds, "try-error")) {
    stop("Error in model prediction", call. = FALSE)
  }

  invisible(gc())

  # # ..................................................................... ###

  # Species richness

  SR_DT <- simplify2array(lapply(X = Preds, FUN = rowSums))

  # # ..................................................................... ###

  # Species probability of occurrences

  Sp_DT <- purrr::map_dfr(
    .x = Preds,
    .f = ~ dplyr::bind_cols(XY, as.data.frame(.x))) %>%
    dplyr::summarise(
      dplyr::across(
        .cols = dplyr::everything(),
        .fns = list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ stats::sd(.x, na.rm = TRUE),
          cov = function(x) {
            mean_x <- mean(x, na.rm = TRUE)
            if (mean_x == 0) return(NA_real_)
            stats::sd(x, na.rm = TRUE) / mean_x
          })),
      .by = c("x", "y")) %>%
    dplyr::bind_cols(
      SR_mean = rowMeans(SR_DT, na.rm = TRUE),
      SR_sd = apply(SR_DT, 1, stats::sd)) %>%
    dplyr::mutate(SR_cov = ifelse(SR_mean == 0, NA_real_, SR_sd / SR_mean)) %>%
    dplyr::select(
      "x", "y", tidyselect::starts_with("SR"),
      tidyselect::everything()) %>%
    sf::st_as_sf(coords = c("x", "y"), crs = 3035)

  # # ..................................................................... ###

  return(Sp_DT)
}
