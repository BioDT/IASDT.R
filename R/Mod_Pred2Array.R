## |------------------------------------------------------------------------| #
# Mod_Pred2Array ----
## |------------------------------------------------------------------------| #

#' Convert results from `Hmsc:::predict.Hmsc` function to array
#'
#' This function converts the results from the (unexported)
#' `Hmsc:::predict.Hmsc` function into an array format. It allows for the direct
#' conversion of pre-computed predictions or the generation of predictions from
#' a model before conversion. It supports parallel computation for generating
#' predictions.
#' @param Predict Logical; if `TRUE`, predictions are made using the Model
#'   parameter before conversion. If `FALSE`, uses pre-computed predictions
#'   provided in `Preds`.
#' @param Model Hmsc object or character for path to a saved model file. If
#'   `Predict` is `True`, the model object will be used to make predictions
#'   first.
#' @param Preds pre-computed predictions, resulted from `Hmsc:::predict.Hmsc`,
#'   to be converted into an array. Required if `Predict` is `FALSE`.
#' @param expected Logical; if `TRUE`, returns the location parameter of the
#'   observation models. If `FALSE`, samples values from the observation models.
#'   See `Hmsc:::predict.Hmsc` for more info.
#' @param NCores Integer; specifies the number of cores to use for parallel
#'   computation of predicted values. Must be a positive integer.
#' @param predictEtaMean Logical; whether to predict the mean value of the
#'   latent variable.
#' @name Mod_Pred2Array
#' @author Ahmed El-Gabbas
#' @return An array of predictions.
#' @examples
#' Preds <- IASDT.R::Mod_Pred2Array(
#'     Predict = TRUE, Model = Hmsc::TD$m, expected = TRUE, NCores = 1)
#' str(Preds)
#'
#' Preds <- IASDT.R::Mod_Pred2Array(
#'     Predict = TRUE, Model = Hmsc::TD$m, expected = FALSE, NCores = 1)
#' str(Preds)
#' @export

Mod_Pred2Array <- function(
    Model, Predict = TRUE, Preds = NULL, expected = TRUE, NCores = 1,
    predictEtaMean = TRUE) {

  if (Predict && is.null(Model)) {
    stop("Model cannot be empty when Predict is TRUE", call. = FALSE)
  }

  if (isFALSE(Predict) && is.null(Preds)) {
    stop("Preds cannot be empty when Predict is FALSE", call. = FALSE)
  }

  if (Predict) {
    if (!is.numeric(NCores) || NCores <= 0) {
      stop("NCores must be a positive integer.", call. = FALSE)
    }

    if (inherits(Model, "character")) {
      Model <- IASDT.R::LoadAs(Model)
    }

    Preds <- stats::predict(
      object = Model, expected = expected, nParallel = NCores,
      predictEtaMean = predictEtaMean)
  }

  array(
    unlist(Preds),
    dim = c(dim(Preds[[1]])[1], dim(Preds[[1]])[2], length(Preds))) %>%
    return()
}
