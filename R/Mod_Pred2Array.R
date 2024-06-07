## |------------------------------------------------------------------------| #
# Pred2Array ----
## |------------------------------------------------------------------------| #

#' Convert results from predict function to array
#'
#' Convert results from predict function to array, with possibility to make predictions first
#'
#' @param Predict Logical. Make predictions first
#' @param Model Hmsc object or path of the saved model
#' @param Preds Pre-computed predictions
#' @param expected Logical. Whether to return the location parameter of the observation models or sample the values from those.
#' @param NCores Integer. Number of parallel computations for computing predicted values.
#' @name Pred2Array
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Pred2Array <- function(
    Predict = TRUE, Model = NULL, Preds = NULL, expected = TRUE, NCores = NULL) {
  if (Predict) {

    if (inherits(Model, "character")) {
      Model <- IASDT.R::LoadAs(Model)
    }

    Preds <- stats::predict(object = Model, expected = expected, nParallel = NCores)
    invisible(gc())
  }

  array(
    unlist(Preds),
    dim = c(dim(Preds[[1]])[1], dim(Preds[[1]])[2], length(Preds))) %>%
    return()
}
