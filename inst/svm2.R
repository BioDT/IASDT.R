# Author: Ahmed El-Gabbas
# Date: 2025-12-18
#
# - Binary classification for PA/PB using e1071::svm with probability outputs.
# - Tuning: Parsimonious grid with 5-fold CV via e1071::tune:
#     kernel = "radial"
#     cost  ∈ {1, 5, 10}
#     gamma ∈ {0.01, 0.05, 0.1}
#   Best model (tune.out$best.model) is returned for prediction.
# - Class weights: Inverse-prevalence weighting (capped at 20) to handle
# imbalance:
#     n0 = count of class "0", n1 = count of class "1"
#     weights = c("0" = 1, "1" = min(n0/n1, 20)) if n0 >= n1
#               c("0" = min(n1/n0, 20), "1" = 1) otherwise
# - Prediction: predict(..., probability = TRUE); returns positive-class
# probability:
#   * Assumes response encoded as 0/1 and present as the left-hand side of the formula.
#   * Predictors and their names in new data must match those used at training.
#   * Grid is intentionally small to remain fast and robust models

methodInfo <- list(
  name = c("svm3", "SVM3", "svm_e1071_3"),
  packages = "e1071",
  modelTypes = c("pa", "pb", "ab", "n"),
  fitParams = list(
    formula = "standard.formula", data = "sdmDataFrame", v = "sdmVariables"),
  fitSettings = list(kernel = "radial"),
  fitFunction = function(formula, data, v, ...) {

    formula <- as.formula(deparse(formula), env = environment())
    resp <- all.vars(formula)[1]
    data[, resp] <- factor(data[, resp], levels = c(0L, 1L))

    # Upweight the minority class
    n0 <- sum(data[, resp] == "0")
    n1 <- sum(data[, resp] == "1")
    max_weight <- 20
    if (n0 >= n1) {
      class.weights <- setNames(
        c(1, min(n0 / max(1, n1), max_weight)),
        c("0", "1"))
    } else {
      class.weights <- setNames(
        c(min(n1 / max(1, n0), max_weight), 1),
        c("0", "1"))
    }

    tune.out <- e1071::tune(
      e1071::svm, train.x = formula, data = data, kernel = "radial",
      ranges = list(cost = c(1, 5, 10), gamma = c(0.01, 0.05, 0.1)),
      class.weights = class.weights, probability = TRUE,
      tunecontrol = e1071::tune.control(cross = 5))

    tune.out$best.model
  },
  settingRules = NULL,
  tuneParams = NULL,
  predictParams = list(
    object = "model", formula = "standard.formula", newx = "sdmDataFrame",
    v = "sdmVariables"),
  predictSettings = list(),
  predictFunction = function(object, formula, newx, v, ...) {
    pred_probs <- predict(
      object = object, newdata = newx, probability = TRUE, ...)
    attr(pred_probs, "probabilities")[, "1"]
  }
)
