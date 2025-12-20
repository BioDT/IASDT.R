# Author: Ahmed El-Gabbas Date: 2025-12-17 Licence GPL v3
#
# - gbm2 fits Boosted Regression Trees (BRTs) using `dismo::gbm.step`, which
# selects an optimal number of trees during training based on a step-wise,
# cross-validated procedure.
# - Predictions are produced with `gbm::predict.gbm` at the best iteration
# chosen by gbm.step (object$gbm.call$best.trees).
#
# The original gbm method in sdm package shows a low maximum predicted values

methodInfo <- list(
  name = c("brt2", "BRT2", "gbm2", "GBM2"),
  packages = c("dismo", "gbm"),
  modelTypes = c("pa", "pb", "ab", "n"),
  fitParams = list(
    formula = "standard.formula", data = "sdmDataFrame",
    v = "sdmVariables"),
  fitSettings = list(
    tree.complexity = 2, learning.rate = 0.01, bag.fraction = 0.5,
    n.folds = 5, max.trees = 20000, step.size = 50),

  fitFunction = function(formula, data, v, ...) {
    fam <- switch(
      v@distribution,
      poisson = "poisson",
      multinomial = "multinomial",
      n = "gaussian",
      "bernoulli")

    dismo::gbm.step(
      data = data, gbm.x = all.vars(formula)[-1], gbm.y = all.vars(formula)[1],
      family = fam, plot.main = FALSE, verbose = FALSE, silent = TRUE, ...)
  },
  settingRules = NULL,
  tuneParams = NULL,
  predictParams = list(
    object = "model", formula = "standard.formula",
    newx = "sdmDataFrame", v = "sdmVariables"),
  predictSettings = list(type = "response"),
  predictFunction = function(object, formula, newx, v, type, ...) {
    gbm::predict.gbm(
      object = object, newdata = newx, n.trees = object$gbm.call$best.trees,
      type = type)
  }
)
