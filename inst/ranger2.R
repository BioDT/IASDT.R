# Author: Ahmed El-Gabbas
# Date: 2025-12-15
# Based on: babaknaimi/sdm inst/methods/sdm/ranger.R
#
# - Convert numeric 0/1 response to factor to ensure classification.

methodInfo <- list(
  name = c("ranger2", "rangerRF2", "rangerForest2"),
  packages = "ranger",
  modelTypes = c("pa", "pb", "ab", "n"),
  fitParams = list(formula = "standard.formula", data = "sdmDataFrame"),
  fitSettings = list(
    num.trees = 1000, mtry = NULL, importance = "none", probability = TRUE,
    quantreg = FALSE, keep.inbag = FALSE, num.threads = 1, verbose = FALSE),

  fitFunction = function(formula, data, ...) {
    df <- as.data.frame(data)
    resp <- all.vars(formula)[1]
    df[[resp]] <- factor(df[[resp]], levels = c(0, 1))
    ranger::ranger(formula = formula, data = df, ...)
  },
  settingRules = NULL,
  tuneParams = NULL,
  predictParams = list(object = "model", data = "sdmDataFrame"),
  predictSettings = list(
    type = "response", num.threads = 1, se.method = "infjack", verbose = FALSE,
    seed = NULL),
  predictFunction = function(object, data, ...) {
    predict(object, data = data, ...)$predictions[, 2]
  }
)
