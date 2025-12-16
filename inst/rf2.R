# Author: Ahmed El-Gabbas
# Date: 2025-12-15
# Based on: babaknaimi/sdm inst/methods/sdm/rf.R
#
# - Enforce classification for PA/PB: convert numeric 0/1 response to a
# two-level factor before fitting. This prevents randomForest from running in
# regression mode.
# - Predict probabilities by default (type = "prob") to align with SDM metrics.

methodInfo <- list(
  name = c("rf2", "RF2", "randomForest2"),
  packages = "randomForest",
  modelTypes = c("pa", "pb", "ab", "n"),
  fitParams = list(formula = "standard.formula", data = "sdmDataFrame"),
  fitSettings = list(ntree = 1000, replace = TRUE, importance = TRUE),
  fitFunction = function(formula, data, ...) {
    df <- as.data.frame(data)
    resp <- all.vars(formula)[1]
    df[[resp]] <- factor(df[[resp]], levels = c(0, 1))
    randomForest::randomForest(formula = as.formula(formula), data = df, ...)
  },

  settingRules = NULL,
  tuneParams = NULL,
  predictParams = list(object = "model", newdata = "sdmDataFrame"),
  predictSettings = list(type = "prob"),
  predictFunction = function(object, newdata, type, ...) {
    px <- predict(object, newdata = newdata, type = "prob", ...)
    if (is.matrix(px) && ncol(px) == 2L) {
      return(px[, 2, drop = TRUE])
    }
    # Fallback: if "prob" not available, return numeric prediction
    as.numeric(predict(object, newdata = newdata, ...))
  }
)
