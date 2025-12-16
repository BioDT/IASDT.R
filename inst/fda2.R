# Author: Ahmed El-Gabbas
# Date: 2025-12-15
# Based on: babaknaimi/sdm inst/methods/sdm/fda.R
#
# - enforce classification for PA/PB: convert numeric 0/1 response to factor
# before fitting.
# - Predict posterior probabilities and return the positive-class column.

methodInfo <- list(
  name = c("fda2", "FDA2"),
  packages = "mda",
  modelTypes = c("pa", "pb", "ab", "n"),
  fitParams = list(formula = "standard.formula", data = "sdmDataFrame"),
  fitSettings = list(method = substitute(polyreg), keep.fitted = FALSE),
  fitFunction = function(formula, data, ...) {
    df <- as.data.frame(data)
    resp <- all.vars(formula)[1]
    df[[resp]] <- factor(df[[resp]], levels = c(0, 1))
    mda::fda(formula = as.formula(formula), data = df, ...)
  },
  settingRules = NULL,
  tuneParams = list(method = substitute(c(polyreg, mars, gen.ridge, bruto))),
  predictParams = list(object = "model", newdata = "sdmDataFrame"),
  predictSettings = list(type = "posterior"),

  predictFunction = function(object, newdata, type, ...) {
    predict(object, newdata = newdata, type = type, ...)[, "1"]
  }
)
