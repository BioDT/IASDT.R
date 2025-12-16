# Author: Ahmed El-Gabbas
# Date: 2025-12-15
# Based on: babaknaimi/sdm inst/methods/sdm/mda.R
#
# - Enforce classification for PA/PB: convert numeric 0/1 response to factor
# before fitting.

methodInfo <- list(
  name = c("mda2", "MDA2"),
  packages = "mda",
  modelTypes = c("pa", "pb", "n"),
  fitParams = list(formula = "standard.formula", data = "sdmDataFrame"),
  fitSettings = list(method = substitute(polyreg), keep.fitted = FALSE),
  fitFunction = function(formula, data, ...) {
    df <- as.data.frame(data)
    resp <- all.vars(formula)[1]
    df[[resp]] <- factor(df[[resp]], levels = c(0, 1))
    mda::mda(formula = as.formula(formula), data = df, ...)
  },
  settingRules = NULL,
  tuneParams = list(method = substitute(c(polyreg, mars, gen.ridge, bruto))),
  predictParams = list(object = "model", newdata = "sdmDataFrame"),
  predictSettings = list(type = "posterior"),
  predictFunction = function(object, newdata, type, ...) {
    predict(object, newdata = newdata, type = type, ...)[, "1"]
  }
)
