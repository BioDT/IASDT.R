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
    formula <- as.formula(deparse(formula), env = environment())
    resp <- all.vars(formula)[1]
    data[, resp] <- factor(data[, resp], levels = c(0L, 1L))
    mda::mda(formula = formula, data = data, ...)
  },
  settingRules = NULL,
  tuneParams = list(method = substitute(c(polyreg, mars, gen.ridge, bruto))),
  predictParams = list(object = "model", newdata = "sdmDataFrame"),
  predictSettings = list(type = "posterior"),
  predictFunction = function(object, newdata, type, ...) {
    predict(object, newdata = newdata, type = type, ...)[, "1"]
  }
)
