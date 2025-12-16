# Author: Ahmed El-Gabbas
# Date: 2025-12-15
# Based on: `babaknaimi/sdm/inst/methods/sdm/mars.R`
#
# The same as original mars method, but with increased maxit in glm control to
# 300 to avoid convergence issues in some cases.
#
# Changing `maxit` via the `modelSettings` argument of `sdm::sdm()` does not
# work due to how sdm handles default parameters internally.

methodInfo <- list(
  name = c("mars2", "MARS2", "earth2"),
  packages = "earth",
  modelTypes = c("pa", "pb", "ab", "n"),
  fitParams = list(formula = "standard.formula", data = "sdmDataFrame"),
  fitSettings = list(
    weights = NULL, pmethod = "none", trace = 0,
    glm = list(family = binomial, control = list(maxit = 300))),
  fitFunction = "earth",
  settingRules = function(x = "sdmVariables", f = "fitSettings") {
    if (x@distribution == "poisson") f[["glm"]] <- list(family = poisson)
    list(fitSettings = f)
  },
  tuneParams = list(glm = c(NULL, "default")),
  predictParams = list(object = "model", newdata = "sdmDataFrame"),
  predictSettings = list(type = "response"),
  predictFunction = "predict"
)
