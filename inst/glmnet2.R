# Author: Ahmed El-Gabbas
# Date: 2025-12-15
# Based on: `babaknaimi/sdm/inst/methods/sdm/glmnet.R`
#
# - New method name only: "glmnet2", "GLMNET2", "glmelastic2", "glmlasso2".
# - Pass ... into `cv.glmnet` and `glmnet` so users can set alpha, maxit, etc.
# (the original method defines alpha in fitSettings = 1, but does not pass it).
# - In predictFunction, predict at the fitted lambda by setting s =
# object$lambda before selecting the first column; this is safer when multiple
# lambdas could be present.

methodInfo <- list(
  name = c("glmnet2", "GLMNET2", "glmelastic2", "glmlasso2"),
  packages = "glmnet",
  modelTypes = c("pa", "pb", "ab", "n"),
  fitParams = list(
    formula = "standard.formula", data = "sdmDataFrame", v = "sdmVariables"),
  fitSettings = list(family = "binomial", alpha = 1),
  fitFunction = function(formula, data, family, v, ...) {
    x <- sdm:::.getData.sdmMatrix(
      formula, data, normalize = TRUE, frame = v@varInfo$numeric)
    y <- sdm:::.getData.sdmY(formula, data)

    if (family == "binomial") {
      glmnet::cv.glmnet(
        x, y, family = "binomial", type.measure = "auc", ...)
    } else {
      glmnet::cv.glmnet(x, y, family = family, ...)
    }

  },
  settingRules = function(x = "sdmVariables", f = "fitSettings") {
    if (x@distribution %in% c("poisson", "multinomial")) {
      f[["family"]] <- x@distribution
    }
    list(fitSettings = f)
  },
  tuneParams = NULL,
  predictParams = list(
    object = "model", formula = "standard.formula",
    newx = "sdmDataFrame", v = "sdmVariables"),
  predictSettings = list(type = "response"),
  predictFunction = function(object, formula, newx, type, v, ...) {
    newx <- sdm:::.getData.sdmMatrix(
      formula, newx, normalize = TRUE, frame = v@varInfo$numeric)
    # Predict at the fitted lambda for safety (minimal change).
    p <- predict(object, newx, type = type, s = "lambda.1se", ...)
    if (is.matrix(p)) p <- p[, 1, drop = TRUE]
    as.vector(p)
  }
)
