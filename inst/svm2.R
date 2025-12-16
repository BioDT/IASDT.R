# Author: Ahmed El-Gabbas
# Date: 2025-12-15
# Based on: `babaknaimi/sdm/inst/methods/sdm/svm.R`
#
# - New method using e1071::svm backend with minimal changes.
# - Enable probability outputs for classification by setting
# `probability = TRUE` when not provided by the user. This matches original
# svm's probability behavior.
#  - Predict function returns positive-class probability when available.
#  - All other arguments are passed via ... without overriding defaults.

methodInfo <- list(
  name = c("svm2", "SVM2", "svm_e1071"),
  packages = "e1071",
  modelTypes = c("pa", "pb", "ab", "n"),
  fitParams = list(
    formula = "standard.formula", data = "sdmDataFrame", v = "sdmVariables"),
  fitSettings = list(kernel = "radial", probability = TRUE),

  fitFunction = function(formula, data, v, ...) {
    x <- sdm:::.getData.sdmMatrix(
      formula, data, normalize = TRUE, frame = v@varInfo$numeric, scale = FALSE)
    y <- sdm:::.getData.sdmY(formula, data)

    # class.weights is used to counter severe class imbalance in PA/PB data.
    #
    # In presence–absence SDMs, absences (0) often vastly outnumber presences
    # (1). Without weighting, the SVM’s loss is dominated by the majority class,
    # leading to poor discrimination (e.g., predicting almost all 0s).

    # Compute counts of absences (n0) and presences (n1)
    n0 <- sum(y == 0, na.rm = TRUE)
    n1 <- sum(y == 1, na.rm = TRUE)

    # Upper bound for weight
    max_weight <- 20

    # - Upweight the minority class so its misclassification cost is comparable
    # to the majority class, using a simple inverse-prevalence rule:
    # minority_weight = n_majority / n_minority
    # - Cap the weight by max_weight (here 20) to avoid numeric instability and
    # overly aggressive rebalancing on extremely imbalanced folds.
    # - If absences are the majority (n0 >= n1), set weight for class "1"
    # (presence) to min(n0 / n1, max_weight), keep class "0" at 1. Otherwise,
    # upweight class "0" similarly when presences are the majority.
    # - Pass class.weights to e1071::svm so the optimization accounts for
    # imbalance while leaving all other defaults unchanged.
    if (n0 >= n1) {
      # More absences
      class.weights <- c("0" = 1, "1" = min(n0 / n1, max_weight))
    } else {
      # More presences
      class.weights <- c("0" = min(n1 / n0, max_weight), "1" = 1)
    }

    e1071::svm(x = x, y = y, scale = TRUE, class.weights = class.weights, ...)
  },
  settingRules = NULL,
  tuneParams = NULL,
  predictParams = list(
    object = "model", formula = "standard.formula", newx = "sdmDataFrame",
    v = "sdmVariables"),
  predictSettings = list(probability = TRUE),
  predictFunction = function(object, formula, newx, v, ...) {
    newx <- sdm:::.getData.sdmMatrix(
      formula, newx, normalize = TRUE,
      frame = v@varInfo$numeric, scale = FALSE)
    predict(object, newx, ...)
  }
)
