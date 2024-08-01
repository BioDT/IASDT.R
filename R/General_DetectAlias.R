## |------------------------------------------------------------------------| #
# DetectAlias ----
## |------------------------------------------------------------------------| #
#'
#' Detect Aliased Variables in a Linear Model
#'
#' This function identifies aliased (linearly dependent) variables in a linear
#' model by adding a constant column to the data frame, fitting a linear model,
#' and then using the alias function to detect aliased variables.
#' @param DT A `data frame` or `tibble` containing the variables to be checked
#'   for aliasing.
#' @param Verbose A logical value indicating whether to print the aliased
#'   variables found (if any). If `TRUE`, aliased variables are printed to the
#'   console. Defaults to `FALSE`.
#' @return Returns a character vector of aliased variable names if any are
#'   found; otherwise, returns `NULL` invisibly. If `Verbose` is `TRUE`, the
#'   function will also print a message to the console.
#' @name DetectAlias
#' @export
#' @examples
#' library("car", warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)
#' x1 <- rnorm(100)
#' x2 <- 2 * x1
#' x3 <- rnorm(100)
#' y <- rnorm(100)
#'
#' model <- lm(y ~ x1 + x2 + x3)
#' summary(model)
#'
#' # there are aliased coefficients in the model
#' try(car::vif(model))
#'
#' # The function identifies the aliased variables
#' DetectAlias(DT = cbind.data.frame(x1, x2, x3))
#'
#' DetectAlias(DT = cbind.data.frame(x1, x2, x3), Verbose = TRUE)
#'
#' # excluding x2 and refit the model
#' model <- lm(y ~ x1 + x3)
#'
#' summary(model)
#'
#' try(car::vif(model))

DetectAlias <- function(DT, Verbose = FALSE) {

  if (is.null(DT)) {
    stop("DT cannot be NULL")
  }

  # Ensure DT is a data.frame or tibble
  if (!is.data.frame(DT)) {
    stop("DT must be a data frame or tibble.")
  }

  # Add a constant column to the data frame
  DT <- cbind.data.frame(XX = rep(1, nrow(DT)), DT)

  # Construct the formula for linear model
  form <- paste(names(DT)[-1], collapse = " + ")
  form <- stats::as.formula(paste("XX", "~", form))

  # Fit the linear model
  fit <- stats::lm(form, data = DT)

  # Detect aliased variables
  Alias <- stats::alias(fit)
  Aliased <- rownames(Alias$Complete)

  # Output aliased variables if any
  if (length(Aliased) > 0) {
    if (Verbose) {
      cat(paste0("Aliased variables: ", paste(Aliased, collapse = ", "), "\n"))
    }
    return(Aliased)
  } else {
    if (Verbose) {
      cat("No aliased variables found\n")
    }
    return(invisible(NULL))
  }
}
