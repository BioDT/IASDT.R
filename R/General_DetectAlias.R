## |------------------------------------------------------------------------| #
# DetectAlias ----
## |------------------------------------------------------------------------| #
#' Detect aliased variables
#'
#' Detect aliased variables
#'
#' @param DT input data as tibble or data frame
#' @param Verbose print message for the output of the function
#' @name DetectAlias
#' @export
#' @examples
#' if (require("car")) {
#'   x1 <- rnorm(100)
#'   x2 <- 2 * x1
#'   x3 <- rnorm(100)
#'   y <- rnorm(100)
#'
#'   model <- lm(y ~ x1 + x2 + x3)
#'   summary(model)
#'
#'   try(vif(model))
#'   DetectAlias(DT = cbind(x1, x2, x3))
#'   DetectAlias(DT = cbind(x1, x2, x3), Verbose = TRUE)
#'
#'   # excluding x2
#'   model <- lm(y ~ x1 + x3)
#'   summary(model)
#'   try(vif(model))
#' }

DetectAlias <- function(DT, Verbose = FALSE) {

  DT <- cbind.data.frame(XX = rep(1, nrow(DT)), DT)
  form <- paste(names(DT)[-c(1)], collapse = " + ")
  form <- stats::formula(paste("XX", "~", form))
  fit <- stats::lm(form, data = DT)
  Alias <- stats::alias(fit)
  Aliased <- rownames(Alias$Complete)
  if (length(Aliased) > 0) {
    if (Verbose) cat(paste0("Aliased variables: ", Aliased, "\n"))
    return(Aliased)
  } else {
    if (Verbose) cat("No aliased variables found\n")
    return(invisible(NULL))
  }
}
