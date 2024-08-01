## |------------------------------------------------------------------------| #
# gelman.preplot ----
## |------------------------------------------------------------------------| #

#' gelman.preplot
#'
#'
#' @param x x
#' @param bin.width bin.width
#' @param max.bins max.bins
#' @param confidence confidence
#' @param transform transform
#' @param autoburnin autoburnin
#' @keywords internal
#' @name gelman.preplot
#' @source This function is copied from the unexported `coda:::gelman.preplot`
#'   function. This function is used internally in the [coda::gelman.plot]
#'   function. I make this function available in the package without exporting
#'   it. See: https://svn.r-project.org/R-packages/trunk/coda/R/gelman.R
#' @noRd

gelman.preplot <- function(
    x, bin.width = bin.width, max.bins = max.bins,
    confidence = confidence, transform = transform,
    autoburnin = autoburnin) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  niter <- thin <- nvar <- varnames <- NULL


  x <- coda::as.mcmc.list(x)
  nbin <- min(floor((niter(x) - 50) / thin(x)), max.bins)
  if (nbin < 1) {
    stop("Insufficient iterations to produce Gelman-Rubin plot")
  }
  binw <- floor((niter(x) - 50) / nbin)
  last.iter <- c(
    seq(from = stats::start(x) + 50 * thin(x), by = binw * thin(x),
        length = nbin),
    stats::end(x))
  shrink <- array(dim = c(nbin + 1, nvar(x), 2))
  dimnames(shrink) <- list(
    last.iter, varnames(x),
    c("median", paste(50 * (confidence + 1), "%", sep = "")))

  for (i in 1:(nbin + 1)) {
    shrink[i, , ] <- coda::gelman.diag(
      stats::window(x, end = last.iter[i]),
      confidence = confidence,
      transform = transform,
      autoburnin = autoburnin,
      multivariate = FALSE)$psrf
  }
  all.na <- apply(is.na(shrink[, , 1, drop = FALSE]), 2, all)
  if (any(all.na)) {
    cat("\n******* Error: *******\n")
    cat("Cannot compute Gelman & Rubin's diagnostic for any chain \n")
    cat("segments for variables", varnames(x)[all.na], "\n")
    cat("This indicates convergence failure\n")
  }
  return(list(shrink = shrink, last.iter = last.iter))
}
