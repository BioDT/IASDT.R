## |------------------------------------------------------------------------| #
# gelman_preplot ----
## |------------------------------------------------------------------------| #

#' gelman_preplot
#'
#'
#' @param x x
#' @param bin.width bin.width
#' @param max.bins max.bins
#' @param confidence confidence
#' @param transform transform
#' @param autoburnin autoburnin
#' @keywords internal
#' @name gelman_preplot
#' @source This function is copied from the unexported `coda:::gelman_preplot`
#'   function. This function is used internally in the [coda::gelman.plot]
#'   function. I make this function available in the package without exporting
#'   it. See: https://svn.r-project.org/R-packages/trunk/coda/R/gelman.R
#' @noRd

gelman_preplot <- function(
    x, bin.width = bin.width, max.bins = max.bins,
    confidence = confidence, transform = transform,
    autoburnin = autoburnin) {

  x <- coda::as.mcmc.list(x)
  nbin <- min(floor((coda::niter(x) - 50) / coda::thin(x)), max.bins)
  if (nbin < 1) {
    IASDT.R::stop_ctx(
      "Insufficient iterations to produce Gelman-Rubin plot", nbin = nbin)
  }
  binw <- floor((coda::niter(x) - 50) / nbin)
  last.iter <- c(
    seq(from = stats::start(x) + 50 * coda::thin(x), by = binw * coda::thin(x),
        length = nbin),
    stats::end(x))
  shrink <- array(dim = c(nbin + 1, coda::nvar(x), 2))
  dimnames(shrink) <- list(
    last.iter, coda::varnames(x),
    c("median", paste0(50 * (confidence + 1), "%")))

  for (i in seq_len(nbin + 1)) {
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
    cat("segments for variables", coda::varnames(x)[all.na], "\n")
    cat("This indicates convergence failure\n")
  }
  return(list(shrink = shrink, last.iter = last.iter))
}
