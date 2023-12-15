# |---------------------------------------------------| #
# integer_breaks ----
# |---------------------------------------------------| #
#
#' integer axis values for ggplot
#'
#' integer axis values for ggplot
#'
#' @param n Desired number of breaks. You may get slightly more or fewer breaks that requested.
#' @param ... other arguments passed on to `scales::pretty()`
#' @source https://joshuacook.netlify.app/post/integer-values-ggplot-axis/

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    breaks
  }
  return(fxn)
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
