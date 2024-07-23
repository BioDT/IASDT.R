## |------------------------------------------------------------------------| #
# integer_breaks ----
## |------------------------------------------------------------------------| #
#
#' Integer Breaks for ggplot Axis
#'
#' This function generates a function that calculates integer axis values for ggplot, ensuring that the axis breaks are integers. It is particularly useful for creating cleaner, more readable plots. The function uses `pretty` from the `scales` package to generate a sequence of pretty breaks and then floors these values to ensure they are integers.
#'
#' @param n integer (default: 5), the desired number of breaks on the axis.  Note that the actual number of breaks may slightly differ from what is requested.
#' @param ... additional arguments passed on to `scales::pretty()`.
#' @export
#' @return A function that takes a numeric vector `x` and returns a vector of integer breaks for the axis, with the names attribute set to the break labels.
#' @references [Click here](https://joshuacook.netlify.app/post/integer-values-ggplot-axis/)
#' @examples
#' ggplot2::ggplot(mtcars, ggplot2::aes(x = drat, y = hp)) +
#'   ggplot2::geom_point() +
#'   ggplot2::scale_x_continuous()

#' ggplot2::ggplot(mtcars, ggplot2::aes(x = disp, y = hp)) +
#'   ggplot2::geom_point() +
#'   ggplot2::scale_x_continuous(breaks = integer_breaks(5))

integer_breaks <- function(n = 5, ...) {
  fxn <- function(x) {
    breaks <- floor(pretty(x, n, ...))
    names(breaks) <- attr(breaks, "labels")
    return(breaks)
  }
  return(fxn)
}
