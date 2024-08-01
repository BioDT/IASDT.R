## |------------------------------------------------------------------------| #
# AddLine ----
## |------------------------------------------------------------------------| #
#
#' Add a horizontal or vertical line to the current plot
#'
#' Add a horizontal or vertical line to the current plot
#'
#' @name AddLine
#' @source The source code of this function was taken from this
#'   [stackoverflow](https://stackoverflow.com/questions/27800307/) question.
#' @export
#' @param at Numeric; the relative location of where the line should be plotted.
#'   Cannot be `NULL`.
#' @param Outer Logical; if `TRUE`, the line is plotted outside of the plotting
#'   area. Default is `FALSE`.
#' @param H Logical; if `TRUE` (default), a horizontal line is added. If
#'   `FALSE`, a vertical line is added.
#' @param ... Additional graphical parameters passed to [graphics::abline].
#' @return Invisible; the function is called for its side effect of drawing on
#'   the current plot.
#' @examples
#' # Horizontal line
#' par(oma = c(1, 1, 1, 1), mar = c(3, 3, 1, 1))
#' plot(1:100)
#' AddLine(at = 0.75)
#' AddLine(at = 0.25, Outer = TRUE, lwd = 2)
#' AddLine(at = 0.5, Outer = TRUE, lwd = 2, col = "red")
#'
#' # ---------------------------------------------
#'
#' # Vertical line
#' plot(1:100)
#' AddLine(H = FALSE, at = 0.75)
#' AddLine(H = FALSE, at = 0.25, Outer = TRUE, lwd = 2)
#' AddLine(H = FALSE, at = 0.5, Outer = TRUE, lwd = 2, col = "red")

AddLine <- function(at = NULL, Outer = FALSE, H = TRUE, ...) {

  if (is.null(at)) {
    stop("at cannot be NULL")
  }

  if (Outer) graphics::par(xpd = TRUE)

  if (H) {
    graphics::abline(h = graphics::grconvertY(at, "npc"), ...)
  } else {
    graphics::abline(v = graphics::grconvertX(at, "npc"), ...)
  }

  if (Outer) {
    graphics::par(xpd = FALSE)
  }
}
