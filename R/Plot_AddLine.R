## |------------------------------------------------------------------------| #
# AddLine ----
## |------------------------------------------------------------------------| #
#
#' Add a line to the current plot
#'
#' Add a line to the current plot
#'
#' @name AddLine
#' @references [Click here](https://stackoverflow.com/questions/27800307/)
#' @export
#' @param at relative location of where the line should be plotted
#' @param Outer plot the line out of plotting area. Default: `FALSE`
#' @param H Horizontal line? if H == FALSE, vertical line will be added. Default: `TRUE`
#' @param ... Other arguments
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
  if (Outer) graphics::par(xpd = TRUE)

  if (H) {
    graphics::abline(h = graphics::grconvertX(at, "npc"), ...)
  } else {
    graphics::abline(v = graphics::grconvertX(at, "npc"), ...)
  }

  if (Outer) {
    graphics::par(xpd = FALSE)
  }
}
