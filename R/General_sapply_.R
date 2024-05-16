# |---------------------------------------------------| #
# sapply_ ----
# |---------------------------------------------------| #

#' Apply a Function over a List or Vector (with no return value)
#'
#' Apply a Function over a List or Vector (with no return value)
#' @param X vector
#' @param FUN function
#' @param simplify should be the output be simplified
#' @param ... additional arguments to sapply
#' @name sapply_
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' par(mfrow = c(1,2), oma = c(0.25, 0.25, 0.25, 0.25), mar = c(3,3,3,1))
#' sapply(list(x = 100:110, y = 110:120), function(V){plot(V, las = 1, main = "sapply")})
#'
#' # -------------------------------------------
#'
#' # nothing returned or printed, only the plotting
#' par(mfrow = c(1,2), oma = c(0.25, 0.25, 0.25, 0.25), mar = c(3,3,3,1))
#' sapply_(list(x = 100:110, y = 110:120), function(V){plot(V, las = 1, main = "sapply_")})
#' @export

sapply_ <- function(X, FUN, simplify = TRUE, ...) {
  invisible(
    sapply(X = X, FUN = FUN, simplify = simplify, ...)
  )
}
