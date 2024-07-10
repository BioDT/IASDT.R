## |------------------------------------------------------------------------| #
# lapply_ ----
## |------------------------------------------------------------------------| #

#' Apply a Function over a List or Vector (with no return value)
#'
#' Apply a Function over a List or Vector (with no return value)
#'
#' @param X vector
#' @param FUN function
#' @param Silent should the output be silenced?
#' @param ... additional arguments to lapply
#' @name lapply_
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' par(mfrow = c(1,2), oma = c(0.25, 0.25, 0.25, 0.25), mar = c(3,3,3,1))
#' lapply(list(x = 100:110, y = 110:120), function(V) {
#'     plot(V, las = 1, main = "lapply")
#' })
#'
#' # -------------------------------------------
#'
#' par(mfrow = c(1,2), oma = c(0.25, 0.25, 0.25, 0.25), mar = c(3,3,3,1))
#' lapply_(list(x = 100:110, y = 110:120), function(V) {
#'     plot(V, las = 1, main = "lapply_")
#' })
#' @export

lapply_ <- function(X, FUN, Silent = TRUE, ...) {
  if (Silent) {
    invisible(lapply(X = X, FUN = FUN, ...))
  } else {
    lapply(X = X, FUN = FUN, ...)
  }
}
