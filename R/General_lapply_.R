## |------------------------------------------------------------------------| #
# lapply_ ----
## |------------------------------------------------------------------------| #

#' Apply a Function over a List or Vector with Optional Silence
#'
#' This function is a wrapper around the base `lapply` function that allows for the application of a function over a list or vector. It optionally allows for the suppression of the function's return value, making it useful for operations where the user is only interested in the side effects of the function.
#'
#' @param X A list or vector. The input data over which `FUN` is to be applied.
#' @param FUN A function to be applied to each element of `X`.
#' @param Silent A logical value. If `TRUE`, the function suppresses the return value of `FUN` and returns `NULL` invisibly. If `FALSE`, the function returns the result of applying `FUN` over `X`.
#' @param ... Additional arguments to be passed to `FUN`.
#' @name lapply_
#' @author Ahmed El-Gabbas
#' @return If `Silent` is `TRUE`, returns `NULL` invisibly, otherwise returns a list of the same length as `X`, where each element is the result of applying `FUN` to the corresponding element of `X`.
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
  
  if (is.null(X) || is.null(FUN)) {
    stop("X or FUN cannot be NULL")
  }

  result <- lapply(X = X, FUN = FUN, ...)
  
  if (Silent) {
    return(invisible(NULL))
  } else {
    return(result)
  }
}
