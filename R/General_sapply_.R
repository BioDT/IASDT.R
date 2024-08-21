## |------------------------------------------------------------------------| #
# sapply_ ----
## |------------------------------------------------------------------------| #

#' Apply a Function over a List or Vector with Optional Silence
#'
#' This function is a wrapper around the base [base::sapply] function, allowing
#' for the application of a function over a list or vector. It extends
#' [base::sapply] by providing an option to suppress the output, effectively
#' allowing for operations where the user may not care about the return value
#' (e.g., plotting).
#'
#' @param X A list or vector over which `FUN` is to be applied.
#' @param FUN The function to be applied to each element of `X`.
#' @param simplify Logical; should the result be simplified to a vector or
#'   matrix if possible? Passed to [base::sapply].
#' @param Silent Logical; if TRUE, the function returns `invisible(NULL)`
#'   instead of the actual result, effectively suppressing the output. This
#'   enhances the base [base::sapply] for cases where the return value is not
#'   necessary and its output is undesired.
#' @param ... Additional arguments to be passed to `FUN`.
#' @name sapply_
#' @author Ahmed El-Gabbas
#' @return If `Silent` is `FALSE`, returns the result of applying `FUN` over
#'   `X`, potentially simplified according to the `simplify` argument. If
#'   `Silent` is TRUE, returns nothing.
#' @examples
#' par(mfrow = c(1,2), oma = c(0.25, 0.25, 0.25, 0.25), mar = c(3,3,3,1))
#' sapply(
#'     list(x = 100:110, y = 110:120),
#'     function(V) {
#'         plot(V, las = 1, main = "sapply")
#'         })
#'
#' # -------------------------------------------
#'
#' # nothing returned or printed, only the plotting
#' par(mfrow = c(1,2), oma = c(0.25, 0.25, 0.25, 0.25), mar = c(3,3,3,1))
#' sapply_(
#'   list(x = 100:110, y = 110:120),
#'   function(V) {
#'     plot(V, las = 1, main = "sapply_")
#'     })
#'
#' @export

sapply_ <- function(X, FUN, simplify = TRUE, Silent = TRUE, ...) {

  if (is.null(X) || is.null(FUN)) {
    stop("X or FUN cannot be NULL", .call = FALSE)
  }

  result <- sapply(X = X, FUN = FUN, simplify = simplify, ...)

  if (Silent) {
    return(invisible(NULL))
  } else {
    return(result)
  }
}
