## |------------------------------------------------------------------------| #
# KeepOnly ----
## |------------------------------------------------------------------------| #

#' Keep only specified objects in the environment, removing all others.
#'
#' This function selectively retains the objects specified in the `Obj`
#' parameter in the current environment, removing all other objects. It is
#' useful for memory management by clearing unnecessary objects from the
#' environment. The function also provides an option to print the names of the
#' kept and removed variables.
#' @name KeepOnly
#' @param Obj A character vector specifying the names of the objects to be kept
#'   in the environment.
#' @param Verbose A logical value indicating whether to print the names of kept
#'   and removed variables. Default to `TRUE`.
#' @return No return value, called for side effects.
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' A <- B <- C <- 15
#' ls()
#'
#' KeepOnly("A")
#'
#' ls()
#' rm(list = ls())
#'
#'
#' A <- B <- C <- 15
#' KeepOnly(c("A","B"))
#' ls()

KeepOnly <- function(Obj, Verbose = TRUE) {

  if (is.null(Obj) || length(Obj) == 0) {
    stop("Obj cannot be NULL or empty.", call. = FALSE)
  }

  if (!is.character(Obj)) {
    stop("Obj must be a character vector.", call. = FALSE)
  }

  AllObjects <- ls(pos = parent.frame())
  RemObjects <- setdiff(AllObjects, Obj)

  if (Verbose) {
    cat(crayon::red(
      paste0("Removed Variables (", length(RemObjects), "): ")),
      crayon::blue(paste0(seq_along(RemObjects), ":", RemObjects,
                          collapse = " ||  ")), sep = "")
  }

  rm(list = RemObjects, pos = parent.frame())

  if (Verbose) {
    cat(crayon::red(
      paste0("\nKept Variables (", length(Obj), "): ")),
      crayon::blue(paste0(seq_along(Obj), ":", Obj, collapse = " ||  ")),
      "\n", sep = "")
  }
}
