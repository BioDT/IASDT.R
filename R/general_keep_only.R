## |------------------------------------------------------------------------| #
# keep_only ----
## |------------------------------------------------------------------------| #

#' Keep only specified objects in the environment, removing all others.
#'
#' This function selectively retains the objects specified in the `objects`
#' parameter in the current environment, removing all other objects. It is
#' useful for memory management by clearing unnecessary objects from the
#' environment. The function also provides an option to print the names of the
#' kept and removed variables.
#'
#' @name keep_only
#' @param objects Character vector. Names of the objects to be kept in the
#'   environment.
#' @param verbose Logical. Whether to print the names of kept and removed
#'   variables. Default to `TRUE`.
#' @return No return value, called for side effects.
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' A <- B <- C <- 15
#' ls()
#'
#' keep_only("A")
#'
#' ls()
#' rm(list = ls())
#'
#'
#' A <- B <- C <- 15
#' keep_only(c("A","B"))
#' ls()

keep_only <- function(objects, verbose = TRUE) {

  if (is.null(objects) || length(objects) == 0) {
    IASDT.R::stop_ctx("`objects` cannot be NULL or empty.", objects = objects)
  }

  if (!is.character(objects)) {
    IASDT.R::stop_ctx(
      "`objects` must be a character vector.", objects = objects)
  }

  AllObjects <- ls(pos = parent.frame())
  RemObjects <- setdiff(AllObjects, objects)

  if (verbose) {
    cat(crayon::red(
      paste0("Removed Variables (", length(RemObjects), "): ")),
      crayon::blue(paste0(seq_along(RemObjects), ":", RemObjects,
                          collapse = " ||  ")), sep = "")
  }

  rm(list = RemObjects, pos = parent.frame())

  if (verbose) {
    cat(crayon::red(
      paste0("\nKept Variables (", length(objects), "): ")),
      crayon::blue(
        paste0(seq_along(objects), ":", objects, collapse = " ||  ")),
      "\n", sep = "")
  }
}
