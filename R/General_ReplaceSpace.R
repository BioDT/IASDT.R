# |---------------------------------------------------| #
# ReplaceSpace ----
# |---------------------------------------------------| #
#
#' Replace space with underscore
#'
#' Replace space with underscore
#' @param x string
#' @name ReplaceSpace
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' ReplaceSpace("Genus species")
#'
#' ReplaceSpace("Genus species subspecies")
#' @export

ReplaceSpace <- function(x) {
  stringr::str_replace_all(x, " ", "_")
}
