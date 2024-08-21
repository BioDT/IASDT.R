## |------------------------------------------------------------------------| #
# ReplaceSpace ----
## |------------------------------------------------------------------------| #
#
#' Replace space with underscore in a string
#'
#' This function takes a string as input and replaces all spaces with
#' underscores. It is useful for formatting strings to be used in contexts where
#' spaces are not allowed or desired.
#'
#' @param x A character string. The string in which spaces will be replaced with
#'   underscores.
#' @name ReplaceSpace
#' @author Ahmed El-Gabbas
#' @return A character string with all spaces replaced by underscores.
#' @examples
#' ReplaceSpace("Genus species")
#'
#' ReplaceSpace("Genus species subspecies")
#' @export

ReplaceSpace <- function(x) {
  if (is.null(x)) {
    stop("x name cannot be NULL", .call = FALSE)
  }

  return(stringr::str_replace_all(as.character(x), " ", "_"))
}
