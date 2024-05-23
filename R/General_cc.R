## |------------------------------------------------------------------------| #
# cc ----
## |------------------------------------------------------------------------| #

#' Concatenate without quotes
#'
#' Concatenate without quotes
#' @param ... one or more string to concatenate
#' @name cc
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @examples
#' cc(A, B, C)
#'
#' # -------------------------------------------
#'
#' cc(A, B, "A and B")
#' \dontrun{
#' # this does not work
#' cc(A, B, "A and B", 10)
#' # this works
#' cc(A, B, "A and B", "10")
#' }

cc <- function(...) {
  rlang::ensyms(...) %>%
    as.character() %>%
    stringr::str_remove_all("`|`")
}
