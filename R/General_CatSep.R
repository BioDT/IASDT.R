## |------------------------------------------------------------------------| #
# CatSep ----
## |------------------------------------------------------------------------| #
#
#' Print separator(s) to the console
#'
#' Print separator(s) to the console
#'
#' @param Rep integer; number of separator lines; default `1` row
#' @param Extra1 integer; number of extra empty lines before the separator; default: `0`
#' @param Extra2 integer; number of extra empty lines after the separator; default: `0`
#' @param Char character; the character to be used as a separator; default "-"
#' @param CharReps integer; number of times the character is repeated; default: 50
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' CatSep()
#'
#' CatSep(2)
#'
#' CatSep(2,2,3)
#'
#' CatSep(2,2,3, Char = "*")
#'
#' CatSep(2,2,3, Char = "*", CharReps = 20)
#' @export

CatSep <- function(Rep = 1, Extra1 = 0, Extra2 = 0, Char = "-", CharReps = 50) {
  if (Extra1 > 0) {
    replicate(n = Extra1, expr = cat("\n"))
  }
  S <- c(rep(Char, CharReps)) %>%
    paste0(collapse = "")
  replicate(n = Rep, expr = cat(S, sep = "\n"))
  if (Extra2 > 0) {
    replicate(n = Extra2, expr = cat("\n"))
  }
  return(invisible(NULL))
}
