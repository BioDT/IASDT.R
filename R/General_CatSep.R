## |------------------------------------------------------------------------| #
# CatSep ----
## |------------------------------------------------------------------------| #

#' Print separator(s) to the console
#'
#' This function prints customizable separator lines to the console, optionally
#' preceded and followed by empty lines. It is useful for improving the
#' readability of console output in R scripts or during interactive sessions.
#' @param Rep integer; the number of separator lines to print. Default is `1`.
#' @param Extra1,Extra2 integer; the number of extra empty lines to print before
#'   and after the separator lines. Default is `0` and `1`, respectively.
#' @param Char character; the character used to construct the separator line.
#'   Default is `"-"`.
#' @param CharReps integer; the number of times the character is repeated to
#'   form a separator line. Default is `50`.
#' @param ... additional arguments to be passed to `cat()`.
#' @name CatSep
#' @inheritParams CatTime
#' @author Ahmed El-Gabbas
#' @return The function is called for its side effect (printing to the console)
#'   and does not return a meaningful value.
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

CatSep <- function(
    Rep = 1L, Extra1 = 0L, Extra2 = 1L, Char = "-",
    CharReps = 50L, Bold = FALSE, Red = FALSE, ...) {

  # Check input arguments
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  NumericArgs <- c("Rep", "Extra1", "Extra2", "CharReps")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = NumericArgs, Type = "numeric")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "Char", Type = "character")

  if (Bold) {
    Char <- crayon::bold(Char)
  }
  if (Red) {
    Char <- crayon::red(Char)
  }

  if (Extra1 > 0) {
    cat(strrep("\n", Extra1), ...)
  }

  cat(paste(rep(Char, CharReps), collapse = ""), ...)

  if (Extra2 > 0) {
    cat(strrep("\n", Extra2), ...)
  }
  return(invisible(NULL))
}
