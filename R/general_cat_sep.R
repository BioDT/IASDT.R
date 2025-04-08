## |------------------------------------------------------------------------| #
# cat_sep ----
## |------------------------------------------------------------------------| #

#' Print separator(s) to the console
#'
#' This function prints customizable separator lines to the console, optionally
#' preceded and followed by empty lines. It is useful for improving the
#' readability of console output in R scripts or during interactive sessions.
#' @param n_separators integer; the number of separator lines to print. Default
#'   is `1`.
#' @param lines_before,lines_after integer; the number of extra empty lines to
#'   print before and after the separator lines. Default is `0` and `1`,
#'   respectively.
#' @param line_char character; the character used to construct the separator
#'   line. Default is `"-"`.
#' @param repetitions integer; the number of times the character is repeated to
#'   form a separator line. Default is `50`.
#' @param ... additional arguments to be passed to `cat()`.
#' @name cat_sep
#' @inheritParams cat_time
#' @author Ahmed El-Gabbas
#' @return The function is called for its side effect (printing to the console)
#'   and does not return a meaningful value.
#' @examples
#' cat_sep()
#'
#' cat_sep(2)
#'
#' cat_sep(2,2,3)
#'
#' cat_sep(2,2,3, line_char = "*")
#'
#' cat_sep(2,2,3, line_char = "*", repetitions = 20)
#' @export

cat_sep <- function(
    n_separators = 1L, lines_before = 0L, lines_after = 1L, line_char = "-",
    repetitions = 50L, bold = FALSE, red = FALSE, ...) {

  # Check input arguments
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  NumericArgs <- c("n_separators", "lines_before", "lines_after", "repetitions")
  IASDT.R::check_args(
    args_all = AllArgs, args_to_check = NumericArgs, args_type = "numeric")
  IASDT.R::check_args(
    args_all = AllArgs, args_to_check = "line_char", args_type = "character")

  if (bold) {
    line_char <- crayon::bold(line_char)
  }
  if (red) {
    line_char <- crayon::red(line_char)
  }

  if (lines_before > 0) {
    cat(strrep("\n", lines_before), ...)
  }

  cat(paste(rep(line_char, repetitions), collapse = ""), ...)

  if (lines_after > 0) {
    cat(strrep("\n", lines_after), ...)
  }
  return(invisible(NULL))
}
