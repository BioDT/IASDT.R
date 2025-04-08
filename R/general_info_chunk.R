## |------------------------------------------------------------------------| #
# info_chunk ----
## |------------------------------------------------------------------------| #

#' Print Information chunk with time stamp
#'
#'
#' This function prints a formatted message with a timestamp, surrounded by
#' separators for better readability in console outputs or logs.
#' @param message Character. The main message to be timestamped. This parameter
#'   is mandatory and cannot be `NULL` or empty.
#' @param date Logical. Whether to include the date in the timestamp. Default is
#'   `FALSE`, meaning only the time is printed. See [IASDT.R::cat_time].
#' @param lines_before,lines_after Integer. Number of extra empty lines to print
#'   before and after the separator lines. See [IASDT.R::cat_sep] for more
#'   details.
#' @param time Logical. Whether to include the time in the timestamp. Default is
#'   `FALSE`.
#' @param ... Additional arguments passed to [IASDT.R::cat_sep] for customizing
#'   the separators.
#' @inheritParams cat_time
#' @author Ahmed El-Gabbas
#' @return The function does not return any value but prints the message and
#'   separators to the console.
#' @name info_chunk
#' @examples
#' info_chunk(message = "Started")
#'
#' info_chunk(message = "finished", line_char = "*", repetitions = 60)
#'
#' info_chunk(message = "Started", bold =  TRUE, red = TRUE)
#'
#' @export

info_chunk <- function(
    message = "", date = TRUE, lines_before = 0L, lines_after = 1L,
    bold = FALSE, red = FALSE, time = FALSE, level = 0L, n_lines = 1L, ...) {

  if (is.null(message)) {
    stop("message cannot be NULL", call. = FALSE)
  }

  IASDT.R::cat_sep(
    ..., lines_before = lines_before + 1, lines_after = lines_after,
    red = red, bold = bold)
  IASDT.R::cat_time(
    text = message, n_lines = n_lines,
    time = time, date = date, level = level, red = red, bold = bold)
  IASDT.R::cat_sep(
    ..., lines_before = lines_before, lines_after = lines_after + 1,
    red = red, bold = bold)

  return(invisible(NULL))
}
