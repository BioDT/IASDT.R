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
#' @param cat_date Logical. Whether to include the date in the timestamp.
#'   Default is `FALSE`, meaning only the time is printed. See
#'   [IASDT.R::cat_time].
#' @param sep_lines_before,sep_lines_after Integer. Number of extra empty lines
#'   to print before and after the separator lines. See [IASDT.R::cat_sep] for
#'   more details.
#' @param cat_timestamp Logical. Whether to include the time in the timestamp.
#'   Default is `FALSE`.
#' @param info_lines_before Integer. Number of extra empty lines to print before
#'   the message. Default is `0L`.
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
#' info_chunk(message = "finished", line_char = "*", line_char_rep = 60)
#'
#' info_chunk(message = "Started", cat_bold =  TRUE, cat_red = TRUE)
#'
#' @export

info_chunk <- function(
    message = "", cat_date = TRUE, sep_lines_before = 0L, sep_lines_after = 1L,
    cat_bold = FALSE, cat_red = FALSE, cat_timestamp = FALSE, level = 0L,
    msg_n_lines = 1L, info_lines_before = 0L, ...) {

  if (is.null(message)) {
    stop("message cannot be NULL", call. = FALSE)
  }

  if (info_lines_before > 1) {
    cat(paste(rep("\n", info_lines_before), collapse = ""), sep = "")
  }

  IASDT.R::cat_sep(
    ..., sep_lines_before = sep_lines_before + 1,
    sep_lines_after = sep_lines_after, cat_red = cat_red, cat_bold = cat_bold)
  IASDT.R::cat_time(
    text = message, msg_n_lines = msg_n_lines,
    cat_timestamp = cat_time, cat_date = cat_date, level = level,
    cat_red = cat_red, cat_bold = cat_bold)
  IASDT.R::cat_sep(
    ..., sep_lines_before = sep_lines_before,
    sep_lines_after = sep_lines_after + 1,
    cat_red = cat_red, cat_bold = cat_bold)

  return(invisible(NULL))
}
