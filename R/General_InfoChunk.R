## |------------------------------------------------------------------------| #
# InfoChunk ----
## |------------------------------------------------------------------------| #

#' Print Information chunk with time stamp
#'
#'
#' This function prints a formatted message with a timestamp, surrounded by
#' separators for better readability in console outputs or logs.
#' @param Message Character. The main message to be timestamped. This parameter
#'   is mandatory and cannot be `NULL` or empty.
#' @param Date Logical. Whether to include the date in the timestamp.
#'   Default is `FALSE`, meaning only the time is printed. See
#'   [IASDT.R::CatTime].
#' @param Extra1,Extra2 Integer. Number of extra empty lines to print before
#'   and after the separator lines. See [IASDT.R::CatSep] for more details.
#' @param Time Logical. Whether to include the time in the timestamp.
#'   Default is `FALSE`.
#' @param ... Additional arguments passed to [IASDT.R::CatSep] for customizing
#'   the separators.
#' @inheritParams CatTime
#' @author Ahmed El-Gabbas
#' @return The function does not return any value but prints the message and
#'   separators to the console.
#' @name InfoChunk
#' @examples
#' InfoChunk(Message = "Started")
#'
#' InfoChunk(Message = "finished", Char = "*", CharReps = 60)
#'
#' InfoChunk(Message = "Started", Bold =  TRUE, Red = TRUE)
#'
#' @export

InfoChunk <- function(
    Message = "", Date = TRUE, Extra1 = 0L, Extra2 = 1L, Bold = FALSE,
    Red = FALSE, Time = FALSE, Level = 0L, NLines = 1L, ...) {

  if (is.null(Message)) {
    stop("Message cannot be NULL", call. = FALSE)
  }

  IASDT.R::CatSep(
    ..., Extra1 = Extra1 + 1, Extra2 = Extra2, Red = Red, Bold = Bold)
  IASDT.R::CatTime(
    Text = Message, NLines = NLines,
    Time = Time, Date = Date, Level = Level, Red = Red, Bold = Bold)
  IASDT.R::CatSep(
    ..., Extra1 = Extra1, Extra2 = Extra2 + 1, Red = Red, Bold = Bold)

  return(invisible(NULL))
}
