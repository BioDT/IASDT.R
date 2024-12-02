## |------------------------------------------------------------------------| #
# InfoChunk ----
## |------------------------------------------------------------------------| #

#' Print Information chunk with time stamp
#'
#'
#' This function prints a formatted message with a timestamp, surrounded by
#' separators for better readability in console outputs or logs.
#' @param Message A character string representing the main message to be
#'   timestamped. This parameter is mandatory and cannot be `NULL` or empty.
#' @param Date logical; whether to include the date in the timestamp. Default is
#'   `FALSE`, meaning only the time is printed. See [IASDT.R::CatTime].
#' @param Extra1,Extra2 integer; the number of extra empty lines to print before
#'   and after the separator lines. See [IASDT.R::CatSep] for more details.
#' @param Time logical; whether to include the time in the timestamp. Default is
#'   `FALSE`.
#' @param ... Additional arguments passed to [IASDT.R::CatSep] for customizing
#'   the separators.
#' @inheritParams CatTime
#' @author Ahmed El-Gabbas
#' @return The function does not return any value but prints the message and
#'   separators to the console.
#' @name InfoChunk
#' @examples
#' InfoChunk(Message = "Started")
#' InfoChunk(Message = "finished", Char = "*", CharReps = 60)
#' @export

InfoChunk <- function(
    Message = "", Date = TRUE, Extra1 = 0L, Extra2 = 1L, Bold = FALSE,
    Red = FALSE, Time = FALSE, Level = 0L, NLines = 1L, ...) {

  if (is.null(Message)) {
    stop("Message cannot be NULL", call. = FALSE)
  }

  IASDT.R::CatSep(
    ..., Extra1 = Extra1, Extra2 = Extra2, Red = Red, Bold = Bold)
  IASDT.R::CatTime(
    Text = Message, NLines = NLines,
    Time = Time, Date = Date, Level = Level, Red = Red, Bold = Bold)
  IASDT.R::CatSep(
    ..., Extra1 = Extra1, Extra2 = Extra2, Red = Red, Bold = Bold)

  return(invisible(NULL))
}

{

InfoChunk(
  paste0("\t", "Making spatial predictions"), Extra1 = 0, Extra2 = 1, Rep = 2,
  Char = "*", CharReps = 70, Red = TRUE, Bold = TRUE, Time = FALSE)

  InfoChunk(
    paste0("\t", "Making spatial predictions"), Extra1 = 0, Extra2 = 1, Rep = 2,
    Char = "*", CharReps = 70, Red = TRUE, Bold = TRUE, Time = FALSE)
cat("A")
}
