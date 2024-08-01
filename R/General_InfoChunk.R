## |------------------------------------------------------------------------| #
# InfoChunk ----
## |------------------------------------------------------------------------| #
#
#' Print Information chunk with time stamp
#'
#'This function prints a formatted message with a timestamp, surrounded by
#'separators for better readability in console outputs or logs.
#'@param Message A character string representing the main message to be
#'  timestamped.  parameter is mandatory and cannot be `NULL` or empty.
#'@param ... Additional arguments passed to [IASDT.R::CatSep] for customizing
#'  the separators.
#'@author Ahmed El-Gabbas
#'@return The function does not return any value but prints the message and
#'  separators to the console.
#'@name InfoChunk
#' @examples
#' InfoChunk(Message = "Started")
#'
#' InfoChunk(Message = "finished", Char = "*", CharReps = 60)
#' @export

InfoChunk <- function(Message = "", ...) {
  if (is.null(Message)) {
    stop("Message cannot be NULL")
  }

  IASDT.R::CatSep(..., Extra1 = 1)
  IASDT.R::CatTime(Message)
  IASDT.R::CatSep(..., Extra2 = 1)

  return(invisible(NULL))
}
