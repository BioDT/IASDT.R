# |---------------------------------------------------| #
# InfoChunk ----
# |---------------------------------------------------| #
#
#' Print Information chunk with time stamp
#'
#' Print Information chunk time stamp
#' @param Message A character string passed to `IASDT.R::CatTime`
#' @param ... additional arguments for the `IASDT.R::CatSep`
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' InfoChunk(Message = "Started")
#'
#' InfoChunk(Message = "finished", Char = "*", CharReps = 60)
#' @export

InfoChunk <- function(Message = "", ...) {
  IASDT.R::CatSep(..., Extra1 = 1)
  IASDT.R::CatTime(Message)
  IASDT.R::CatSep(..., Extra2 = 1)
}
