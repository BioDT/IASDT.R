## |------------------------------------------------------------------------| #
# SplitDF2Chunks ----
## |------------------------------------------------------------------------| #

#' Split a data.frame into smaller chunks
#'
#' This function divides a data.frame into smaller chunks based on the specified
#' number of rows per chunk (`ChunkSize`) or the specified number of chunks
#' (`NChunks`). If neither is provided, it defaults to splitting the data.frame
#' into  a minimum of 5 chunks or less if the data.frame has fewer than 5 rows.
#' The function ensures that the data is evenly distributed among the chunks as
#' much as possible.
#' @param DF data.frame. The data.frame to be split into chunks.
#' @param ChunkSize integer (optional). The desired number of rows each chunk
#'   should contain. It must be a positive integer and less than the number of
#'   rows in `DF`.
#' @param NChunks integer (optional). The desired number of chunks to split the
#'   data.frame into. It must be a positive integer.
#' @param Prefix Character. A string value that will be used as a prefix for the
#'   names of the chunks. Default is "Chunk".
#' @name SplitDF2Chunks
#' @author Ahmed El-Gabbas
#' @return A list of data.frames, where each data.frame represents a chunk of
#'   the original data.frame. The names of the list elements are constructed
#'   using the `Prefix` parameter followed by an underscore and the chunk number
#'   (e.g., "Chunk_1", "Chunk_2", ...).
#' @export
#' @examples
#' SplitDF2Chunks(mtcars, ChunkSize = 16)
#'
#' # -------------------------------------------
#'
#' SplitDF2Chunks(mtcars, NChunks = 3)
#'
#' # -------------------------------------------
#'
#' SplitDF2Chunks(mtcars, NChunks = 3, Prefix = "T")

SplitDF2Chunks <- function(
    DF = NULL, ChunkSize = NULL, NChunks = NULL, Prefix = "Chunk") {

  if (is.null(DF)) {
    stop("DF cannot be NULL", .call = FALSE)
  }

  if (!is.null(ChunkSize) && (ChunkSize < 1 || !is.numeric(ChunkSize))) {
    stop("ChunkSize must be numeric and larger than 1", .call = FALSE)
  }

  if (is.null(ChunkSize) && is.null(NChunks)) {
    NChunks <- min(5, nrow(DF))
    cat(paste0(crayon::green(
      paste0("ChunkSize and NChunks are not determined by user. ",
             "Defaulting to split into ")), NChunks, " chunks.\n"))
  }

  if (!is.null(ChunkSize) && nrow(DF) <= ChunkSize) {
    stop(paste0("ChunkSize is larger than the number of rows in the data ",
                "frame!\nPlease use a smaller ChunkSize."), .call = FALSE)
  }

  if (is.null(ChunkSize)) {
    ChunkSize <- ceiling(nrow(DF) / NChunks)
  }

  DF <- tibble::as_tibble(DF)
  Out <- split(DF, (seq_len(nrow(DF)) - 1) %/% ChunkSize)
  names(Out) <- paste0(Prefix, "_", seq_along(Out))

  return(Out)
}
