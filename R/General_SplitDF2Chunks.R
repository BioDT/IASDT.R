## |------------------------------------------------------------------------| #
# SplitDF2Chunks ----
## |------------------------------------------------------------------------| #

#' Split data.frame into smaller chunks
#'
#' Split data.frame into smaller chunks
#' @param DF The data frame to split
#' @param ChunkSize Integer. Number of rows per chunk
#' @param NChunks Integer. Number of chunks (only if `ChunkSize` not provided)
#' @param Prefix Prefix
#' @name SplitDF2Chunks
#' @author Ahmed El-Gabbas
#' @return NULL
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
    DF = NULL, ChunkSize = NULL,
    NChunks = NULL, Prefix = "Chunk") {

  if (is.null(DF)) {
    cat(paste0(crayon::red("DF should be a loaded data frame"), "\n"))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }

  DefaultChunkSize <- DefaultNChunks <- FALSE
  if (is.null(NChunks)) {
    NChunks <- min(5, nrow(DF))
    DefaultNChunks <- TRUE
  }

  if (is.null(ChunkSize)) {
    ChunkSize <- ceiling((nrow(DF) / NChunks))
    DefaultChunkSize <- TRUE
  }

  if (any(!is.numeric(ChunkSize), !(ChunkSize > 1))) {
    cat(paste0(crayon::red("ChunkSize should be a numeric vector of length of 1"), "\n"))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }

  if (all(DefaultChunkSize, DefaultNChunks)) {
    cat(paste0(crayon::green("ChunkSize is not determined by user. The default split into 5 chunks is implemented"), "\n"))
  }


  if (nrow(DF) <= ChunkSize) {
    cat(paste0(crayon::red("ChunkSize is larger than the number of rows in the data frame!\nPlease use a smaller ChunkSize."), "\n"))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }

  DF <- tibble::as_tibble(DF)
  Out <- split(DF, (seq_len(nrow(DF)) - 1) %/% ChunkSize)
  names(Out) <- paste0(Prefix, "_", seq_along(Out))
  Out
}
