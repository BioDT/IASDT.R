## |------------------------------------------------------------------------| #
# SplitVector ----
## |------------------------------------------------------------------------| #

#' Split a vector into smaller chunks
#'
#' This function divides a given vector into a specified number of smaller
#' chunks. It is useful for partitioning data into more manageable pieces or for
#' parallel processing tasks.
#' @param Vector A numeric or character vector that you want to split.
#' @param NSplit An integer specifying the number of chunks to split the vector
#'   into. It must not exceed the length of the vector.
#' @param Prefix A string value that serves as the prefix for the names of the
#'   chunks in the returned list. Defaults to "Chunk".
#' @name SplitVector
#' @author Ahmed El-Gabbas
#' @return A list of vectors, where each vector represents a chunk of the
#'   original vector. The names of the list elements are generated using the
#'   specified prefix followed by an underscore and the chunk number.
#' @examples
#' SplitVector(Vector = seq_len(100), NSplit = 3)
#'
#' # -------------------------------------------
#'
#' SplitVector(Vector = seq_len(100), NSplit = 2, Prefix = "T")
#'
#' @export

SplitVector <- function(Vector = NULL, NSplit = NULL, Prefix = "Chunk") {

  if (is.null(Vector) || is.null(NSplit)) {
    stop("Vector and NSplit cannot be NULL", call. = FALSE)
  }

  if (NSplit > length(Vector)) {
    stop("NSplit cannot be greater than the length of Vector", call. = FALSE)
  }

  split(
    Vector,
    cut(seq_along(Vector), breaks = NSplit,
        labels = paste0(Prefix, "_", seq_len(NSplit)),
        include.lowest = TRUE)) %>%
    return()
}
