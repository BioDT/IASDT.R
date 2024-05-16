# |---------------------------------------------------| #
# SplitVector ----
# |---------------------------------------------------| #

#' Split a vector into smaller chunks
#'
#' Split a vector into smaller chunks
#' @param Vector vector to split
#' @param NSplit number of splits
#' @param Prefix prefix string for the name of the resulted list
#' @name SplitVector
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' SplitVector(1:100, 3)
#'
#' # -------------------------------------------
#'
#' SplitVector(1:100, 2, "T")
#'
#' # -------------------------------------------
#'
#' \dontrun{
#' SplitVector(1:100)
#' }
#' @export

SplitVector <- function(Vector = NULL, NSplit = NULL, Prefix = "Chunk") {
  if (inherits(Vector, "NULL") || inherits(NSplit, "NULL")) {
    stop("Vector and NSplit parameters can not be NULL")
  }
  rep(1:NSplit, length.out = length(Vector)) %>%
    sort() %>%
    as.factor() %>%
    split(x = Vector) %>%
    stats::setNames(paste0(Prefix, "_", 1:NSplit))
}
