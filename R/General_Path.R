## |------------------------------------------------------------------------| #
# Path ----
## |------------------------------------------------------------------------| #

#' Construct path to a file or directory
#'
#' A wrapper around `fs::path()` that constructs file paths from input
#' components and returns them as a character string instead of an `"fs_path"`
#' object.
#'
#' @return A character vector representing the constructed file path(s).
#' @inheritParams fs::path
#' @export
#' @name Path
#' @rdname Path
#' @examples
#' # Basic usage
#' IASDT.R::Path("datasets", "processed", "model_fitting")
#'
#' # Adding a file extension
#' IASDT.R::Path("results", "output", ext = "csv")
#'
#' # Handling multiple components
#' IASDT.R::Path("home", "user", "documents", "report", ext = "pdf")
#'
#' # Using vectorized input
#' IASDT.R::Path("folder", c("file1", "file2"), ext = "txt")

Path <- function(..., ext = "") {
  as.character(fs::path(..., ext = ext))
}
