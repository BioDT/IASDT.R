## |------------------------------------------------------------------------| #
# normalize_path ----
## |------------------------------------------------------------------------| #

#' Normalize and quote file paths
#'
#' This function ensures that file paths are expressed in a consistent and
#' canonical form. It first converts paths to absolute form using
#' `fs::path_abs()`, then tidies them with `fs::path_tidy()`, and finally quotes
#' them correctly based on the operating system. By default, `normalizePath()`
#' behaves differently on Windows and Linux when a file does not exist. On
#' Windows, it tries to construct an absolute path, while on Linux, it returns
#' the input path as-is (relative). To maintain consistency across platforms,
#' this function uses `fs::path_abs()` instead of `normalizePath()`.
#' @param path Character vector. file path(s).
#' @param must_work Logical; if `TRUE`, the function errors for non-existing
#'   paths.
#' @return A character vector of absolute, tidied, and shell-quoted paths.
#' @export
#' @author Ahmed El-Gabbas

normalize_path <- function(path, must_work = FALSE) {

  # Validate input
  if (is.null(path)) {
    stop("Error: `path` cannot be NULL.", call. = FALSE)
  }
  if (!is.character(path)) {
    stop("Error: `path` must be a character vector.", call. = FALSE)
  }
  if (length(path) == 0) {
    stop("Error: `path` cannot be an empty character vector.", call. = FALSE)
  }

  # Check path existence before transformation (if must_work = TRUE)
  if (must_work) {
    Exists <- dplyr::if_else(
      fs::is_dir(path), fs::dir_exists(path), fs::file_exists(path))

    if (isFALSE(Exists)) {
      stop("`path` does not exist: ", path, call. = FALSE)
    }
  }

  # Process and return normalized path
  Out <- fs::path_abs(path) %>%
    fs::path_tidy()

  return(Out)
}
