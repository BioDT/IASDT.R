## |------------------------------------------------------------------------| #
# NormalizePath ----
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
#' @param Path A character vector of file paths.
#' @param MustWork Logical; if `TRUE`, the function errors for non-existing
#'   paths.
#' @return A character vector of absolute, tidied, and shell-quoted paths.
#' @export
#' @author Ahmed El-Gabbas

NormalizePath <- function(Path, MustWork = FALSE) {

  # Validate input
  if (is.null(Path)) {
    stop("Error: 'Path' cannot be NULL.", call. = FALSE)
  }
  if (!is.character(Path)) {
    stop("Error: 'Path' must be a character vector.", call. = FALSE)
  }
  if (length(Path) == 0) {
    stop("Error: 'Path' cannot be an empty character vector.", call. = FALSE)
  }


  # Check path existence before transformation (if MustWork = TRUE)
  if (MustWork) {
    Exists <- dplyr::if_else(
      fs::is_dir(Path), fs::dir_exists(Path),
      fs::file_exists(Path))

    if (isFALSE(Exists)) {
      stop(
        paste0("Path does not exist: ", Path), call. = FALSE)
    }
  }

  # Determine shell type
  shell_type <- dplyr::if_else(
    nzchar(Sys.getenv("ComSpec", unset = "")), "cmd" , "sh")

  # Process and return normalized path
  fs::path_abs(Path) %>%
    fs::path_tidy() %>%
    shQuote(type = shell_type) %>%
    return()
}
