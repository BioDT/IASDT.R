## |------------------------------------------------------------------------| #
# CheckFeather ----
## |------------------------------------------------------------------------| #

#' Check the Integrity of an `feather` File
#'
#' This function checks if the specified file is an `feather` file and attempts
#' to read it using [arrow::read_feather]. It returns `TRUE` if the file is an
#' `feather` file and contains a non-null object, otherwise `FALSE.`
#'
#' @param File A character string specifying the path to the file to be checked.
#'   This can not be empty. The function will attempt to determine the file type
#'   based on the file extension. If the file extension is not `feather`, the
#'   function will return `FALSE`.
#' @param warning logical. If `TRUE` (default), warnings will be printed if the
#'   file does not exist.
#' @return A logical value: `TRUE` if the file is an `feather` file and contains
#'   a non-null object, `FALSE` otherwise.
#' @name CheckFeather
#' @author Ahmed El-Gabbas
#' @export

CheckFeather <- function(File, warning = TRUE) {

  if (!file.exists(File) || is.null(File) || !nzchar(File)) {
    if (warning) {
      warning(paste0("The provided file does not exist: `", File, "`"))
    }
    return(FALSE)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(File))

  if (Extension == "feather") {

    Obj <- try(arrow::read_feather(File), silent = TRUE)

    if (inherits(Obj, "try-error")) {
      return(FALSE)
    }

    if (exists("Obj") && !is.null(Obj)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    if (warning) {
      warning(
        "Unsupported file type. Please provide a feather file.", call. = FALSE)
    }
    return(FALSE)

  }
}
