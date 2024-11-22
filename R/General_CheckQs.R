## |------------------------------------------------------------------------| #
# CheckQs ----
## |------------------------------------------------------------------------| #

#' Check the Integrity of an `qs` File
#'
#' This function checks if the specified file is an `qs` file and attempts to
#' read it using [qs::qread]. It returns `TRUE` if the file is an `qs` file and
#' contains a non-null object, otherwise `FALSE.`
#'
#' @param File A character string specifying the path to the file to be checked.
#'   This can not be empty. The function will attempt to determine the file type
#'   based on the file extension. If the file extension is not `qs`, the
#'   function will return `FALSE`.
#' @param warning logical. If `TRUE` (default), warnings will be printed if the
#'   file does not exist.
#' @param qs_nthreads integer. The number of threads to use when reading the
#'   `qs` file. Default is 5.
#' @return A logical value: `TRUE` if the file is an qs file and contains a
#'   non-null object, `FALSE` otherwise.
#' @name CheckQs
#' @author Ahmed El-Gabbas
#' @export

CheckQs <- function(File, warning = TRUE, qs_nthreads = 5) {

  if (!file.exists(File) || is.null(File) || !nzchar(File)) {
    if (warning) {
      warning(paste0("The provided file does not exist: `", File, "`"))
    }
    return(FALSE)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(File))

  if (Extension == "qs") {

    Obj <- try(qs::qread(File, nthreads = qs_nthreads), silent = TRUE)

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
      warning("Unsupported file type. Please provide a qs file.", call. = FALSE)
    }
    return(FALSE)

  }
}
