## |------------------------------------------------------------------------| #
# CheckRData ----
## |------------------------------------------------------------------------| #

#' Check the Integrity of an `RData` File
#'
#' This function checks if a given file is a valid RData file by verifying its
#' extension and attempting to load it. It uses the [LoadAs] to load the file.
#' If the file loads successfully and the object exists, it returns `TRUE`. If
#' the file cannot be loaded or the object does not exist, it returns `FALSE`.
#' @param File character. The file path of the file to be checked. This can not
#'   be empty.
#' @param Warning logical. If `TRUE` (default), warnings will be printed if the
#'   file does not exist.
#' @return A logical value indicating if the file is a valid RData file
#' @name CheckRData
#' @author Ahmed El-Gabbas
#' @export

CheckRData <- function(File, Warning = TRUE) {

  if (!file.exists(File) || is.null(File) || !nzchar(File)) {
    if (Warning) {
      warning(paste0("The provided file does not exist: `", File, "`"))
    }
    return(FALSE)
  }

  if (tools::file_ext(File) == "RData") {

    Obj <- try(IASDT.R::LoadAs(File), silent = TRUE)

    if (inherits(Obj, "try-error")) {
      return(FALSE)
    }

    if (exists("Obj") && !is.null(Obj)) {
      return(TRUE)
    } else {
      return(FALSE)
    }

  } else {
    warning("The provided file is not an RData file")
    return(FALSE)
  }
}
