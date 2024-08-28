## |------------------------------------------------------------------------| #
# CheckRData ----
## |------------------------------------------------------------------------| #

#' Check if a file is a valid RData file
#'
#' This function checks if a given file is a valid RData file by verifying its
#' extension and attempting to load it. It uses the [LoadAs] to load the file.
#' If the file loads successfully and the object exists, it returns `TRUE`. If
#' the file cannot be loaded or the object does not exist, it returns `FALSE`.
#' @param File character. The file path of the file to be checked. This can not
#'   be empty.
#' @return A logical value indicating if the file is a valid RData file
#' @name CheckRData
#' @author Ahmed El-Gabbas
#' @export

CheckRData <- function(File) {

  if (!file.exists(File) || is.null(File) || nzchar(File)) {
    warning(paste0("The provided file does not exist: `", File, "`"))
    return(FALSE)
  }

  if (tools::file_ext(File) == "RData") {
    Obj <- IASDT.R::LoadAs(File)
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
