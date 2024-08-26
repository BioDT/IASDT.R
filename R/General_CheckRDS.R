#' Check if an RDS file can be successfully read
#'
#' This function checks if the specified file is an RDS file and attempts to
#' read it. It returns `TRUE` if the file is an RDS file and contains a non-null
#' object, otherwise `FALSE.` If the file is not an RDS file, the function stops
#' with an error message.
#'
#' @param File A character string specifying the path to the file to be checked.
#'
#' @return A logical value: `TRUE` if the file is an RDS file and contains a
#'   non-null object, `FALSE` otherwise.
#' @name CheckRDS
#' @author Ahmed El-Gabbas
#' @export

CheckRDS <- function(File) {

  if (!file.exists(File)) {
    stop(
      paste0("The provided file: `", File, "` does not exist"), call. = FALSE)
  }

  tryCatch({
    if (tools::file_ext(File) == "rds") {
      Obj <- readRDS(File)
      if (exists("Obj") && !is.null(Obj)) {
        return(TRUE)
      } else {
        return(FALSE)
      }
    } else {
      stop("Unsupported file type. Please provide an rds file.", call. = FALSE)
    }
  }, error = function(e) {
    return(FALSE)
  })
}
