## |------------------------------------------------------------------------| #
# SaveAs ----
## |------------------------------------------------------------------------| #

#' Save an object to a file with a new name
#'
#' This function saves an R object to a specified file path with a potentially
#' new name. It is useful for renaming objects during the save process.
#'
#' @param InObj The input object to be saved. This can be an actual R object or
#'   a character string representing the name of an object.
#' @param OutObj A character string specifying the new name for the saved
#'   object. This name is used when the object is loaded back into R.
#' @param OutPath A character string specifying the file path (`*.RData`) where
#'   the object should be saved. This includes the directory and the file name.
#' @name SaveAs
#' @author Ahmed El-Gabbas
#' @return The function does not return a value but saves an object to the
#'   specified file path.
#' @export
#' @examples
#' TMP_Folder <- file.path(tempdir(), stringi::stri_rand_strings(1, 5))
#' fs::dir_create(TMP_Folder)
#' list.files(TMP_Folder)
#'
#' # save iris data in `iris2.RData` with `iris2` object name
#' SaveAs(iris, "iris2", file.path(TMP_Folder, "iris2.RData"))
#' list.files(TMP_Folder, pattern = "^.+.RData")
#'
#' (load(file.path(TMP_Folder, "iris2.RData")))
#'
#' tibble::tibble(iris2)

SaveAs <- function(InObj, OutObj, OutPath) {

  if (is.null(InObj) || is.null(OutObj) || is.null(OutPath)) {
    stop("InObj, OutObj, OutPath cannot be NULL", .call = FALSE)
  }

  if (inherits(InObj, "character")) {
    InObj <- get(InObj)
  }

  # Create directory if not available
  fs::dir_create(dirname(OutPath))

  OutObj <- eval(OutObj)
  assign(OutObj, InObj)
  save(list = OutObj, file = OutPath)
}
