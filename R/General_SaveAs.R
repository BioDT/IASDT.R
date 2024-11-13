## |------------------------------------------------------------------------| #
# SaveAs ----
## |------------------------------------------------------------------------| #

#' Save an object to a file with a new name
#'
#' This function saves an R object to a specified file path with a potentially
#' new name. It is useful for renaming objects during the save process. The
#' function also supports saving objects in the `qs` format.
#'
#' @param InObj The input object to be saved. This can be an actual R object or
#'   a character string representing the name of an object.
#' @param OutObj A character string specifying the new name for the saved
#'   object. This name is used when the object is loaded back into R.
#' @param OutPath A character string specifying the file path (`*.RData`) where
#'   the object should be saved. This includes the directory and the file name.
#' @param qs_preset A character string specifying the preset to use when saving
#'   the object in the `qs` format. The default is "fast". See [qs::qsave].
#' @param ... Additional arguments to be passed to the `save` function.
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

SaveAs <- function(InObj, OutObj, OutPath, qs_preset = "fast", ...) {

  if (is.null(InObj) || is.null(OutPath)) {
    stop("`InObj` and `OutPath` cannot be NULL", call. = FALSE)
  }

  if (inherits(InObj, "character")) {
    InObj <- get(InObj)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(OutPath))

  if (!Extension %in% c("qs", "rdata")) {
    stop(
      "Extension of `OutPath` must be either 'qs' or 'RData'.", .call = FALSE)
  }

  if (Extension == "rdata" && is.null(OutObj)) {
    stop("`OutObj` cannot be NULL for saving RData files", call. = FALSE)
  }

  OutObj <- eval(OutObj)
  assign(OutObj, InObj)

  # Create directory if not available
  fs::dir_create(dirname(OutPath))

  if (Extension == "rdata") {
    save(list = OutObj, file = OutPath, ...)
  } else {
    qs::qsave(x = InObj, file = OutPath, preset = qs_preset)
  }

  return(invisible())
}
