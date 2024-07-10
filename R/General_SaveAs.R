## |------------------------------------------------------------------------| #
# SaveAs ----
## |------------------------------------------------------------------------| #

#' Save an object with a different name
#'
#' Save an object with a different name
#'
#' @param InObj input object
#' @param OutObj output object
#' @param OutPath save path
#' @name SaveAs
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @examples
#' TMP_Folder <- file.path(tempdir(), stringi::stri_rand_strings(1, 5))
#' fs::dir_create(TMP_Folder)
#'
#' # save iris data in `iris2.RData` with `iris2` object name
#' SaveAs(iris, "iris2", file.path(TMP_Folder, "iris2.RData"))
#'
#' list.files(TMP_Folder, pattern = "^.+.RData")
#'
#' (load(file.path(TMP_Folder, "iris2.RData")))
#'
#' tibble::tibble(iris2)

SaveAs <- function(InObj, OutObj, OutPath) {
  if (inherits(InObj, "character")) {
    InObj <- get(InObj)
  }
  OutObj <- eval(OutObj)
  assign(OutObj, InObj)
  save(list = OutObj, file = OutPath)
}
