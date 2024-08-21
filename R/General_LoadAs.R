## |------------------------------------------------------------------------| #
# LoadAs ----
## |------------------------------------------------------------------------| #
#
#' Load RData file and return its contents as a list or a single object
#'
#' This function loads an `RData` file specified by the `File` parameter. If the
#' `RData` file contains a single object, that object is returned directly. If
#' the file contains multiple objects, they are returned as a list with each
#' object accessible by its name. This allows for flexible handling of loaded
#' data without needing to know the names of the objects stored within the RData
#' file ahead of time.
#' @param File character; the file path of the `.RData` file to be loaded.
#' @author Ahmed El-Gabbas
#' @return Depending on the contents of the `RData` file, this function returns
#'   either a single R object or a named list of R objects. The names of the
#'   list elements (if a list is returned) correspond to the names of the
#'   objects stored within the RData file.
#' @export
#' @name LoadAs
#' @examples
#'
#' File <- system.file("testdata", "culcita_dat.RData", package = "lme4")
#'
#' # ---------------------------------------------------------
#' # loading RData using base library
#' # ---------------------------------------------------------
#' (load(File))
#'
#' ls()
#'
#' tibble::tibble(culcita_dat)
#'
#' # ---------------------------------------------------------
#' # Loading as custom object name
#' # ---------------------------------------------------------
#' NewObj <- LoadAs(File = File)
#'
#' ls()
#'
#' print(tibble::tibble(NewObj))
#'
#' # ---------------------------------------------------------
#' # Loading multiple objects stored in single RData file
#' # ---------------------------------------------------------
#' # store three objects to single RData file
#' mtcars2 <- mtcars3 <- mtcars
#' TempFile <- tempfile(pattern = "mtcars_", fileext = ".RData")
#'
#' save(mtcars2, mtcars3, mtcars, file = TempFile)
#' mtcars_all <- LoadAs(TempFile)
#'
#' # overwrite the file with different order of objects
#' save(mtcars, mtcars2, mtcars3, file = TempFile)
#' mtcars_all2 <- LoadAs(TempFile)
#'
#' # single list object with 3 items, keeping original object names and order
#' names(mtcars_all)
#' names(mtcars_all2)

LoadAs <- function(File) {

  if (is.null(File)) {
    stop("File cannot be NULL", .call = FALSE)
  }

  # Load the .RData file and capture the names of loaded objects
  InFile0 <- load(File)

  if (length(InFile0) == 1) {
    OutFile <- get(paste0(InFile0))
  } else {
    OutFile <- lapply(InFile0, function(x) {
      get(paste0(x))
    })
    names(OutFile) <- InFile0
  }
  return(OutFile)
}
