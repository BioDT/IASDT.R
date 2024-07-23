## |------------------------------------------------------------------------| #
# LoadAs ----
## |------------------------------------------------------------------------| #
#
#' Load RData file, ignoring original object name
#'
#' Load RData file, ignoring original object name
#'
#' @param File character; path of `.RData` file
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
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
#' # single list object with three items, keeping original object names and order
#' names(mtcars_all)
#' names(mtcars_all2)

LoadAs <- function(File) {

  if (is.null(File)) {
    stop("File cannot be NULL")
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
