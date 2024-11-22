## |------------------------------------------------------------------------| #
# LoadAs ----
## |------------------------------------------------------------------------| #

#' Load objects from `RData` / `qs2` / `rds` / `feather` file
#'
#' This function loads an `RData` file specified by the `File` parameter. If the
#' `RData` file contains a single object, that object is returned directly. If
#' the file contains multiple objects, they are returned as a list with each
#' object accessible by its name. This allows for flexible handling of loaded
#' data without needing to know the names of the objects stored within the RData
#' file ahead of time. The function also supports loading `feather`, `qs2` and
#' `rds` files.
#' @param File character; the file path of the file to be loaded.
#' @param nthreads Number of threads to use when reading `qs2` files. Default 5;
#'   see [qs2::qs_read].
#' @param ... Additional arguments to be passed to the respective load
#'   functions. [base::load] for `RData` files; [qs2::qs_read] for `qs2` files;
#'   [arrow::read_feather] for `feather` files; and [base::readRDS] for `rds`
#'   files.
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

LoadAs <- function(File, nthreads = 5, ...) {

  if (is.null(File)) {
    stop("File cannot be NULL", call. = FALSE)
  }

  Extension <- tools::file_ext(File)

  OutFile <- switch(
    Extension,
    qs2 = qs2::qs_read(file = File, nthreads = nthreads, ...),
    RData = {
      # Load the .RData file and capture the names of loaded objects
      InFile0 <- load(File, ...)

      if (length(InFile0) == 1) {
        OutFile <- get(paste0(InFile0))
      } else {
        OutFile <- lapply(InFile0, function(x) {
          get(paste0(x))
        })
        names(OutFile) <- InFile0
      }
      return(OutFile)
    },
    rds = readRDS(File, ...),
    feather = arrow::read_feather(file = File, ...),
    stop("Unknown file extension", call. = FALSE))

  return(OutFile)
}
