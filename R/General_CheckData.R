## |------------------------------------------------------------------------| #
# CheckData ----
## |------------------------------------------------------------------------| #

#' Check the Integrity of `RData` / `qs` / `rds` / `feather` Files
#'
#' @param File character. The file path of the file to be checked. This can not
#'   be empty. The function will attempt to determine the file type based on the
#'   file extension. If the file extension is not recognized, the function will
#'   return `FALSE`. The recognized file extensions are `RData` (using
#'   [IASDT.R::CheckRData]), `qs` ([IASDT.R::CheckQs]), `rds`
#'   ([IASDT.R::CheckRDS]), and `feather` ([IASDT.R::CheckFeather]).
#' @param warning logical. If `TRUE` (default), warnings will be printed if the
#'   file does not exist.
#' @param qs_nthreads integer. The number of threads to use when reading the
#'   `qs` file. Default is 5.
#' @return A logical value indicating if the file is a valid RData file
#' @name CheckData
#' @author Ahmed El-Gabbas
#' @export

CheckData <- function(File, warning = TRUE, qs_nthreads = 5) {

  if (!file.exists(File) || is.null(File) || !nzchar(File)) {
    if (warning) {
      warning(paste0("The provided file does not exist: `", File, "`"))
    }
    return(FALSE)
  }

  Extension <- stringr::str_to_lower(tools::file_ext(File))

  OutFile <- switch(
    Extension,
    qs = IASDT.R::CheckQs(File, qs_nthreads = qs_nthreads, warning = warning),
    rdata = IASDT.R::CheckRData(File, warning = warning),
    rds = IASDT.R::CheckRDS(File, warning = warning),
    feather = IASDT.R::CheckFeather(File, warning = warning),
    FALSE)

  return(OutFile)
}
