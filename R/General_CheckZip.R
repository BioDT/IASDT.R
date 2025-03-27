## |------------------------------------------------------------------------| #
# CheckZip ----
## |------------------------------------------------------------------------| #

#' Check the Integrity of a ZIP File
#'
#' Tests the integrity of a ZIP file using the `unzip -t` command. Verifies that
#' the file exists, is non-empty, and has no detectable errors in its compressed
#' data. Returns `FALSE` with a message if the file is invalid or if `unzip` is
#' unavailable.
#'
#' @author Ahmed El-Gabbas
#' @param File Character. The path to the ZIP file to check. Must be a single,
#'   non-empty string.
#' @return Logical: `TRUE` if the file exists, is non-empty, and passes the
#'   integrity check; `FALSE` otherwise.
#' @name CheckZip
#' @export

CheckZip <- function(File) {

  if (isFALSE(IASDT.R::CheckCommands("unzip"))) {
    stop("The 'unzip' command is not available", call. = FALSE)
  }

  if (length(File) != 1 || !inherits(File, "character") || !nzchar(File)) {
    stop("`File` must be a single non-empty character string", call. = FALSE)
  }

  # Verify the file exists
  if (!file.exists(File)) {
    message("File does not exist: ", File)
    return(FALSE)
  }

  # Verify the file is not empty
  if (file.info(File)$size == 0) {
    message("File is empty: ", File)
    return(FALSE)
  }

  # Validate the ZIP file
  FileOkay <- tryCatch(
    expr = {
      IASDT.R::System(stringr::str_glue("unzip -t {File}")) %>%
        stringr::str_detect("No errors detected in compressed data") %>%
        any()
    },
    warning = function(w) {
      message("Warning during file validation: ", conditionMessage(w))
      return(FALSE)
    },
    error = function(e) {
      message("Error during file validation: ", conditionMessage(e))
      return(FALSE)
    })

  # Ensure the result is a logical value
  return(inherits(FileOkay, "logical") && FileOkay && file.info(File)$size > 0)
}
