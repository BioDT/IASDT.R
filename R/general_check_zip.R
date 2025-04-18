## |------------------------------------------------------------------------| #
# check_zip ----
## |------------------------------------------------------------------------| #

#' Check the Integrity of a ZIP file
#'
#' Tests the integrity of a ZIP file using the `unzip -t` command. Verifies that
#' the file exists, is non-empty, and has no detectable errors in its compressed
#' data. Returns `FALSE` with a message if the file is invalid or if `unzip` is
#' unavailable.
#'
#' @author Ahmed El-Gabbas
#' @param file Character. The path to the ZIP file to check. Must be a single,
#'   non-empty string.
#' @return Logical: `TRUE` if the file exists, is non-empty, and passes the
#'   integrity check; `FALSE` otherwise.
#' @name check_zip
#' @export

check_zip <- function(file) {

  if (isFALSE(IASDT.R::check_system_command("unzip"))) {
    stop("The 'unzip' command is not available", call. = FALSE)
  }

  if (length(file) != 1 || !inherits(file, "character") || !nzchar(file)) {
    stop("`file` must be a single non-empty character string", call. = FALSE)
  }

  # Verify the file exists
  if (!file.exists(file)) {
    message("file does not exist: ", file)
    return(FALSE)
  }

  # Verify the file is not empty
  if (file.info(file)$size == 0) {
    message("file is empty: ", file)
    return(FALSE)
  }

  # Validate the ZIP file
  FileOkay <- tryCatch(
    expr = {
      IASDT.R::system_command(stringr::str_glue("unzip -t {file}")) %>%
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
  return(inherits(FileOkay, "logical") && FileOkay && file.info(file)$size > 0)
}
