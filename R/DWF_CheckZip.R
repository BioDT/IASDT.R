## |------------------------------------------------------------------------| #
# CheckZip ----
## |------------------------------------------------------------------------| #

#' Check the Integrity of a ZIP File
#'
#' This function checks if a given ZIP file exists and whether it is intact
#' without any errors. The function uses the system's `unzip` command to test
#' the file's integrity.
#'
#' @author Ahmed El-Gabbas
#' @param File Character. The path to the ZIP file that needs to be checked.
#' @return Logical. Returns `TRUE` if the file exists and passes the integrity
#'   check, and `FALSE` otherwise.
#' @name CheckZip
#' @export

CheckZip <- function(File) {

  # Check if the unzip command is available
  IASDT.R::CheckCommands("unzip")

  # Verify the file exists
  if (!file.exists(File)) {
    message(paste0("File does not exist: ", File))
    return(FALSE)
  }

  # Validate the ZIP file
  FileOkay <- tryCatch(
    expr = {
      system2(
        command = "unzip", args = c("-t", File),
        stdout = TRUE, stderr = TRUE) %>%
        stringr::str_detect("No errors detected in compressed data")
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
  return(inherits(FileOkay, "logical") && FileOkay)
}
