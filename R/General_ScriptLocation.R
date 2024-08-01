## |------------------------------------------------------------------------| #
# ScriptLocation ----
## |------------------------------------------------------------------------| #
#
#' Retrieve the location of the current R script.
#'
#' This function attempts to find the location of the currently running R
#' script. It first tries to identify the script's location based on the command
#' line arguments used to start the script. If the script is being run in an
#' interactive session within RStudio, it falls back to using the `rstudioapi`
#' to find the file path of the script in the source editor. If the location
#' cannot be determined, it returns `NA`.
#' @return A character string representing the file path of the current R
#'   script, or `NA` if the path cannot be determined.
#' @name ScriptLocation
#' @source The source code of this function was taken from this
#'   [stackoverflow](https://stackoverflow.com/questions/47044068/) question.
#' @importFrom rlang .data
#' @export
#' @examples
#' # ScriptLocation()

ScriptLocation <- function() {
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(
      col = .data$value, into = c("key", "value"), sep = "=",
      fill = "right") %>%
    dplyr::filter(.data$key == "--file") %>%
    dplyr::pull(.data$value)

  if (length(this_file) == 0) {
    this_file <- rstudioapi::getSourceEditorContext()$path
  } else {
    this_file <- NA_character_
  }

  return(this_file)
}
