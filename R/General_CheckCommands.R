## |------------------------------------------------------------------------| #
# CheckCommands ----
## |------------------------------------------------------------------------| #

#' Check system commands availability
#'
#' This function checks if a list of system commands are available on the user's
#' PATH. If any commands are missing, it stops execution and returns an
#' informative error message.
#' @param Commands A character vector of system command names to check (e.g.,
#'   `c("git", "Rscript", "unzip")`).
#' @return The function returns `TRUE` if all specified commands are available
#'   on the system, `FALSE` if any is not available.
#' @export
#' @name CheckCommands
#' @author Ahmed El-Gabbas
#' @examples
#' CheckCommands(c("unzip", "head", "curl"))
#' @export

CheckCommands <- function(Commands) {
  # Check the availability of each command in the list
  Okay <- purrr::map_lgl(.x = Commands, .f = ~ nzchar(Sys.which(.x)))

  # Identify missing tools
  MissingTools <- Commands[!Okay]

  # If any tools are missing, stop with an informative error message
  if (length(MissingTools) > 0) {
    warning(
      paste0(
        "The following tool(s) are missing: ",
        paste0(MissingTools, collapse = ", ")),
      call. = FALSE)
    return(FALSE)
  } else {
    return(TRUE)
  }
}
