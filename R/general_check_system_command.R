## |------------------------------------------------------------------------| #
# check_system_command ----
## |------------------------------------------------------------------------| #

#' Check system commands availability
#'
#' This function checks if a list of system commands are available on the user's
#' PATH. If any commands are missing, it stops execution and returns an
#' informative error message.
#' @param commands A character vector of system command names to check (e.g.,
#'   `c("git", "Rscript", "unzip")`).
#' @return The function returns `TRUE` if all specified commands are available
#'   on the system, `FALSE` if any is not available.
#' @export
#' @name check_system_command
#' @author Ahmed El-Gabbas
#' @examples
#' check_system_command(c("unzip", "head", "curl", "missing"))
#' @export

check_system_command <- function(commands) {

  # Check the availability of each command in the list
  Okay <- purrr::map_lgl(.x = commands, .f = ~ nzchar(Sys.which(.x)))

  # Identify missing tools
  MissingTools <- commands[!Okay]

  # If any tools are missing, stop with an informative error message
  if (length(MissingTools) > 0) {
    warning(
      "The following tool(s) are missing: ", toString(MissingTools),
      call. = FALSE)
    return(FALSE)
  } else {
    return(TRUE)
  }
}
