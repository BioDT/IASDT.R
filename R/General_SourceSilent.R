# |---------------------------------------------------| #
# SourceSilent ----
# |---------------------------------------------------| #

#' Silently source R script
#'
#' Silently source R script
#'
#' @name SourceSilent
#' @param File String. Path of the file to be sourced
#' @param ... Additional arguments passed to `source` function
#' @param Messages Show messages; default: `TRUE`
#' @param Warnings Show warnings; default: `TRUE`
#' @author Ahmed El-Gabbas
#' @export

SourceSilent <- function(File, Messages = TRUE, Warnings = TRUE, ...) {

  if (Messages && Warnings) {
    File %>%
      source(...) %>%
      utils::capture.output(file = nullfile())
  }

  if (!Messages && !Warnings) {
    File %>%
      source(...) %>%
      utils::capture.output(file = nullfile()) %>%
      suppressMessages() %>%
      suppressWarnings()
  }

  if (Messages && !Warnings) {
    File %>%
      source(...) %>%
      utils::capture.output(file = nullfile()) %>%
      suppressWarnings()
  }

  if (!Messages && Warnings) {
    File %>%
      source(...) %>%
      utils::capture.output(file = nullfile()) %>%
      suppressMessages()
  }
}
