## |------------------------------------------------------------------------| #
# system ----
## |------------------------------------------------------------------------| #

#' Run a system command in a cross-platform manner
#'
#' This function executes a system command, using either `shell` on Windows or
#' `system` on Linux. It allows the output of the command to be captured into an
#' R object.
#' @param command Character. The bash command to be executed.
#' @param R_object Logical. Whether to capture the output of the command as an R
#'   object. If `TRUE` (Default), the output is captured; if `FALSE`, the output
#'   is printed to the console.
#' @param ... Additional arguments passed to either `shell` or `system`
#'   function, depending on the operating system.
#' @name system_command
#' @author Ahmed El-Gabbas
#' @return Depending on the value of `R_object`, either the output of the
#'   executed command as an R object or `NULL` if `R_object` is `FALSE` and the
#'   output is printed to the console.
#' @examples
#' # print working directory
#' IASDT.R::system_command("pwd")
#'
#' # first 5 files on the working directory
#' (A <- IASDT.R::system_command("ls | head -n 5"))
#'
#' (A <- IASDT.R::system_command("ls | head -n 5", R_object = FALSE))
#' @export

system_command <- function(command, R_object = TRUE, ...) {

  # Ensure that command is not NULL
  if (is.null(command)) {
    IASDT.R::stop_ctx("`command` cannot be NULL", command = command)
  }

  if (IASDT.R::OS() == "Windows") {
    Out <- shell(cmd = command, intern = R_object, ...)
  }
  if (IASDT.R::OS() == "Linux") {
    Out <- system(command = command, intern = R_object, ...)
  }

  return(Out)
}
