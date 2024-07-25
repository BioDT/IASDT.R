## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1")  utils::globalVariables(c("."))


#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
#' @param lhs A value or the magrittr placeholder.
#' @param rhs A function call using the magrittr semantics.
#' @return The result of calling `rhs(lhs)`.
NULL


#' @noRd
.onAttach <- function(...) {

  # Retrieve the package version and date dynamically
  packageVersionInfo <- utils::packageVersion("IASDT.R")
  packageDateInfo <- utils::packageDescription("IASDT.R")$Date

  # Display the startup message
  packageStartupMessage(paste0("IASDT.R v", packageVersionInfo, " - Last updated on ", packageDateInfo))

}
