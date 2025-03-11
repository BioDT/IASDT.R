## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1") {
  utils::globalVariables(".")
}

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

#' Null-coalescing operator
#'
#' See \code{rlang::\link[rlang:op-null-default]{\%||\%}} for details.
#'
#' @name %||%
#' @rdname null_default
#' @keywords internal
#' @export
#' @importFrom rlang %||%
#' @usage lhs \%||\% rhs
#' @param lhs Any R object, typically checked for `NULL`.
#' @param rhs A default value to return if `lhs` is `NULL`.
#' @return `lhs` if it is not `NULL`, otherwise `rhs`.
NULL


#' @noRd
.onAttach <- function(...) {

  # Retrieve the package version and date dynamically
  packageVersionInfo <- utils::packageVersion("IASDT.R")
  packageDateInfo <- utils::packageDescription("IASDT.R")$Date

  # Display the startup message
  packageStartupMessage(
    "IASDT.R v", packageVersionInfo,
    " - Last updated on ", packageDateInfo)

}


#' @useDynLib IASDT.R, .registration=TRUE
NULL
