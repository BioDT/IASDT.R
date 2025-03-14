## |------------------------------------------------------------------------| #
# Package_RemoteSha ----
## |------------------------------------------------------------------------| #

#' Get the remote SHA of R packages
#'
#' This function retrieves the remote SHA (Secure Hash Algorithm) reference for
#' one or more specified R packages from their remote repositories. The SHA
#' uniquely identifies the version of a package's source code, useful for
#' reproducibility and version tracking.
#'
#' @name Package_RemoteSha
#' @param ... Character. Names of one or more R packages (quoted or unquoted)
#'   for which to retrieve the remote SHA.
#' @return A named character vector where names are the package names and values
#'   are the corresponding remote SHAs. If a package is not found or has no
#'   remote SHA, the value will be `NA`.
#' @export
#' @author Ahmed El-Gabbas
#' @details The function uses `pak::lib_status()` to query the status of
#'   installed packages and extract their remote SHAs. It supports packages
#'   installed from GitHub, GitLab, or other remote sources via `pak`. If a
#'   package is installed from CRAN or locally without a remote SHA, the result
#'   will be `NA`.
#' @examples
#' Package_RemoteSha(Hmsc, IASDT.R, "nonexistent")

Package_RemoteSha <- function(...) {

  # Capture package names as symbols and convert to character strings
  Pk <- rlang::ensyms(...) %>%
    purrr::map_chr(.f = rlang::as_string)

  # Retrieve library status once for efficiency and map over packages
  lib_status <- IASDT.R::AddMissingCols(
    pak::lib_status(), NA_character_, "remotesha")

  Out <- purrr::map_chr(
    .x = Pk,
    .f = ~{

      # Filter library status for the current package
      sha <- dplyr::filter(lib_status, package == .x) %>%
        dplyr::pull("remotesha")

      # Return the SHA if found, otherwise NA
      if (length(sha) > 0) {
        sha
      } else {
        NA_character_
      }

    }) %>%
    # Name the output vector with package names
    stats::setNames(Pk)

  # Return the named vector of SHAs
  return(Out)
}
