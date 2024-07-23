## |------------------------------------------------------------------------| #
# Package_RemoteSha ----
## |------------------------------------------------------------------------| #

#' Get the remote SHA of specified R packages
#'
#' This function retrieves the remote SHA (Secure Hash Algorithm) reference for one or more specified R packages. It is useful for tracking the exact version of a package's source code.
#'
#' @name Package_RemoteSha
#' @param ... character. Names of one or more R packages for which to retrieve the remote SHA.
#' @return A named character vector where names are the package names and values are the corresponding remote SHAs.
#' @export
#' @examples
#' # Package_RemoteSha(IASDT.R, devtools)

Package_RemoteSha <- function(...) {
  Pk <- rlang::ensyms(...)  %>%
    purrr::map_chr(.f = rlang::as_string)

  purrr::map_chr(
    .x = Pk,
    .f = ~{
      pak::lib_status() %>%
        dplyr::filter(package == .x) %>%
        dplyr::pull("remotesha")
    }) %>%
    stats::setNames(Pk) %>%
    return()
}
