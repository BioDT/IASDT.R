## |------------------------------------------------------------------------| #
# Package_RemoteSha ----
## |------------------------------------------------------------------------| #

#' get remote sha of R packages
#'
#' get remote sha of R packages
#'
#' @name Package_RemoteSha
#' @param ... name of one or more R packages
#' @export
#' @examples
#' # Package_RemoteSha(IASDT.R, devtools)

Package_RemoteSha <- function(...) {
  Pk <- rlang::ensyms(...) %>%
    as.character()
  Pk %>%
    purrr::map_chr(~{
      pak::lib_status() %>%
        dplyr::filter(package == .x) %>%
        dplyr::pull("remotesha")
    }) %>%
    stats::setNames(Pk)
}
