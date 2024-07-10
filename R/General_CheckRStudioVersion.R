## |------------------------------------------------------------------------| #
# CheckRStudioVersion ----
## |------------------------------------------------------------------------| #
#' Check if `RStudio` should be updated
#'
#' Check if `RStudio` should be updated
#'
#' @name CheckRStudioVersion
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' \dontrun{
#' CheckRStudioVersion()
#' }

CheckRStudioVersion <- function() {
  XPath <- ".flex-inhe:nth-child(8)"
  OnlineVersion <- "https://posit.co/download/rstudio-desktop/" %>%
    xml2::read_html() %>%
    rvest::html_node(XPath) %>%
    rvest::html_text2() %>%
    stringr::str_remove_all("RStudio-|.exe") %>%
    stringr::str_replace_all("-", ".")

  InstalledVersion <- rstudioapi::versionInfo() %>%
    "[["("long_version") %>%
    stringr::str_replace_all("\\+", "\\.")

  if (identical(OnlineVersion, InstalledVersion) == FALSE) {
    cat(
      crayon::blue(
        "R-Studio version:",
        crayon::red(crayon::bold(OnlineVersion)),
        "is available.\nInstalled R-studio version:",
        crayon::red(crayon::bold(InstalledVersion)),
        "\nPlease consider updating R-Studio.\n"))
  } else {
    cat(
      crayon::blue(
        "You are using the most recent version of R-Studio: v",
        crayon::red(crayon::bold(InstalledVersion)), ".",
        sep = ""))
  }
}
