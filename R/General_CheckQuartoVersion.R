## |------------------------------------------------------------------------| #
# CheckQuartoVersion ----
## |------------------------------------------------------------------------| #

#' Check if the installed Quarto version is up to date
#'
#' This function compares the installed Quarto version on the user's system with
#' the latest version available online. If the versions differ, it suggests the
#' user to update Quarto. It uses web scraping to find the latest version
#' available on the Quarto GitHub releases page and the system command to find
#' the installed version.
#' @name CheckQuartoVersion
#' @author Ahmed El-Gabbas
#' @return A message indicating whether the installed Quarto version is up to
#'   date or suggesting an update if it is not.
#' @export
#' @examples
#' CheckQuartoVersion()

CheckQuartoVersion <- function() {

  OnlineVersion <- "https://github.com/quarto-dev/quarto-cli/releases/" %>%
    xml2::read_html() %>%
    rvest::html_nodes(".Link--primary") %>%
    rvest::html_text2() %>%
    stringr::str_remove_all("v") %>%
    gtools::mixedsort(decreasing = TRUE) %>%
    magrittr::extract(1)

  InstalledVersion <- system("quarto --version", intern = TRUE)

  if (identical(OnlineVersion, InstalledVersion) == FALSE) {
    cat(
      crayon::blue(
        "Quarto version:",
        crayon::red(crayon::bold(OnlineVersion)),
        "is available.\nInstalled Quarto version:",
        crayon::red(crayon::bold(InstalledVersion)),
        "\nPlease consider updating Quarto.\n"))
  } else {
    cat(
      crayon::blue(
        "You are using the most recent version of Quarto: v",
        crayon::red(crayon::bold(InstalledVersion)), ".",
        sep = ""))
  }
}
