# |---------------------------------------------------| #
# CheckQuartoVersion ----
# |---------------------------------------------------| #

#' Check if `Quarto` should be updated
#'
#' Check if `Quarto` should be updated
#' @name CheckQuartoVersion
#' @author Ahmed El-Gabbas
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
    "["(1)

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
