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
#' @param Prerelease Logical. Whether to check for pre-release versions. Default
#'   is `FALSE`.
#' @return A message indicating whether the installed Quarto version is up to
#'   date or suggesting an update if it is not.
#' @export
#' @examples
#' CheckQuartoVersion()

CheckQuartoVersion <- function(Prerelease = FALSE) {

  OnlineVersions <- "https://github.com/quarto-dev/quarto-cli/releases/" %>%
    xml2::read_html() %>%
    rvest::html_nodes(".Link--primary") %>%
    rvest::html_text2() %>%
    stringr::str_remove_all("v")
  Version_Latest <- OnlineVersions[1]

  # # ..................................................................... ###


  # Check `quarto` system command
  if (IASDT.R::CheckCommands("quarto")) {
    InstalledVersion <- system("quarto --version", intern = TRUE)
  } else {
    cat(crayon::blue("Quarto is not available in the system.\n"))
    InstalledVersion <- NA_character_
  }

  if (isFALSE(identical(Version_Latest, InstalledVersion))) {

    if (Prerelease) {
      Version_PreRelease <- OnlineVersions %>%
        gtools::mixedsort(decreasing = TRUE) %>%
        magrittr::extract(1)
      cat(
        crayon::blue(
          paste0(
            "Available pre-release version is: ",
            crayon::red(crayon::bold(Version_PreRelease)), " [installed: ",
            crayon::red(crayon::bold(InstalledVersion)), "]\n")))
    } else {
      cat(
        crayon::blue(
          paste0(
            "Latest quarto version is ",
            crayon::red(crayon::bold(Version_Latest)), " [installed: ",
            crayon::red(crayon::bold(InstalledVersion)), "]\n")))

    }
  } else {
    cat(
      crayon::blue(
        "You are using the most recent version of Quarto: v",
        crayon::red(crayon::bold(InstalledVersion)), ".",
        sep = ""))
  }

  return(invisible(NULL))
}
