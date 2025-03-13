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

  Version <- Labels <- NULL

  # URL of the Quarto releases page
  release_blocks <- "https://github.com/quarto-dev/quarto-cli/releases/" %>%
    # Read the HTML content of the page
    rvest::read_html() %>%
    # Extract all release blocks (adjust selector based on current structure)
    rvest::html_nodes("div.Box")

  # Extract version and label for each release
  Releases <- release_blocks %>%
    lapply(function(block) {
      # Get the version number from the <h2><a> or <a> tag
      version <- block %>%
        # Targets release tag links or titles
        rvest::html_node("h2 a, a[href*='/tag/v']") %>%
        rvest::html_text(trim = TRUE)

      # Get all labels (<span> tags with class "Label" or similar)
      labels <- block %>%
        rvest::html_nodes(
          "span.Label, span.Label--success, span.Label--orange") %>%
        rvest::html_text(trim = TRUE) %>%
        unique() %>%
        # Combine multiple labels if present
        paste(collapse = ", ")

      # Return a data frame row
      if (is.na(version)) {
        # Skip if no version found
        return(NULL)
      } else {
        return(
          tibble::tibble(
            Version = version,
            Labels = dplyr::if_else(nzchar(labels), labels, "None")))
      }
    }) %>%
    # Combine into a single data frame
    dplyr::bind_rows() %>%
    dplyr::filter(stringr::str_detect(Version, "^v"))

  Version_Latest <- Releases %>%
    dplyr::filter(stringr::str_detect(Labels, "Latest")) %>%
    dplyr::pull("Version")

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

      Version_PreRelease <- Releases %>%
        dplyr::slice(gtools::mixedorder(Version)) %>%
        dplyr::slice_tail(n = 1) %>%
        dplyr::pull("Version")

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
