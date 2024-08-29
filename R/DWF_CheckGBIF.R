# # |------------------------------------------------------------------------| #
# Check_GBIF ----
## |------------------------------------------------------------------------| #

#' Check if GBIF login credentials exist
#'
#' @param Renviron Character. The path to the `.Renviron` file containing GBIF
#'   login credentials: `GBIF_EMAIL`, `GBIF_USER`, and `GBIF_PWD`. Defaults to
#'   the current working directory's `.Renviron` file. The credentials must be
#'   in the format:
#'   - `GBIF_EMAIL=your_email`
#'   - `GBIF_USER=your_username`
#'   - `GBIF_PWD=your_password`
#' @name Check_GBIF
#' @author Ahmed El-Gabbas
#' @export

Check_GBIF <- function(Renviron = ".Renviron") {
  # Check if the accessing information already read
  Missing_Account <- purrr::map_lgl(
    .x = c("GBIF_EMAIL", "GBIF_PWD", "GBIF_USER"),
    .f = ~ Sys.getenv(.x) == "") %>%
    any()

  # If access information not read, try to read it from the specified
  # `.Renviron` file
  if (Missing_Account) {
    if (!file.exists(Renviron)) {
      stop(paste0("`.Renviron` file does not exist: ", Renviron), call. = FALSE)
    }
    readRenviron(Renviron)
  }

  Missing_Account <- any(
    purrr::map_lgl(
      .x = c("GBIF_EMAIL", "GBIF_PWD", "GBIF_USER"),
      .f = ~ Sys.getenv(.x) == ""))

  if (Missing_Account) {
    stop(
      paste0(
        "GBIF access information is not available or read from ",
        "the `.Renviron` file"),
      call. = FALSE)
  }
  return(invisible(NULL))
}
