# # |------------------------------------------------------------------------| #
# GBIF_Check ----
## |------------------------------------------------------------------------| #

#' @author Ahmed El-Gabbas
#' @export
#' @name GBIF_data
#' @rdname GBIF_data
#' @order 2

GBIF_Check <- function(Renviron = ".Renviron") {

  # Check if the accessing information already read
  Missing_Account <- purrr::map_lgl(
    .x = c("GBIF_EMAIL", "GBIF_PWD", "GBIF_USER"),
    .f = ~ !nzchar(Sys.getenv(.x))) %>%
    any()

  # If access information not read, try to read it from the specified
  # `.Renviron` file
  if (Missing_Account) {
    if (!file.exists(Renviron)) {
      stop("`.Renviron` file does not exist: ", Renviron, call. = FALSE)
    }
    readRenviron(Renviron)
  }

  Missing_Account <- any(
    purrr::map_lgl(
      .x = c("GBIF_EMAIL", "GBIF_PWD", "GBIF_USER"),
      .f = ~ !nzchar(Sys.getenv(.x))))

  if (Missing_Account) {
    stop(
      "GBIF access information is not available or read from ",
      "the `.Renviron` file", call. = FALSE)
  }

  return(invisible(NULL))
}
