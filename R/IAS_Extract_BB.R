## |------------------------------------------------------------------------| #
# Extract_BB ----
## |------------------------------------------------------------------------| #

#' Extract a specific column from the output of GBIF standardization
#'
#' Extract a specific column from the output of GBIF standardization
#'
#' @name Extract_BB
#' @param x GBIF results of GBIF standardization `rgbif::name_backbone`
#' @param var Column name to extract
#' @author Ahmed El-Gabbas
#' @return Scientific name of the accepted taxa
#' @export
#' @details
#' Extract a specific column from the output of GBIF
#' @examples
#' rgbif::name_backbone(name = "Helianthus annuus", kingdom = "plants") %>%
#'    Extract_BB(status)
#'
#' # ------------------------
#'
#' c("Helianthus annuus", "Tagetes patula L.") %>%
#'    tibble::tibble(Taxa = .) %>%
#'    dplyr::mutate(
#'       BB = purrr::map(Taxa, rgbif::name_backbone),
#'       status = purrr::map_chr(BB, Extract_BB, status),
#'       SpKey = purrr::map_int(BB, Extract_BB, speciesKey))

Extract_BB <- function(x = ., var) {
  var <- rlang::ensyms(var) %>%
    as.character()
  if (var %in% names(x)) {
    x[[var]]
  } else {
    NA
  }
}
