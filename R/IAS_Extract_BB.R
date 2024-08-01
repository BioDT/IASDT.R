## |------------------------------------------------------------------------| #
# Extract_BB ----
## |------------------------------------------------------------------------| #

#' Extract a specific column from the output of GBIF standardization
#'
#' This function is designed to extract a specific column from the output of
#' GBIF standardization, typically obtained via the [rgbif::name_backbone]
#' function. It is useful for retrieving specific information such as the
#' scientific name of the accepted taxa.
#' @name Extract_BB
#' @param x A list or data frame containing the results of GBIF standardization
#'   ([rgbif::name_backbone]).
#' @param var A character string specifying the name of the column to extract
#'   from `x`.
#' @author Ahmed El-Gabbas
#' @return Returns the content of the specified column `var` from the input `x`.
#'   If `var` is not found within `x`, the function returns `NA_character_`.
#' @export
#' @details
#' Extract a specific column from the output of GBIF
#' @examples
#'
#' \dontrun{
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
#' }

Extract_BB <- function(x, var) {

  if (is.null(x) || is.null(var)) {
    stop("x, and var cannot be NULL")
  }

  var <- as.character(rlang::ensyms(var))

  if (var %in% names(x)) {
    return(x[[var]])
  } else {
    return(NA_character_)
  }
}
