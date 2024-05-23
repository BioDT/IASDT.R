## |------------------------------------------------------------------------| #
# sort_ ----
## |------------------------------------------------------------------------| #

#' Sort alphanumeric strings
#'
#' Sort alphanumeric strings
#' @param x Vector to be sorted.
#' @param decreasing logical. Should the sort be increasing or decreasing? Note that descending=TRUE reverses the meanings of na.last and blanks.last. Default: `FALSE`
#' @param na.last for controlling the treatment of NA values. If TRUE, missing values in the data are put last; if FALSE, they are put first; if NA, they are removed. Default: `TRUE`
#' @param blank.last for controlling the treatment of blank values. If TRUE, blank values in the data are put last; if FALSE, they are put first; if NA, they are removed. Default: `FALSE`
#' @param numeric.type either "decimal" (default) or "roman". Are numeric values represented as decimal numbers (numeric.type="decimal") or as Roman numerals (numeric.type="roman")?
#' @param roman.case one of "upper", "lower", or "both". Are roman numerals represented using only capital letters ('IX') or lower-case letters ('ix') or both?
#' @name sort_
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' # example code
#' (AA <- paste0("V", 1:12))
#'
#' sort(AA)
#'
#' sort_(AA)
#'
#' @export

sort_ <- function(
    x, decreasing = FALSE, na.last = TRUE,
    blank.last = FALSE, numeric.type = c("decimal", "roman"),
    roman.case = c("upper", "lower", "both")) {

  gtools::mixedsort(
    x = x, decreasing = decreasing, na.last = na.last,
    blank.last = blank.last,
    numeric.type = numeric.type,
    roman.case = roman.case)
}
