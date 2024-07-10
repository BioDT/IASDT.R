## |------------------------------------------------------------------------| #
# CatTime ----
## |------------------------------------------------------------------------| #

#' Print text with time stamp
#'
#' Print text with time stamp
#'
#' @param Text character; the text to print: default empty string (print time only)
#' @param NLines number of empty lines after the printing; default: 1
#' @param Date Also print date? Default value `FALSE`
#' @param TZ time zone (default: CET)
#' @param ... other arguments passed to `cat`
#' @name CatTime
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @examples
#' CatTime()
#' CatTime(Date = TRUE)
#'
#' CatTime("Time now")
#' CatTime("Time now", Date = TRUE)

CatTime <- function(Text = "", NLines = 1, Date = FALSE, TZ = "CET", ...) {
  DateFormat <- dplyr::if_else(Date, "%d/%m/%Y %X", "%X")
  Now <- lubridate::now(tzone = TZ)
  if (Text == "") {
    cat(format(Now, DateFormat), ...)
    cat(rep("\n", NLines))
  } else {
    Text <- rlang::quo_name(rlang::enquo(Text))
    cat(paste0(Text, " - ", format(Now, DateFormat)), ...)
    cat(rep("\n", NLines))
  }
}
