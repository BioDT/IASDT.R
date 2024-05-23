## |------------------------------------------------------------------------| #
# Valid_URL ----
## |------------------------------------------------------------------------| #
#
#' Check the validity of a URL
#'
#' Check the validity of a URL
#' @param url_in URL path
#' @param t Timeout
#' @name Valid_URL
#' @references https://stackoverflow.com/questions/52911812/check-if-url-exists-in-r
#' @return NULL
#' @examples
#' urls <- c("http://www.amazon.com", "http://this.isafakelink.biz", "https://stackoverflow.com")
#' sapply(urls, Valid_URL)
#' @export

Valid_URL <- function(url_in, t = 2) {
  con <- url(url_in)
  check <- suppressWarnings(try(open.connection(con, open = "rt", timeout = t), silent = TRUE)[1])
  suppressWarnings(try(close.connection(con), silent = TRUE))
  ifelse(is.null(check), TRUE, FALSE)
}
