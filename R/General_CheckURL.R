## |------------------------------------------------------------------------| #
# CheckURL ----
## |------------------------------------------------------------------------| #
#
#' Check the validity of a URL
#'
#' This function opens a connection to the specified URL to check its validity.
#' It returns `TRUE` if the URL is valid (i.e., the connection can be opened),
#' and `FALSE` otherwise.
#' @param URL Character. The URL to be checked.
#' @param t Numeric. Timeout in seconds for the connection attempt. Default is
#'   2 seconds.
#' @name CheckURL
#' @source The source code of this function was taken from this
#'   [stackoverflow](https://stackoverflow.com/q/52911812) discussion.
#' @return A logical value: `TRUE` if the URL is valid, `FALSE` if not.
#' @examples
#' urls <- c(
#'      "http://www.amazon.com", "http://this.isafakelink.biz",
#'      "https://stackoverflow.com", "https://stack-overflow.com")
#' sapply(urls, CheckURL)
#' @export

CheckURL <- function(URL, t = 2) {

  if (is.null(URL)) {
    stop("URL cannot be NULL", call. = FALSE)
  }

  con <- url(URL)
  check <- suppressWarnings(
    try(open.connection(con, open = "rt", timeout = t), silent = TRUE)[1])

  suppressWarnings(try(close.connection(con), silent = TRUE))

  return(is.null(check))
}
