## |------------------------------------------------------------------------| #
# Valid_URL ----
## |------------------------------------------------------------------------| #
#
#' Check the validity of a URL
#'
#' This function opens a connection to the specified URL to check its validity.
#' It returns `TRUE` if the URL is valid (i.e., the connection can be opened),
#' and `FALSE` otherwise.
#' @param url_in A character string specifying the URL to be checked.
#' @param t A numeric value specifying the timeout in seconds for the connection
#'   attempt. Default is 2 seconds.
#' @name Valid_URL
#' @source The source code of this function was taken from this
#'   [stackoverflow](https://stackoverflow.com/q/52911812) discussion.
#' @return A logical value: `TRUE` if the URL is valid, `FALSE` if not.
#' @examples
#' urls <- c(
#'      "http://www.amazon.com", "http://this.isafakelink.biz",
#'      "https://stackoverflow.com", "https://stack-overflow.com")
#' sapply(urls, Valid_URL)
#' @export

Valid_URL <- function(url_in, t = 2) {

  if (is.null(url_in)) {
    stop("url_in cannot be NULL")
  }

  con <- url(url_in)
  check <- suppressWarnings(
    try(open.connection(con, open = "rt", timeout = t), silent = TRUE)[1])

  suppressWarnings(try(close.connection(con), silent = TRUE))

  return(ifelse(is.null(check), TRUE, FALSE))
}
