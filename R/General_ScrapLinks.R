## |------------------------------------------------------------------------| #
# ScrapLinks ----
## |------------------------------------------------------------------------| #
#
#' Extracts link texts and URLs from a web page
#'
#' This function scrapes a web page for all links (`<a>` tags) and extracts both
#' the URLs and the link text.
#' @param URL Character. The URL of the web page to scrape. This URL is also
#'   used to resolve relative links to absolute URLs.
#' @name ScrapLinks
#' @return A tibble with two columns: `link_text` containing the text of each
#'   link, and `URL` containing  the absolute URL of each link. The tibble is
#'   sorted by URL and then by link text, and only unique links are included.
#' @importFrom rlang .data
#' @source
#' The source code of this function was taken from this
#'   [gist](https://gist.github.com/paulrougieux/e1ee769577b40cd9ed9db7f75e9a2cc2).
#' @examples
#' ScrapLinks("https://github.com/")
#' @export

ScrapLinks <- function(URL) {

  if (is.null(URL)) {
    stop("URL cannot be NULL", call. = FALSE)
  }

  if (isFALSE(IASDT.R::CheckURL(URL))) {
    stop("Invalid URL", call. = FALSE)
  }

  # Create an html document from the URL
  webpage <- xml2::read_html(URL) %>%
    rvest::html_nodes("a")

  # Extract the URLs
  url_ <- webpage %>%
    rvest::html_attr("href") %>%
    stringr::str_c(URL, ., sep = "") %>%
    stringr::str_replace_all(paste0(URL, URL), URL) %>%
    stringr::str_replace_all("//", "/")

  # Extract the link text
  link_ <- webpage %>%
    rvest::html_text() %>%
    stringr::str_replace_all("\n", "") %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_trim()

  tibble::tibble(link_text = link_, URL = url_) %>%
    dplyr::arrange(.data$URL, .data$link_text) %>%
    dplyr::distinct() %>%
    return()
}
