## |------------------------------------------------------------------------| #
# ScrapLinks ----
## |------------------------------------------------------------------------| #
#
#' Extract link texts and urls from a web page
#'
#' Extract link texts and urls from a web page
#'
#' @param url the url
#' @name ScrapLinks
#' @return NULL
#' @importFrom rlang .data
#' @references https://gist.github.com/paulrougieux/e1ee769577b40cd9ed9db7f75e9a2cc2
#' @examples
#' ScrapLinks("https://github.com/")
#' @export

ScrapLinks <- function(url) {
  # Create an html document from the url
  webpage <- xml2::read_html(url) %>%
    rvest::html_nodes("a")

  # Extract the URLs
  url_ <- webpage %>%
    rvest::html_attr("href") %>%
    stringr::str_c(url, ., sep = "") %>%
    stringr::str_replace_all(paste0(url, url), url) %>%
    stringr::str_replace_all("//", "/")

  # Extract the link text
  link_ <- webpage %>%
    rvest::html_text() %>%
    stringr::str_replace_all("\n", "") %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_trim()

  tibble::tibble(link_text = link_, url = url_) %>%
    dplyr::arrange(.data$url, .data$link_text) %>%
    dplyr::distinct() %>%
    return()
}
