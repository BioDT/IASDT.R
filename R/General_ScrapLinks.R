## |------------------------------------------------------------------------| #
# ScrapLinks ----
## |------------------------------------------------------------------------| #

#' Extracts link texts and URLs from a web page
#'
#' This function scrapes a web page for all links (`<a>` tags) and extracts both
#' the URLs and the link text.
#' @param URL Character. The URL of the web page to scrape. This URL is also
#'   used to resolve relative links to absolute URLs.
#' @param Arrange Character vector of length 1 or 2. The columns to arrange the
#'   output by. The default is c("Link", "Link_text"). The first column is the
#'   URL of the link, and the second column is the text of the link. The
#'   function will arrange the output in ascending order by the column(s)
#'   specified in this argument.
#' @name ScrapLinks
#' @return A tibble with two columns: `Link_text` containing the text of each
#'   link, and `URL` containing  the absolute URL of each link. The tibble is
#'   sorted by URL and then by link text, and only unique links are included.
#' @importFrom rlang .data
#' @examples
#'

#' head(
#' ScrapLinks(URL = "https://github.com/BioDT/IASDT.R"))
#'
#' head(
#'   ScrapLinks(
#'     URL = "https://github.com/BioDT/IASDT.R", Arrange = "Link_text"))
#'
#' # This will give an "Invalid URL" error
#' \dontrun{
#'  ScrapLinks(URL = "https://github50.com")
#' }
#' @export

ScrapLinks <- function(URL, Arrange = c("Link", "Link_text")) {

  Link <- Link_text <- NULL

  # Ensure that Arrange is a character vector of length 1 or 2
  if (!is.character(Arrange) || length(Arrange) > 2 ||
      length(Arrange) < 1) {
    stop(
      "`Arrange` must be a character vector of length 1 or 2",
      call. = FALSE)
  }

  # Ensure that all values of Arrange are in c("Link", "Link_text")
  if (!all(Arrange %in% c("Link", "Link_text"))) {
    stop(
      "`Arrange` must contain only 'Link' and 'Link_text'", call. = FALSE)
  }

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
  Out <- purrr::map(
    .x = webpage,
    .f = ~ tibble::tibble(
      Link = stringr::str_trim(rvest::html_attr(.x, "href")),
      Link_text = stringr::str_trim(rvest::html_text(.x)))) %>%
    dplyr::bind_rows() %>%
    # Remove empty or anchor links
    dplyr::filter(
      !is.na(Link) & !stringr::str_starts(Link, "#"),
      Link != "..", nzchar(Link_text)) %>%
    dplyr::mutate(
      Link = dplyr::if_else(
        stringr::str_starts(Link, "http"), Link,
        IASDT.R::Path(
          stringr::str_remove(URL, "/$"), stringr::str_remove(Link, "^/")
        )),

      Link_text = {
        stringr::str_remove_all(Link_text, "\n") %>%
          stringr::str_replace_all("\\s+", " ") %>%
          stringr::str_trim()
      }) %>%
    dplyr::distinct() %>%
    dplyr::select(Link_text, Link) %>%
    dplyr::arrange(dplyr::across(tidyselect::all_of(Arrange)))

  return(Out)
}
