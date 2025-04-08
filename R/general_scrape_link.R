## |------------------------------------------------------------------------| #
# scrape_link ----
## |------------------------------------------------------------------------| #

#' Extracts link texts and URLs from a web page
#'
#' This function scrapes a web page for all links (`<a>` tags) and extracts both
#' the URLs and the link text.
#' @param URL Character. The URL of the web page to scrape. This URL is also
#'   used to resolve relative links to absolute URLs.
#' @param sort_by Character vector of length 1 or 2. The columns to arrange the
#'   output by. The default is c("Link", "Link_text"). The first column is the
#'   URL of the link, and the second column is the text of the link. The
#'   function will arrange the output in ascending order by the column(s)
#'   specified in this argument.
#' @name scrape_link
#' @return A tibble with two columns: `Link_text` containing the text of each
#'   link, and `URL` containing  the absolute URL of each link. The tibble is
#'   sorted by URL and then by link text, and only unique links are included.
#' @importFrom rlang .data
#' @examples
#'
#' head(
#' scrape_link(URL = "https://github.com/BioDT/IASDT.R"))
#'
#' head(
#'   scrape_link(
#'     URL = "https://github.com/BioDT/IASDT.R", sort_by = "Link_text"))
#'
#' # This will give an "Invalid URL" error
#' \dontrun{
#'  scrape_link(URL = "https://github50.com")
#' }
#' @export

scrape_link <- function(URL, sort_by = c("Link", "Link_text")) {

  Link <- Link_text <- NULL

  # Ensure that sort_by is a character vector of length 1 or 2
  if (!is.character(sort_by) || length(sort_by) > 2 ||
      length(sort_by) < 1) {
    stop(
      "`sort_by` must be a character vector of length 1 or 2",
      call. = FALSE)
  }

  # Ensure that all values of sort_by are in c("Link", "Link_text")
  if (!all(sort_by %in% c("Link", "Link_text"))) {
    stop(
      "`sort_by` must contain only 'Link' and 'Link_text'", call. = FALSE)
  }

  if (is.null(URL)) {
    stop("URL cannot be NULL", call. = FALSE)
  }

  if (isFALSE(IASDT.R::check_URL(URL))) {
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
        IASDT.R::path(
          stringr::str_remove(URL, "/$"), stringr::str_remove(Link, "^/")
        )),

      Link_text = {
        stringr::str_remove_all(Link_text, "\n") %>%
          stringr::str_replace_all("\\s+", " ") %>%
          stringr::str_trim()
      }) %>%
    dplyr::distinct() %>%
    dplyr::select(Link_text, Link) %>%
    dplyr::arrange(dplyr::across(tidyselect::all_of(sort_by)))

  return(Out)
}
