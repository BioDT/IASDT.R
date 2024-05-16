# |---------------------------------------------------| #
# SaveSessionInfo ----
# |---------------------------------------------------| #
#
#' Save all objects (except functions) of the global environment as list items
#'
#' Save all objects (except functions) of the global environment as list items
#'
#' @param Path Path of where to save the output RData file
#' @param SessionObj List of objects and their sizes (typically a a result of `IASDT::SaveSession` function)
#' @param Prefix file prefix
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

SaveSessionInfo <- function(Path = getwd(), SessionObj = NULL, Prefix = "S") {

  FileName <- lubridate::now(tzone = "CET") %>%
    purrr::map_chr(
      .f = ~{
        c(lubridate::year(.x), lubridate::month(.x),
          lubridate::day(.x), "__",
          lubridate::hour(.x), lubridate::minute(.x)) %>%
          sapply(stringr::str_pad, width = 2, pad = "0") %>%
          stringr::str_c(collapse = "") %>%
          stringr::str_replace_all("__", "_") %>%
          stringr::str_c(Prefix, "_", ., collapse = "_")
      }) %>%
    stringr::str_c(Path, "/", ., ".txt")

  IASDT.R::InfoChunk("Session Info") %>%
    utils::capture.output(file = FileName, append = TRUE)
  sessioninfo::session_info() %>%
    utils::capture.output(file = FileName, append = TRUE)

  if (magrittr::not(is.null(SessionObj))) {
    utils::capture.output(
      IASDT.R::InfoChunk("Objects in the current session (except functions and\npre-selected objects; Size in megabytes)"),
      file = FileName, append = TRUE)

    sink(FileName, append = TRUE)
    print.data.frame(tibble::tibble(SessionObj), row.names = FALSE)
    sink()
  }

  return(invisible(NULL))
}
