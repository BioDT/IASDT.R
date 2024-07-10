## |------------------------------------------------------------------------| #
# CatDiff ----
## |------------------------------------------------------------------------| #
#
#' Print time difference
#'
#' Print time difference
#'
#' @name CatDiff
#' @author Ahmed El-Gabbas
#' @return NULL
#' @param InitTime day and time of the time event
#' @param Prefix Prefix character to print
#' @param CatInfo logical; also print `IASDT.R::InfoChunk("Session summary")`
#' @export

CatDiff <- function(InitTime, Prefix = "Completed in ", CatInfo = TRUE) {

  if (CatInfo) {
    IASDT.R::InfoChunk("Session summary")
  }

  (lubridate::now(tzone = "CET") - InitTime) %>%
    lubridate::time_length(unit = "min") %>%
    round(2) %>%
    paste0(Prefix, ., " minutes") %>%
    IASDT.R::CatTime(... = "\n")

  return(invisible(NULL))
}
