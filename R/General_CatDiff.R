## |------------------------------------------------------------------------| #
# CatDiff ----
## |------------------------------------------------------------------------| #
#
#' Print time difference
#'
#' This function calculates the time difference from a given initial time to the current time and prints it with a specified prefix. Optionally, it can also print a session summary.
#'
#' @name CatDiff
#' @author Ahmed El-Gabbas
#' @return NULL
#' @param InitTime POSIXct; The initial time from which the difference is calculated.
#' @param Prefix character; A prefix string to prepend to the printed time difference. Defaults to "Completed in ".
#' @param CatInfo logical; If `TRUE`, prints a session summary using [IASDT.R::InfoChunk] ("Session summary"). Defaults to `TRUE`.
#' @return NULL; The function is used for its side effect of printing to the console and does not return any value.
#' @examples
#' # Assuming the current time, it prints the time difference from one minute ago.
#' CatDiff(Sys.time() - 60)
#' @export

CatDiff <- function(InitTime, Prefix = "Completed in ", CatInfo = TRUE) {
  
  if (is.null(InitTime)) {
    stop("InitTime cannot be NULL")
  }

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
