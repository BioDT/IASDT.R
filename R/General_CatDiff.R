## |------------------------------------------------------------------------| #
# CatDiff ----
## |------------------------------------------------------------------------| #
#
#' Print time difference
#'
#' This function calculates the time difference from a given initial time to the
#' current time and prints it with a specified prefix. Optionally, it can also
#' print a session summary.
#'
#' @name CatDiff
#' @author Ahmed El-Gabbas
#' @param InitTime `POSIXct`. The initial time from which the difference is
#'   calculated.
#' @param ChunkText Character. The message printed as chunk info. Default
#'   value: `Session summary`. See: [InfoChunk] for more information.
#' @param Prefix Character. Prefix to prepend to the printed time difference.
#'   Defaults to "Completed in ".
#' @param CatInfo Logical. If `TRUE`, prints a session summary using
#'   [IASDT.R::InfoChunk] ("Session summary"). Defaults to `FALSE`.
#' @param ... Additional arguments for [CatTime].
#' @return The function is used for its side effect of printing to the console
#'   and does not return any value.
#' @inheritParams CatTime
#' @export
#' @examples
#' RefTime <- (lubridate::now() - lubridate::seconds(45))
#' CatDiff(RefTime)
#'
#' RefTime <- (lubridate::now() -
#'     (lubridate::minutes(50) + lubridate::seconds(45)))
#' CatDiff(RefTime)
#'
#' RefTime <- (lubridate::now() - lubridate::minutes(50))
#' CatDiff(RefTime)
#'
#' RefTime <- (lubridate::now() - lubridate::minutes(70))
#' CatDiff(RefTime)
#'
#' RefTime <- (lubridate::now() - lubridate::hours(4))
#' CatDiff(RefTime)
#'
#' RefTime <- lubridate::now() -
#'   (lubridate::hours(4) + lubridate::minutes(50) + lubridate::seconds(45))
#'   CatDiff(RefTime)

CatDiff <- function(
    InitTime, ChunkText = "Session summary", Prefix = "Completed in ",
    CatInfo = FALSE, Level = 0, Time = FALSE, ...) {

  if (is.null(InitTime)) {
    stop("InitTime cannot be NULL", call. = FALSE)
  }

  if (CatInfo) {
    IASDT.R::InfoChunk(Message = ChunkText)
    Prefix <- paste0("\n", Prefix)
  }

  lubridate::time_length(lubridate::now(tzone = "CET") - InitTime) %>%
    lubridate::seconds_to_period() %>%
    {
      paste0(
        stringr::str_pad(
          (lubridate::hour(.) + 24 * lubridate::day(.)), width = 2, pad = "0"),
        ":", stringr::str_pad(lubridate::minute(.), width = 2, pad = "0"), ":",
        stringr::str_pad(round(lubridate::second(.)), width = 2, pad = "0"))
    } %>%
    paste0(Prefix, .) %>%
    IASDT.R::CatTime(Level = Level, Time = Time, ...)

  return(invisible(NULL))
}
