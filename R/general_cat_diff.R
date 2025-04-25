## |------------------------------------------------------------------------| #
# cat_diff ----
## |------------------------------------------------------------------------| #
#
#' Print time difference
#'
#' This function calculates the time difference from a given initial time to the
#' current time and prints it with a specified prefix. Optionally, it can also
#' print a session summary.
#'
#' @name cat_diff
#' @author Ahmed El-Gabbas
#' @param init_time `POSIXct`. The initial time from which the difference is
#'   calculated.
#' @param chunk_text Character. The message printed as chunk info. Default
#'   value: `Session summary`. See: [info_chunk] for more information.
#' @param prefix Character. prefix to prepend to the printed time difference.
#'   Defaults to "Completed in ".
#' @param cat_info Logical. If `TRUE`, prints a session summary using
#'   [IASDT.R::info_chunk] ("Session summary"). Defaults to `FALSE`.
#' @param ... Additional arguments for [cat_time].
#' @return The function is used for its side effect of printing to the console
#'   and does not return any value.
#' @inheritParams cat_time
#' @export
#' @examples
#' RefTime <- (lubridate::now() - lubridate::seconds(45))
#' cat_diff(RefTime)
#'
#' RefTime <- (lubridate::now() -
#'     (lubridate::minutes(50) + lubridate::seconds(45)))
#' cat_diff(RefTime)
#'
#' RefTime <- (lubridate::now() - lubridate::minutes(50))
#' cat_diff(RefTime)
#'
#' RefTime <- (lubridate::now() - lubridate::minutes(70))
#' cat_diff(RefTime)
#'
#' RefTime <- (lubridate::now() - lubridate::hours(4))
#' cat_diff(RefTime)
#'
#' RefTime <- lubridate::now() -
#'   (lubridate::hours(4) + lubridate::minutes(50) + lubridate::seconds(45))
#' cat_diff(RefTime)

cat_diff <- function(
    init_time, chunk_text = "Session summary", prefix = "Completed in ",
    cat_info = FALSE, level = 0L, cat_timestamp = FALSE, ...) {

  if (is.null(init_time)) {
    IASDT.R::stop_ctx("`init_time` cannot be NULL", init_time = init_time)
  }

  if (cat_info) {
    IASDT.R::info_chunk(message = chunk_text)
    prefix <- paste0("\n", prefix)
  }

  Period <- lubridate::time_length(
    lubridate::now(tzone = "CET") - init_time) %>%
    lubridate::seconds_to_period()
  Period_hours <- stringr::str_pad(
    (lubridate::hour(Period) + 24 * lubridate::day(Period)),
    width = 2, pad = "0")
  Period_minutes <- stringr::str_pad(
    lubridate::minute(Period), width = 2, pad = "0")
  Period_seconds <- stringr::str_pad(
    round(lubridate::second(Period)), width = 2, pad = "0")

  paste0(Period_hours, ":", Period_minutes, ":", Period_seconds) %>%
    paste0(prefix, .) %>%
    IASDT.R::cat_time(level = level, cat_timestamp = cat_timestamp, ...)

  return(invisible(NULL))
}
