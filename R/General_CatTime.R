## |------------------------------------------------------------------------| #
# CatTime ----
## |------------------------------------------------------------------------| #

#' Print text with time stamp
#'
#' This function prints a given text followed by the current time (and
#' optionally the date) to the console. It allows for customization of the time
#' zone, the inclusion of the date, and the number of newline characters to
#' print after the message.
#' @param Text character; the text to print before the timestamp. If empty
#'   (default), only the timestamp is printed.
#' @param NLines integer; the number of newline characters to print after the
#'   message. Default is 1.
#' @param Time logical; whether to include the time in the timestamp. Default is
#'   `TRUE`. If `FALSE`, only the text is printed.
#' @param Date logical; whether to include the date in the timestamp. Only
#'   effective if `Time` is `TRUE`. Default is `FALSE`, meaning only the time is
#'   printed. If `TRUE`, the date is printed in the format "%d/%m/%Y %X".
#' @param TZ character; the time zone to use for the timestamp. Default is
#'   `CET`.
#' @param Level integer; the level at which the message will be printed. If e.g.
#'   Level = 1, the following string will be printed at the beginning of the
#'   message: "   >>>   ". Default is `0`.
#' @param Bold logical; whether to print the text in bold. Default is `FALSE`.
#' @param Red logical; whether to print the text in red. Default is `FALSE`.
#' @param ... additional arguments passed to `cat`.
#' @name CatTime
#' @author Ahmed El-Gabbas
#' @return The function is called for its side effect of printing to the
#'   console.
#' @export
#' @examples
#' CatTime()
#'
#' CatTime(Date = TRUE)
#'
#' CatTime("Time now")
#'
#' CatTime("Time now", Date = TRUE)
#'
#' # The use of levels
#' {
#'   CatTime("Task 1")
#'   CatTime("subtask L1", Level = 1)
#'   CatTime("subtask L2", Level = 2)
#'   CatTime("subtask L3", Level = 3)
#' }

CatTime <- function(
    Text = "", NLines = 1L, Time = TRUE, Bold = FALSE, Red = FALSE,
    Date = FALSE, TZ = "CET", Level = 0L, ...) {

  if (Time) {
    DateFormat <- dplyr::if_else(Date, "%d/%m/%Y %X", "%X")
    Now <- lubridate::now(tzone = TZ) %>%
      format(DateFormat)
    Now2 <- paste(" - ", Now)
  } else {
    Now <- Now2 <- ""
  }

  if (Text == "") {
    cat(Now, ...)
    cat(rep("\n", NLines))
  } else {
    if (Bold) {
      Text <- crayon::bold(Text)
    }
    if (Red) {
      Text <- crayon::red(Text)
    }
    if (Level > 0) {
      Prefix <- rep("  >>>", each = Level) %>%
        paste0(collapse = "") %>%
        paste0("  ")
      Text <- paste0(Prefix, Text)
    }

    cat(paste0(Text, Now2), ...)
    cat(rep("\n", NLines))
  }
}
