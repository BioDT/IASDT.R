## |------------------------------------------------------------------------| #
# cat_time ----
## |------------------------------------------------------------------------| #

#' Print text with time stamp
#'
#' This function prints a given text followed by the current time (and
#' optionally the date) to the console. It allows for customization of the time
#' zone, the inclusion of the date, and the number of newline characters to
#' print after the message.
#' @param text character; the text to print before the timestamp. If empty
#'   (default), only the timestamp is printed.
#' @param n_lines integer; the number of newline characters to print after the
#'   message. Default is 1.
#' @param time logical; whether to include the time in the timestamp. Default is
#'   `TRUE`. If `FALSE`, only the text is printed.
#' @param date logical; whether to include the date in the timestamp. Only
#'   effective if `time` is `TRUE`. Default is `FALSE`, meaning only the time is
#'   printed. If `TRUE`, the date is printed in the format "%d/%m/%Y %X".
#' @param time_zone character; the time zone to use for the timestamp. Default
#'   is `CET`.
#' @param level integer; the level at which the message will be printed. If e.g.
#'   level = 1, the following string will be printed at the beginning of the
#'   message: "   >>>   ". Default is `0`.
#' @param bold logical; whether to print the text in bold. Default is `FALSE`.
#' @param red logical; whether to print the text in red. Default is `FALSE`.
#' @param ... additional arguments passed to `cat`.
#' @name cat_time
#' @author Ahmed El-Gabbas
#' @return The function is called for its side effect of printing to the
#'   console.
#' @export
#' @examples
#' cat_time()
#'
#' cat_time(date = TRUE)
#'
#' cat_time("time now")
#'
#' cat_time("\n\nTime now", n_lines = 2L, level = 1L)
#'
#' cat_time(
#'   "\ntime now", date = TRUE, bold = TRUE, red = TRUE,
#'   n_lines = 2L, level = 1L)
#'
#' # The use of levels
#' {
#'   cat_time("Task 1")
#'   cat_time("subtask L1", level = 1L)
#'   cat_time("subtask L2", level = 2L)
#'   cat_time("subtask L3", level = 3L)
#' }

cat_time <- function(
    text = "", n_lines = 1L, time = TRUE, bold = FALSE,
    red = FALSE, date = FALSE, time_zone = "CET", level = 0L, ...) {

  # Validate inputs
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "logical",
    args_to_check = c("time", "bold", "red", "date"))
  IASDT.R::check_args(
    args_all = AllArgs, args_type = "numeric",
    args_to_check = c("n_lines", "level"))
  rm(AllArgs, envir = environment())

  # Current time
  Now <- lubridate::now(tzone = time_zone)

  # Format date / time
  if (date && time) {
    Now <- format(Now, "%d/%m/%Y %X")
    Now2 <- paste0(" - ", Now)
  } else if (date) {
    Now <- format(Now, "%d/%m/%Y")
    Now2 <- paste0(" - ", Now)
  } else if (time) {
    Now <- format(Now, "%X")
    Now2 <- paste0(" - ", Now)
  } else {
    Now <- Now2 <- ""
  }

  NLinesBefore <- stringr::str_extract(text, "^\\n+") %>%
    stringr::str_count("\n")
  if (is.na(NLinesBefore)) {
    NLinesBefore <- 0
  }

  text <- stringr::str_remove(text, "^\\n+")

  if (text == "") {
    if (NLinesBefore > 0) {
      text <- paste0(strrep("\n", NLinesBefore), text)
    }
    text <- paste0(text, Now)
  } else {
    if (level > 0) {
      Prefix <- rep("  >>>", each = level) %>%
        paste(collapse = "") %>%
        paste0("  ")
      text <- paste0(Prefix, text)
    }

    if (NLinesBefore > 0) {
      text <- paste0(strrep("\n", NLinesBefore), text)
    }
    text <- paste0(text, Now2)

  }

  if (bold) {
    text <- crayon::bold(text)
  }
  if (red) {
    text <- crayon::red(text)
  }

  cat(text, ...)
  cat(rep("\n", n_lines))

  return(invisible(NULL))
}
