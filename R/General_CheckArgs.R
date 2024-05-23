## |------------------------------------------------------------------------| #
# CheckArgs ----
## |------------------------------------------------------------------------| #

#' Check function arguments
#'
#' Check function arguments
#' @name CheckArgs
#' @author Ahmed El-Gabbas
#' @return NULL
#' @param AllArgs String. Vector for the input parameters of the function. Usually as a result of `formals()` function
#' @param Args Arguments to check
#' @param Type Class of the arguments to check. This has to be one of "character", "logical", or "numeric".
#' @export

CheckArgs <- function(AllArgs, Args, Type) {

  match.arg(arg = Type, choices = c("character", "logical", "numeric"))

  if (Type == "character") {
    MissingArgs <- AllArgs[Args] %>%
      purrr::map(~inherits(.x, "character") && nchar(.x) > 0) %>%
      purrr::discard(.p = isTRUE) %>%
      names() %>%
      sort()

    if (length(MissingArgs) > 0) {
      MSG <- paste0(
        "The following character argument(s) must be provided\n  >>  ",
        paste0(MissingArgs, collapse = " | ")) %>%
        stop(call. = FALSE)
    }
  }

  if (Type == "logical") {
    MissingArgs <- AllArgs[Args] %>%
      purrr::map(~inherits(.x, "logical")) %>%
      purrr::discard(.p = isTRUE) %>%
      names() %>%
      sort()

    if (length(MissingArgs) > 0) {
      MSG <- paste0(
        "The following argument(s) must be logical\n  >>  ",
        paste0(MissingArgs, collapse = " | ")) %>%
        stop(call. = FALSE)
    }
  }

  if (Type == "numeric") {
    MissingArgs <- AllArgs[Args] %>%
      purrr::map(~inherits(.x, "numeric")) %>%
      purrr::discard(.p = isTRUE) %>%
      names() %>%
      sort()

    if (length(MissingArgs) > 0) {
      MSG <- paste0(
        "The following argument(s) must be numeric (integer)\n  >>  ",
        paste0(MissingArgs, collapse = " | ")) %>%
        stop(call. = FALSE)
    }
  }
}
