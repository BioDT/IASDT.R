## |------------------------------------------------------------------------| #
# CheckArgs ----
## |------------------------------------------------------------------------| #

#' Check function arguments for specific types
#'
#' This function checks if the specified arguments of a function match the
#' expected type. It is useful for validating function inputs.
#'
#' @name CheckArgs
#' @author Ahmed El-Gabbas
#' @param AllArgs Character vector. Input parameters of the function.
#'   Usually as a result of `formals()` function
#' @param Args Character vector. Names of the arguments to be checked.
#' @param Type Character. The expected type of the arguments. Must be
#'   one of "character", "logical", or "numeric".
#' @return The function does not return a value but will stop execution and
#'   throw an error if any of the specified arguments do not match the expected
#'   type.
#' @export

CheckArgs <- function(AllArgs, Args, Type) {

  if (is.null(AllArgs) || is.null(Args) || is.null(Type)) {
    stop("AllArgs, Args, or Type cannot be NULL", call. = FALSE)
  }

  Type <- match.arg(arg = Type, choices = c("character", "logical", "numeric"))

  if (Type == "character") {
    MissingArgs <- AllArgs[Args] %>%
      purrr::map(~inherits(.x, "character") && all(nzchar(.x))) %>%
      purrr::keep(.p = Negate(isTRUE)) %>%
      names() %>%
      sort()

    if (length(MissingArgs) > 0) {
      stop(
        "The following character argument(s) must be provided\n  >>  ",
        paste(MissingArgs, collapse = " | "), call. = FALSE)
    }
  }

  if (Type == "logical") {
    MissingArgs <- AllArgs[Args] %>%
      purrr::map(.f = inherits, what = "logical") %>%
      purrr::keep(.p = Negate(isTRUE)) %>%
      names() %>%
      sort()

    if (length(MissingArgs) > 0) {
      stop(
        "The following argument(s) must be logical\n  >>  ",
        paste(MissingArgs, collapse = " | "), call. = FALSE)
    }
  }

  if (Type == "numeric") {
    MissingArgs <- AllArgs[Args] %>%
      purrr::map(~(inherits(.x, "numeric") || inherits(.x, "integer"))) %>%
      purrr::keep(.p = Negate(isTRUE)) %>%
      names() %>%
      sort()

    if (length(MissingArgs) > 0) {
      stop(
        "The following argument(s) must be numeric or integer\n  >>  ",
        paste(MissingArgs, collapse = " | "), call. = FALSE)
    }
  }

  return(invisible(NULL))
}
