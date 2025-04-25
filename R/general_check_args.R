## |------------------------------------------------------------------------| #
# check_args ----
## |------------------------------------------------------------------------| #

#' Check function arguments for specific types
#'
#' This function checks if the specified arguments of a function match the
#' expected type. It is useful for validating function inputs.
#'
#' @name check_args
#' @author Ahmed El-Gabbas
#' @param args_all Character vector. Input parameters of the function.
#'   Usually as a result of `formals()` function
#' @param args_to_check Character vector. Names of the arguments to be
#'   checked.
#' @param args_type Character. The expected type of the arguments. Must be
#'   one of "character", "logical", or "numeric".
#' @return The function does not return a value but will stop execution and
#'   throw an error if any of the specified arguments do not match the expected
#'   type.
#' @export

check_args <- function(args_all, args_to_check, args_type) {

  if (is.null(args_all) || is.null(args_to_check) || is.null(args_type)) {
    IASDT.R::stop_ctx(
      "`args_all`, `args_to_check`, or `args_type` cannot be NULL",
      args_all = args_all, args_to_check = args_to_check,
      args_type = args_type)
  }

  args_type <- match.arg(
    arg = args_type, choices = c("character", "logical", "numeric"))

  if (args_type == "character") {
    MissingArgs <- args_all[args_to_check] %>%
      purrr::map(~inherits(.x, "character") && all(nzchar(.x))) %>%
      purrr::keep(.p = Negate(isTRUE)) %>%
      names() %>%
      sort()

    if (length(MissingArgs) > 0) {
      IASDT.R::stop_ctx(
        paste0(
          "The following character argument(s) must be provided\n  >>  ",
          paste(MissingArgs, collapse = " | ")),
        length_MissingArgs = length(MissingArgs))
    }
  }

  if (args_type == "logical") {
    MissingArgs <- args_all[args_to_check] %>%
      purrr::map(.f = inherits, what = "logical") %>%
      purrr::keep(.p = Negate(isTRUE)) %>%
      names() %>%
      sort()

    if (length(MissingArgs) > 0) {
      IASDT.R::stop_ctx(
        paste0(
          "The following argument(s) must be logical\n  >>  ",
          paste(MissingArgs, collapse = " | ")),
        length_MissingArgs = length(MissingArgs))
    }
  }

  if (args_type == "numeric") {
    MissingArgs <- args_all[args_to_check] %>%
      purrr::map(~(inherits(.x, "numeric") || inherits(.x, "integer"))) %>%
      purrr::keep(.p = Negate(isTRUE)) %>%
      names() %>%
      sort()

    if (length(MissingArgs) > 0) {
      IASDT.R::stop_ctx(
        paste0(
          "The following argument(s) must be numeric or integer\n  >>  ",
          paste(MissingArgs, collapse = " | ")),
        length_MissingArgs = length(MissingArgs))
    }
  }

  return(invisible(NULL))
}
