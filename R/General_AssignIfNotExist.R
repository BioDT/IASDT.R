## |------------------------------------------------------------------------| #
# AssignIfNotExist ----
## |------------------------------------------------------------------------| #
#
#' Assign a value to a variable, only if not existing in the global environment
#'
#' Assign a value to a variable, only if not existing in the global environment
#' @param Variable Variable name
#' @param Value Value to be assigned
#' @param Env environment to assign value to
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @examples
#' AssignIfNotExist(x, TRUE)
#' print(x)
#'
#' # --------------------------------------------------
#'
#' y <- 10
#'
#' # y exists and thus its value was not changed
#' AssignIfNotExist(y, TRUE)
#' print(y)

AssignIfNotExist <- function(Variable, Value, Env = globalenv()) {

  Variable <- as.character(rlang::ensyms(Variable))

  if (!Variable %in% ls(envir = Env)) {
    assign(x = Variable, value = Value, envir = Env)
  } else {
    "The `{Variable}` object already exists in the environment. Current value is:" %>%
      stringr::str_glue() %>%
      crayon::blue() %>%
      cat()

    print(rlang::env_get(Env, paste0(Variable)))
  }
}
