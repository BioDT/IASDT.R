## |------------------------------------------------------------------------| #
# AssignIfNotExist ----
## |------------------------------------------------------------------------| #
#
#' Assign a value to a variable if it does not already exist in the specified environment
#'
#' This function checks if a given variable exists in the specified environment (global environment by default). If the variable does not exist, it assigns a given value to it. If the variable already exists, it prints the current value of the variable. The function is designed to prevent overwriting existing variables unintentionally.
#'
#' @param Variable Character; the name of the variable to be checked and potentially assigned a value.
#' @param Value any; the value to be assigned to the variable if it does not already exist.
#' @param Env environment; the environment in which to check for the existence of the variable and potentially assign the value. Defaults to the global environment.
#' @author Ahmed El-Gabbas
#' @return NULL; The function explicitly returns NULL, but its primary effect is the side-effect of assigning a value to a variable in an environment or printing the current value of an existing variable.
#' @export
#' @examples
#' AssignIfNotExist(Variable = "x", Value = TRUE)
#' print(x)
#'
#' # --------------------------------------------------
#'
#' y <- 10
#'
#' # y exists and thus its value was not changed
#' AssignIfNotExist(Variable = "y", Value = TRUE)
#' print(y)

AssignIfNotExist <- function(Variable, Value, Env = globalenv()) {

  if (is.null(Variable) || is.null(Value)) {
    stop("Variable and Value cannot be NULL")
  }

  Variable <- as.character(rlang::ensyms(Variable))

  if (!exists(Variable, envir = Env)) {
    assign(x = Variable, value = Value, envir = Env)
  } else {
    "The `{Variable}` object already exists in the environment. Current value is:" %>%
      stringr::str_glue() %>%
      crayon::blue() %>%
      cat()

    print(rlang::env_get(Env, paste0(Variable)))
  }
}
