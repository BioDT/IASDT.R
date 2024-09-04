## |------------------------------------------------------------------------| #
# CatFunArgs ----
## |------------------------------------------------------------------------| #

#' Print function Arguments
#'
#' This function takes another function as input and prints its arguments in the
#' format `ArgumentName = DefaultValue`. Each argument is printed on a new line.
#' The function can optionally assign the formatted arguments to the global
#' environment and can load a specified package before processing.
#' @param Function A function whose arguments you want to print. Must be a valid
#'   R function.
#' @param assign Logical. Whether to assign the arguments as variables in the
#'   global environment. Defaults to `FALSE`.
#' @param Package Character string. The name of the R package to be loaded
#'   before processing the function. Default is NULL.
#' @export
#' @name CatFunArgs
#' @author Ahmed El-Gabbas
#' @return The function prints the formatted arguments to the console. If
#'   `Assign` is TRUE, it also assigns arguments to the global environment.
#' @examples
#' formals(stats::setNames)
#' CatFunArgs(stats::setNames)

CatFunArgs <- function(Function, Assign = FALSE, Package = NULL) {

  if (!is.function(Function)) {
    stop("The provided function is not a function", call. = FALSE)
  }

  if (!is.null(Package)) {
    IASDT.R::LoadPackages(List = Package)
  }

  # Extract the formal arguments of the function
  args_list <- formals(Function)

  formatted_args <- purrr::map_chr(
    .x = seq_along(args_list),
    .f = function(i) {
      arg_name <- names(args_list)[i]

      if (is.name(args_list[[i]])) {
        arg_value <- NULL
      } else {
        arg_value <- args_list[[i]]
      }

      PostProcess <- !any(
        c(
          is.numeric(arg_value),
          is.logical(arg_value),
          is.null(arg_value)
        )
      )

      if (PostProcess) {
        arg_value <- deparse(args_list[[i]]) %>%
          stringr::str_replace_all('^\"|\"$', "") %>%
          stringr::str_replace_all('\"', '"')
      }

      Output <- if (is.null(arg_value)) {
        paste0(arg_name, " = NULL")
      } else if (is.character(arg_value)) {
        if (stringr::str_detect(arg_value, "\\(|\\)")) {
          paste0(arg_name, " = ", arg_value)
        } else {
          paste0(arg_name, " = \"", arg_value, "\"")
        }
      } else {
        paste0(arg_name, " = ", as.character(arg_value))
      }

      if (Assign == TRUE) {
        eval(expr = parse(text = Output), envir = .GlobalEnv)
      }
      return(Output)
    }
  )

  cat(formatted_args, sep = "\n")
  return(invisible(NULL))
}
