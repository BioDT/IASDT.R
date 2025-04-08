## |------------------------------------------------------------------------| #
# assign_env_vars ------
## |------------------------------------------------------------------------| #

#' Assign environment variables from a .env file
#'
#' This function reads environment variables from a specified `.env` file and
#' assigns them to variables in the R environment based on a data frame
#' containing variable names, values, and checks for directories and files. It
#' is designed to facilitate the management of environment variables in a
#' structured and reproducible manner.
#' @name assign_env_vars
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param env_variables_data `data.frame`. A data frame or tibble containing the
#'   columns `VarName`, `Value`, `CheckDir`, and `CheckFile`. Each row specifies
#'   an environment variable, the name to assign it to, and whether to check if
#'   it is a directory or file. This structure allows for additional validation
#'   checks on the variables being imported.
#' @author Ahmed El-Gabbas
#' @export
#' @return This function is used for its side effects of setting environment
#'   variables in the R environment. It assigns each variable from the `.env`
#'   file to the R environment with the name specified in the
#'   `env_variables_data` data frame.

assign_env_vars <- function(env_file = ".env", env_variables_data = NULL) {

  if (is.null(env_file) || is.null(env_variables_data)) {
    stop(
      "`env_file` and `env_variables_data` can not be empty",
      call. = FALSE)
  }

  if (!file.exists((env_file))) {
    stop(env_file, " file does not exist", call. = FALSE)
  }

  if (!inherits(env_variables_data, "data.frame")) {
    stop(
      "The provided env_variables_data object should be either tibble or ",
      "data frame",
      call. = FALSE)
  }

  MatchNames <- setdiff(
    c("VarName", "Value", "CheckDir", "CheckFile"),
    names(env_variables_data))

  if (length(MatchNames) > 0) {
    stop(
      "The following columns are missing from `env_variables_data`: ",
      paste(MatchNames, collapse = "; "), call. = FALSE)
  }

  InClasses <- purrr::map_chr(env_variables_data, class)
  ExpClasses <- c("character", "character", "logical", "logical")
  MatchClass <- all(InClasses == ExpClasses)

  if (isFALSE(MatchClass)) {
    Diff <- which(InClasses != ExpClasses)
    purrr::map_chr(
      .x = Diff,
      .f = ~{
        paste0(
          '"', names(env_variables_data)[.x], '" is ', InClasses[.x],
          " not ", ExpClasses[.x])
      }) %>%
      stringr::str_c(collapse = "\n") %>%
      stringr::str_c("\n", ., collapse = "\n") %>%
      stop(call. = FALSE)
  }

  readRenviron(env_file)

  purrr::walk(
    .x = seq_len(nrow(env_variables_data)),
    .f = ~{
      EV_Name <- env_variables_data$Value[.x]
      Var_Name <- env_variables_data$VarName[.x]
      CheckDir <- env_variables_data$CheckDir[.x]
      CheckFile <- env_variables_data$CheckFile[.x]

      Val <- Sys.getenv(EV_Name)

      if (CheckDir && CheckFile) {
        stop(
          Val, " should be checked as either file or directory, not both",
          call. = FALSE)
      }

      if (!nzchar(Val)) {
        stop(
          "`", EV_Name,
          "` environment variable was not set in the .env file", call. = FALSE)
      }

      if (CheckDir && !dir.exists(Val)) {
        stop("`", Val, "` directory does not exist", call. = FALSE)
      }

      if (CheckFile && !file.exists(Val)) {
        stop("`", Val, "` file does not exist", call. = FALSE)
      }

      assign(x = Var_Name, value = Val, envir = parent.frame(5))
    })

  return(invisible(NULL))
}
