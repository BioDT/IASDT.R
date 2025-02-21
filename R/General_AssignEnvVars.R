## |------------------------------------------------------------------------| #
# AssignEnvVars ------
## |------------------------------------------------------------------------| #

#' Assign environment variables from a .env file
#'
#' This function reads environment variables from a specified `.env` file and
#' assigns them to variables in the R environment based on a data frame
#' containing variable names, values, and checks for directories and files. It
#' is designed to facilitate the management of environment variables in a
#' structured and reproducible manner.
#' @name AssignEnvVars
#' @param EnvFile Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param EnvVarDT `data.frame`. A data frame or tibble containing the columns
#'   `VarName`, `Value`, `CheckDir`, and `CheckFile`. Each row specifies an
#'   environment variable, the name to assign it to, and whether to check if it
#'   is a directory or file. This structure allows for additional validation
#'   checks on the variables being imported.
#' @author Ahmed El-Gabbas
#' @export
#' @return This function is used for its side effects of setting environment
#'   variables in the R environment. It assigns each variable from the `.env`
#'   file to the R environment with the name specified in the `EnvVarDT` data
#'   frame.

AssignEnvVars <- function(EnvFile = ".env", EnvVarDT = NULL) {

  if (is.null(EnvFile) || is.null(EnvVarDT)) {
    stop("EnvFile and EnvVarDT can not be empty", call. = FALSE)
  }

  if (!file.exists((EnvFile))) {
    stop(paste0(EnvFile, " file does not exist"), call. = FALSE)
  }

  if (!inherits(EnvVarDT, "data.frame")) {
    stop("The provided EnvVarDT object should be either tibble or data frame")
  }

  MatchNames <- setdiff(
    c("VarName", "Value", "CheckDir", "CheckFile"),
    names(EnvVarDT))
  if (length(MatchNames) > 0) {
    stop(paste0(
      "The following columns are missing from the EnvVarDT object: ",
      paste0(MatchNames, collapse = "; ")), call. = FALSE)
  }

  InClasses <- purrr::map_chr(EnvVarDT, class)
  ExpClasses <- c("character", "character", "logical", "logical")
  MatchClass <- all(InClasses == ExpClasses)

  if (isFALSE(MatchClass)) {
    Diff <- which(InClasses != ExpClasses)
    purrr::map_chr(
      .x = Diff,
      .f = ~{
        paste0(
          '"', names(EnvVarDT)[.x], '" is ', InClasses[.x],
          " not ", ExpClasses[.x])
      }) %>%
      stringr::str_c(collapse = "\n") %>%
      stringr::str_c("\n", ., collapse = "\n") %>%
      stop(call. = FALSE)
  }

  readRenviron(EnvFile)

  purrr::walk(
    .x = seq_len(nrow(EnvVarDT)),
    .f = ~{
      EV_Name <- EnvVarDT$Value[.x]
      Var_Name <- EnvVarDT$VarName[.x]
      CheckDir <- EnvVarDT$CheckDir[.x]
      CheckFile <- EnvVarDT$CheckFile[.x]

      Val <- Sys.getenv(EV_Name)

      if (CheckDir && CheckFile) {
        stop(
          paste0(
            Val, " should be checked as either file or directory, not both"),
          call. = FALSE)
      }

      if (nchar(Val) == 0) {
        stop(
          paste0(
            "`", EV_Name,
            "` environment variable was not set in the .env file"),
          call. = FALSE)
      }

      if (CheckDir && !dir.exists(Val)) {
        stop(paste0("`", Val, "` directory does not exist"), call. = FALSE)
      }

      if (CheckFile && !file.exists(Val)) {
        stop(paste0("`", Val, "` file does not exist"), call. = FALSE)
      }

      assign(x = Var_Name, value = Val, envir = parent.frame(5))
    })

  return(invisible(NULL))
}
