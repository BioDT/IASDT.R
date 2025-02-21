## |------------------------------------------------------------------------| #
# git_log ----
## |------------------------------------------------------------------------| #

#' Print or return a detailed `git log` of the git repository located in the
#' specified directory.
#'
#' This function checks if the specified directory is a Git repository and, if
#' so, executes a `git log` command to either print the log to the console or
#' return it. It supports execution on Windows and Linux operating systems and
#' provides a visually appealing graph format of the log, showing the commit
#' hash, references, commit message, relative commit date, and author name.
#' @param Path Character. Path to the directory to check. Defaults to the
#'   current working directory ".". If the path does not exist, the function
#'   will stop and throw an error. If the path is not a git repository, the
#'   function will throw a warning.
#' @param Num Integer. Number of recent commits to display. If `NULL` (the
#'   default), the complete log is shown. If `Num` is not `NULL` or a positive
#'   number, the function will stop and throw an error. This parameter is
#'   ignored if `ReturnLog` is `TRUE`.
#' @param ReturnLog Logical. Whether to return the log (`TRUE`) or print it to
#'   the console (`FALSE`, default). If `TRUE`, the function returns a character
#'   vector containing the log lines.
#' @return If `ReturnLog` is `TRUE`, returns a character vector containing the
#'   git log lines. If `ReturnLog` is `FALSE`, the function is called for its
#'   side effect of printing to the console.
#' @note The function will stop and throw an error if the specified path does
#'   not exist, the operating system is not supported, the directory is not a
#'   Git repository, Git is not installed, or if the `Num` parameter is not
#'   `NULL` or a positive number.
#' @name git_log
#' @author Ahmed El-Gabbas
#' @export
#' @seealso [IASDT.R::System]
#' @examples
#' # not a git repo
#' git_log(Path = "C:/")
#'
#' # Show the most recent commit
#' git_log(Num = 1)
#'
#' # Show the most recent 5 commits
#' git_log(Num = 5)
#'
#' # Return the log as a character vector
#' Log <- git_log(ReturnLog = TRUE)
#'
#' length(Log)
#'
#' head(Log, 8)

git_log <- function(Path = ".", Num = NULL, ReturnLog = FALSE) {

  # Check system commands
  Commands <- c("git")
  CommandsAvail <- purrr::map_lgl(Commands, IASDT.R::CheckCommands)
  if (!all(CommandsAvail)) {
    Missing <- paste0(Commands[!CommandsAvail], collapse = " + ")
    stop(
      paste0("The following command(s) are not available: ", Missing),
      call. = FALSE)
  }

  if (!dir.exists(Path)) {
    stop("The provided path does not exist.", call. = FALSE)
  }

  # Determine the OS-specific command
  os <- IASDT.R::CurrOS()

  if (!os %in% c("Windows", "Linux")) {
    stop(
      "Unsupported OS. This function supports only Windows and Linux.",
      call. = FALSE)
  }

  # Construct the command to check if the directory is a Git repo
  git_check_command <- if (os == "Windows") {
    paste0(
      'cmd.exe /c "cd /d ', IASDT.R::NormalizePath(Path),
      ' && git rev-parse --is-inside-work-tree"')
  } else {
    paste0(
      'sh -c "cd ', IASDT.R::NormalizePath(Path),
      ' && git rev-parse --is-inside-work-tree"')
  }

  # Check if the directory is a Git repo
  is_git <- tryCatch({
    result <- system(
      command = git_check_command, intern = TRUE, ignore.stderr = TRUE) %>%
      suppressWarnings()

    if (length(result) == 0) {
      message("The git_check_command returned an empty result.")
      return(invisible(NULL))
    }
    result
  },
  error = function(e) {
    message("Error checking if directory is a Git repository: ", e$message)
    return(invisible(NULL))
  })

  if (is.null(is_git)) {
    stop("Failed to determine if the directory is a Git repository.")
  }

  if (is_git == "true") {
    # Construct the command to get the Git log
    log_command <- paste0(
      "git -C ", IASDT.R::NormalizePath(Path),
      ' log --graph --pretty=format:"%Cred%h%Creset ',
      "-%C(yellow)%d%Creset %s %Cgreen(%cr) ",
      '%C(bold blue)<%an>%Creset" --abbrev-commit')

    # Execute the command and capture the output
    log_output <- tryCatch({
      IASDT.R::System(log_command, RObj = TRUE)
    },
    error = function(e) {
      stop(
        paste0(
          "Failed to retrieve Git log. Ensure Git is installed and the ",
          "directory is a valid Git repository."),
        call. = FALSE)
    })

    if (isFALSE(ReturnLog)) {
      # Print the log output
      if (is.null(Num)) {
        cat(log_output, sep = "\n")
      } else {
        if (is.numeric(Num) && Num > 0) {
          cat(utils::head(log_output, n = Num), sep = "\n")
        } else {
          stop(paste0(
            "The 'Num' argument can be either NULL to show the ",
            "complete log or a positive numeric value to show the ",
            "most recent commits."), call. = FALSE)
        }
      }
    }
  } else {
    warning("The provided directory is not a Git repository.", call. = FALSE)
  }

  if (ReturnLog) {
    return(log_output)
  } else {
    return(invisible(NULL))
  }
}
