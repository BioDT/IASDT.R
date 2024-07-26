## |------------------------------------------------------------------------| #
# git_log ----
## |------------------------------------------------------------------------| #

#' Print a detailed `git log` of the git repository located in the current working directory.
#'
#' This function executes a `git log` command with specific formatting options to display a detailed and colored log. It shows the commit hash, references (like branches or tags), commit message, relative commit date, and author name in a visually appealing graph format.
#'
#' @return This function is called for its side effect of printing to the console.
#' @name git_log
#' @export
#' @examples
#' git_log()

git_log <- function() {
  'git log --graph --pretty=format:"%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr) %C(bold blue)<%an>%Creset" --abbrev-commit' %>%
    system(intern = TRUE) %>%
    cat(sep = "\n")
}
