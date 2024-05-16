# |---------------------------------------------------| #
# git_log ----
# |---------------------------------------------------| #

#' print detailed `git log` of the git repo located in the current working directory
#'
#' print detailed `git log` of the git repo located in the current working directory
#'
#' @name git_log
#' @export
#' @examples
#' git_log()

git_log <- function() {
  'git log --graph --pretty=format:"%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr) %C(bold blue)<%an>%Creset" --abbrev-commit' %>%
    IASDT.R::System() %>%
    cat(sep = "\n")
}
