# |---------------------------------------------------| #
# GetBashVars ----
# |---------------------------------------------------| #

#' Loading arguments passed to Rscript bash command
#'
#' Loading arguments passed to Rscript bash command
#' @name GetBashVars
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

GetBashVars <- function() {
  commandArgs(trailingOnly = TRUE) %>%
    stringr::str_split("=", simplify = TRUE) %>%
    as.data.frame() %>%
    stats::setNames(c("Var", "Value")) %>%
    tibble::tibble() %>%
    dplyr::mutate(Eval = purrr::walk2(
      Var, Value, assign, envir = globalenv())) %>%
    invisible()
}
