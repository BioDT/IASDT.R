## |------------------------------------------------------------------------| #
# GetBashVars ----
## |------------------------------------------------------------------------| #

#' Loading arguments passed to Rscript bash command
#'
#' Loading arguments passed to Rscript bash command
#'
#' @name GetBashVars
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

GetBashVars <- function() {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Var <- Value <- NULL

  commandArgs(trailingOnly = TRUE) %>%
    stringr::str_split("=", simplify = TRUE) %>%
    as.data.frame() %>%
    stats::setNames(c("Var", "Value")) %>%
    tibble::tibble() %>%
    dplyr::mutate(Eval = purrr::walk2(
      Var, Value, assign, envir = globalenv())) %>%
    invisible()
}
