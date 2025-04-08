## |------------------------------------------------------------------------| #
# bash_variables ----
## |------------------------------------------------------------------------| #

#' Read command line arguments passed to an R script
#'
#' This function reads command line arguments passed to an R script executed via
#' the bash command line and assigns them into the global environment of the R
#' session. This allows for dynamic passing of variables from a bash script to
#' an R script.
#' @name bash_variables
#' @author Ahmed El-Gabbas
#' @return This function does not return anything but has the side effect of
#'   assigning variables into the global environment.
#' @export

bash_variables <- function() {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Var <- Value <- NULL

  CommandArgs <- commandArgs(trailingOnly = TRUE)

  if (length(CommandArgs) > 0) {
    # Convert command line arguments into a data frame
    CommandArgs %>%
      stringr::str_split("=", simplify = TRUE) %>%
      as.data.frame(stringsAsFactors = FALSE) %>%
      stats::setNames(c("Var", "Value")) %>%
      tibble::tibble() %>%
      dplyr::mutate(
        Eval = purrr::walk2(Var, Value, assign, envir = globalenv()))
  }
  return(invisible(NULL))
}
