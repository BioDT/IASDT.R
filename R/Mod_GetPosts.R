## |------------------------------------------------------------------------| #
# GetPosts ----
## |------------------------------------------------------------------------| #

#' Combining posteriors exported by `Hmsc-HPC` into Hmsc object
#'
#' Combining posteriors exported by `Hmsc-HPC` into Hmsc object
#'
#' @param FilePath String. Vector of the file paths of the `rds` files created by `Hmsc-HPC`
#' @param FromJSON Logical. Convert loaded models to JSON format before reading
#' @name GetPosts
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

GetPosts <- function(FilePath, FromJSON = FALSE) {

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "FilePath", Type = "character")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "FromJSON", Type = "logical")
  rm(AllArgs)

  if (FromJSON) {
    readRDS(file = FilePath) %>%
      magrittr::extract2(1) %>%
      jsonify::from_json() %>%
      magrittr::extract2(1) %>%
      return()
  } else {
    readRDS(file = FilePath) %>%
      magrittr::extract2(1) %>%
      # This needs to be checked again
      # jsonify::from_json() %>%
      # magrittr::extract2(1) %>%
      return()
  }
}
