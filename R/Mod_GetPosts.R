## |------------------------------------------------------------------------| #
# GetPosts ----
## |------------------------------------------------------------------------| #

#' Combines posteriors exported by `Hmsc-HPC` into an Hmsc object
#'
#' This function converts posteriors exported by `Hmsc-HPC` into an Hmsc object.
#' It can either read the data directly from RDS files or convert it from JSON
#' format if specified.
#' @param FilePath A character string specifying the path to the RDS files
#'   containing the exported posteriors.
#' @param FromJSON A logical flag indicating whether the loaded models should be
#'   converted from JSON format. Defaults to `FALSE`, meaning the data will be
#'   read directly from an RDS file without conversion.
#' @name GetPosts
#' @author Ahmed El-Gabbas
#' @return Depending on the `FromJSON` parameter, returns an Hmsc object either
#'   directly from the RDS file or after converting it from JSON format.
#' @export

GetPosts <- function(FilePath, FromJSON = FALSE) {

  if (is.null(FilePath)) {
    stop("FilePath cannot be empty", call. = FALSE)
  }

  # Checking arguments
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "FilePath", Type = "character")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "FromJSON", Type = "logical")

  if (FromJSON) {
    readRDS(file = FilePath) %>%
      magrittr::extract2(1) %>%
      jsonify::from_json() %>%
      magrittr::extract2(1) %>%
      return()
  } else {
    readRDS(file = FilePath) %>%
      magrittr::extract2(1) %>%
      magrittr::extract2(1) %>%
      return()
  }
}
