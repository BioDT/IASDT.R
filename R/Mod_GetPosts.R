## |------------------------------------------------------------------------| #
# Mod_GetPosts ----
## |------------------------------------------------------------------------| #

#' Combines posteriors exported by `Hmsc-HPC` into an Hmsc object
#'
#' This function converts posterior files exported by `Hmsc-HPC` into an Hmsc
#' object. It can either read the data directly from RDS files or convert it
#' from JSON format if specified.
#' @param Path_Post Character vector. Path to the RDS files containing the
#'   exported posterior files. This argument is mandatory and cannot be empty.
#' @param FromJSON Logical. Whether the loaded models should be converted from
#'   `JSON` format. Defaults to `FALSE`, meaning the data will be read directly
#'   from RDS files without conversion.
#' @name Mod_GetPosts
#' @author Ahmed El-Gabbas
#' @return Depending on the `FromJSON` parameter, returns an Hmsc object either
#'   directly from the RDS files or after converting it from JSON format.
#' @export

Mod_GetPosts <- function(Path_Post = NULL, FromJSON = FALSE) {

  # Check if Path_Post is empty
  if (is.null(Path_Post)) {
    stop("Path_Post cannot be empty", call. = FALSE)
  }

  # Checking arguments
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "Path_Post", Type = "character")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "FromJSON", Type = "logical")

  if (FromJSON) {
    Out <- readRDS(file = Path_Post) %>%
      magrittr::extract2(1) %>%
      jsonify::from_json() %>%
      magrittr::extract2(1)
    return(Out)
  } else {
    Out <- readRDS(file = Path_Post) %>%
      magrittr::extract2(1) %>%
      magrittr::extract2(1)
    return(Out)
  }
}
