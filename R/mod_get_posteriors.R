## |------------------------------------------------------------------------| #
# mod_get_posteriors ----
## |------------------------------------------------------------------------| #

#' Combines posteriors exported by `Hmsc-HPC` into an Hmsc object
#'
#' This function converts posterior files exported by `Hmsc-HPC` into an Hmsc
#' object. It can either read the data directly from RDS files or convert it
#' from JSON format if specified.
#' @param path_posterior Character vector. Path to the RDS files containing the
#'   exported posterior files. This argument is mandatory and cannot be empty.
#' @param from_JSON Logical. Whether the loaded models should be converted from
#'   `JSON` format. Defaults to `FALSE`, meaning the data will be read directly
#'   from RDS files without conversion.
#' @name mod_get_posteriors
#' @author Ahmed El-Gabbas
#' @return Depending on the `from_JSON` parameter, returns an Hmsc object either
#'   directly from the RDS files or after converting it from JSON format.
#' @export

mod_get_posteriors <- function(path_posterior = NULL, from_JSON = FALSE) {

  # Check if path_posterior is empty
  if (is.null(path_posterior)) {
    IASDT.R::stop_ctx(
      "`path_posterior` cannot be empty", path_posterior = path_posterior)
  }

  # Checking arguments
  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::check_args(
    args_all = AllArgs, args_to_check = "path_posterior",
    args_type = "character")
  IASDT.R::check_args(
    args_all = AllArgs, args_to_check = "from_JSON", args_type = "logical")

  if (from_JSON) {
    Out <- readRDS(file = path_posterior) %>%
      magrittr::extract2(1) %>%
      jsonify::from_json() %>%
      magrittr::extract2(1)
    return(Out)
  } else {
    Out <- readRDS(file = path_posterior) %>%
      magrittr::extract2(1) %>%
      magrittr::extract2(1)
    return(Out)
  }
}
