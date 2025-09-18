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
#' @param from_json Logical. Whether the loaded models should be converted from
#'   `JSON` format. Defaults to `FALSE`, meaning the data will be read directly
#'   from RDS files without conversion.
#' @name mod_get_posteriors
#' @author Ahmed El-Gabbas
#' @return Depending on the `from_json` parameter, returns an Hmsc object either
#'   directly from the RDS files or after converting it from JSON format.
#' @export

mod_get_posteriors <- function(path_posterior = NULL, from_json = FALSE) {

  # Checking arguments
  ecokit::check_args(
    args_to_check = "path_posterior", args_type = "character")
  ecokit::check_args(args_to_check = "from_json", args_type = "logical")

  if (!ecokit::check_data(path_posterior)) {
    ecokit::stop_ctx(
      "path_posterior does not exist or invalid",
      path_posterior = path_posterior)
  }

  if (from_json) {
    out <- readRDS(file = path_posterior) %>%
      magrittr::extract2(1) %>%
      jsonify::from_json() %>%
      magrittr::extract2(1)
  } else {
    out <- readRDS(file = path_posterior) %>%
      magrittr::extract2(1) %>%
      magrittr::extract2(1)
  }
  return(out)
}
