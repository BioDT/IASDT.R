## |------------------------------------------------------------------------| #
# chelsa_variables ----
## |------------------------------------------------------------------------| #

#' Detailed information on CHELSA climate variables
#'
#' @return A tibble containing detailed information about various climate
#'   variables in `CHELSA` climate data. Each row in the tibble represents a
#'   different climate variable, with columns providing additional details
#'   such as the long name, unit, scale, offset, and an explanation of the
#'   variable
#'
#' @examples
#' ecokit::load_packages(dplyr, gt)
#'
#' data(chelsa_variables)
#' gt::gt(chelsa_variables)

"chelsa_variables"
