#' @keywords internal
#' @noRd
.validate_n_cores <- function(n_cores) {

  if (is.null(n_cores)) ecokit::stop_ctx("n_cores cannot be NULL")

  if (!is.numeric(n_cores) || length(n_cores) != 1L ||
      n_cores < 1L || is.na(n_cores)) {
    ecokit::stop_ctx(
      "n_cores must be a positive integer of length 1",
      n_cores = n_cores, class_n_cores = class(n_cores))
  }

  n_cores <- as.integer(n_cores)
  max_cores <- parallelly::availableCores()
  if (n_cores > max_cores) {
    warning(
      stringr::str_glue(
        "`n_cores` exceeds available cores: {n_cores}. Using all available",
        " cores: {max_cores}"),
      call. = FALSE)
    n_cores <- max_cores
  }
  as.integer(n_cores)
}

# # ..................................................................... ###

#' @keywords internal
#' @noRd
.validate_strategy <- function(strategy) {

  if (is.null(strategy)) ecokit::stop_ctx("strategy cannot be NULL")

  if (!is.character(strategy) || length(strategy) != 1L) {
    ecokit::stop_ctx(
      "strategy must be a character of length 1",
      strategy = strategy, class_strategy = class(strategy))
  }

  valid_strategies <- c("sequential", "multisession", "multicore", "cluster")
  if (!strategy %in% valid_strategies) {
    ecokit::stop_ctx(
      "Invalid strategy name",
      strategy = strategy, valid_strategies = valid_strategies)
  }

  strategy
}

# # ..................................................................... ###

#' @keywords internal
#' @noRd
.validate_hab_abb <- function(hab_abb) {

  if (is.null(hab_abb)) ecokit::stop_ctx("hab_abb cannot be NULL")

  if (!is.character(hab_abb) || length(hab_abb) != 1L) {
    ecokit::stop_ctx(
      "hab_abb must be a character of length 1",
      hab_abb = hab_abb, class_hab_abb = class(hab_abb))
  }

  if (!nzchar(hab_abb)) {
    ecokit::stop_ctx(
      "hab_abb must be a non-empty character string",
      hab_abb = hab_abb, class_hab_abb = class(hab_abb))
  }

  valid_hab_abbs <- as.character(c(0:3, "4a", "4b", 10, "12a", "12b"))
  if (!hab_abb %in% valid_hab_abbs) {
    ecokit::stop_ctx(
      "Invalid habitat abbreviation",
      hab_abb = hab_abb, valid_hab_abbs = valid_hab_abbs)
  }

  tolower(hab_abb)
}


# # ..................................................................... ###

#' @keywords internal
#' @noRd

.validate_cv_name <- function(cv_name) {

  if (is.null(cv_name)) ecokit::stop_ctx("cv_name cannot be NULL")

  if (!is.character(cv_name) || length(cv_name) != 1L || !nzchar(cv_name)) {
    ecokit::stop_ctx(
      "cv_name must be a character of length 1",
      cv_name = cv_name, class_cv_name = class(cv_name))
  }

  valid_cv_names <- c("CV_Dist", "CV_Large", "CV_SAC")
  if (!cv_name %in% valid_cv_names) {
    ecokit::stop_ctx(
      "Invalid cv_name", cv_name = cv_name, valid_cv_names = valid_cv_names)
  }

  cv_name
}

# # ..................................................................... ###

#' Validate SLURM job runtime string.
#'
#' Checks if the runtime string follows SLURM standards: accepted formats are
#' "minutes", "hours:minutes:seconds", "days-hours", etc.
#'
#' @param runtime Character string, e.g. "2:00:00" or "1-12:00:00".
#' @keywords internal
#' @noRd

# # Valid SLURM runtimes
# .validate_slurm_runtime("30")
# .validate_slurm_runtime("2:00:00")
# .validate_slurm_runtime("0-12:00:00")
# .validate_slurm_runtime("5-23:59")
# .validate_slurm_runtime("2:00")

# # Invalid SLURM runtimes
# .validate_slurm_runtime("2:70:00")
# .validate_slurm_runtime("12-24:00:00")
# .validate_slurm_runtime("abc")

.validate_slurm_runtime <- function(runtime, warning = TRUE) {

  if (is.null(runtime)) {
    ecokit::stop_ctx("runtime cannot be NULL")
  }

  if (length(runtime) != 1L || !is.character(runtime) || !nzchar(runtime)) {
    ecokit::stop_ctx(
      "runtime must be a non-empty character string of length 1",
      runtime = runtime, class_ram = class(runtime),
      length_runtime = length(runtime), include_backtrace = TRUE)
  }

  pattern <- paste0(
    "^(\\d+-)?([0-9]|[0-1][0-9]|2[0-3])",
    ":[0-5][0-9](:[0-5][0-9])?$|^\\d+$")

  if (stringr::str_detect(runtime, pattern, negate = TRUE)) {
    ecokit::stop_ctx(
      "Invalid SLURM runtime format: ",
      runtime = runtime, include_backtrace = TRUE)
  }

  return(runtime)
}

# # ..................................................................... ###

#' Validate SLURM requested RAM string.
#'
#' Checks if the RAM string follows SLURM standards: accepted units are
#' uppercase K, M, G, T, or their long forms (e.g. MB, GB). Only uppercase units
#' are accepted.
#'
#' @param ram Character string, e.g. "8G", "1024M", "32GB".
#' @keywords internal
#' @noRd

# # Example usage for .validate_slurm_ram
# .validate_slurm_ram("8G")
# .validate_slurm_ram("1024M")
# .validate_slurm_ram("32GB")
# .validate_slurm_ram("1T")

# # Invalid SLURM RAM requests
# .validate_slurm_ram("500")
# .validate_slurm_ram("500m")
# .validate_slurm_ram("4.5G")
# .validate_slurm_ram("100MBs")
# .validate_slurm_ram("GB32")

.validate_slurm_ram <- function(ram) {

  if (is.null(ram)) {
    ecokit::stop_ctx("ram cannot be NULL")
  }

  if (length(ram) != 1L || !is.character(ram) || !nzchar(ram)) {
    ecokit::stop_ctx(
      "ram must be a non-empty character string of length 1",
      ram = ram, class_ram = class(ram),
      length_ram = length(ram), include_backtrace = TRUE)
  }

  pattern <- "^\\d+(K|M|G|T|KB|MB|GB|TB)$"
  if (stringr::str_detect(ram, pattern, negate = TRUE)) {
    ecokit::stop_ctx(
      "Invalid SLURM RAM request: ", ram = ram, include_backtrace = TRUE)
  }

  return(ram)
}



# # ..................................................................... ###

#' Reduce a ggplot object to a simplified ggplot grob
#'
#' This function takes a ggplot object and converts it into a simplified ggplot
#' grob using `ggplot2::ggplotGrob` and wraps it as a ggpubr ggplot object.
#'
#' @param gg A ggplot object to be reduced.
#' @return A ggpubr ggplot object containing the simplified grob.
#'
#' @keywords internal
#' @noRd

ggplot_reduce <- function(gg) {
  if (!inherits(gg, "gg")) {
    ecokit::stop_ctx("Input must be a ggplot object")
  }
  ggpubr::as_ggplot(ggplot2::ggplotGrob(gg))
}
