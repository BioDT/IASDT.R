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
  n_cores
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

  hab_abb
}
