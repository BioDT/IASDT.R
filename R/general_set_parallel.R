## |------------------------------------------------------------------------| #
# set_parallel ----
## |------------------------------------------------------------------------| #

#' Set up or stop parallel processing plan
#'
#' Configures parallel processing with [future::plan()] or stops an existing
#' plan. When stopping, it resets to sequential mode.
#'
#' @param n_cores Integer. Number of cores to use. If `NULL`, defaults to
#'   sequential mode. Default is `1`.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are `future::sequential` (sequential), `future::multisession`
#'   (default), `future::multicore` (not supported on Windows), and
#'   `future::cluster`. If `strategy` is not one of the valid options or if
#'   `future::multicore` on Windows PC, it defaults to `future::multisession`.
#'   See [future::plan()] for more details.
#' @param stop Logical. If `TRUE`, stops any parallel cluster and resets to
#'   sequential mode. If `FALSE` (default), sets up a new plan.
#' @param show_log Logical. If `TRUE` (default), logs messages via
#'   [IASDT.R::cat_time()].
#' @param future_max_size Numeric. Maximum allowed total size (in megabytes) of
#'   global variables identified. See `future.globals.maxSize` argument of
#'   [future::future.options] for more details. Default is `8`.
#' @param ... Additional arguments to pass to [cat_time].
#' @export
#' @name set_parallel
#' @author Ahmed El-Gabbas
#' @examples
#' # Prepare working in parallel
#' IASDT.R::set_parallel(n_cores = 2)
#' future::plan()
#'
#' # ---------------------------------------------
#'
#' # Stopping parallel processing
#' IASDT.R::set_parallel(stop = TRUE)
#' future::plan()
#'
#' # ---------------------------------------------
#'
#' # Prepare working in parallel using `future::cluster`
#' IASDT.R::set_parallel(n_cores = 2, strategy = "future::cluster")
#' future::plan()
#'
#' # Stopping parallel processing
#' IASDT.R::set_parallel(stop = TRUE)
#' future::plan()

set_parallel <- function(
    n_cores = 1L, strategy = "future::multisession", stop = FALSE,
    show_log = TRUE, future_max_size = 8L, ...) {

  # Validate n_cores input
  n_cores <- ifelse((is.null(n_cores) || n_cores < 1), 1L, as.integer(n_cores))

  # n_cores can not be more than the available cores
  AvailableCores <- parallelly::availableCores()
  n_cores <- ifelse(
    n_cores > AvailableCores,
    {
      warning(
        "`n_cores` > number of available cores. ",
        "It was reset to the number of available cores: ", AvailableCores,
        call. = FALSE)
      AvailableCores
    },
    n_cores)

  if (stop) {
    if (show_log) {
      IASDT.R::cat_time("Stopping parallel processing", ...)
    }

    # stop any running future plan and reset to sequential
    future::plan("future::sequential", gc = TRUE)

  } else {

    # strategy can not be NULL
    strategy <- ifelse(
      is.null(strategy),
      {
        message(
          "`strategy` cannot be NULL. It was reset to `future::multisession`",
          call. = FALSE)
        "future::multisession"
      },
      strategy)

    # strategy should be a character vector of length 1
    if (length(strategy) != 1) {
      stop("`strategy` must be a character vector of length 1", call. = FALSE)
    }

    # strategy can be only one of the following: "future::sequential",
    # "future::multisession", "future::multicore", "future::cluster".
    ValidStrategy <- c(
      "future::sequential", "future::multisession", "future::multicore",
      "future::cluster")
    strategy <- ifelse(
      (strategy %in% ValidStrategy),
      strategy,
      {
        warning(
          "`strategy` must be one of the following: `",
          paste(ValidStrategy, collapse = "`, `"),
          "` It was reset to `future::multisession`", call. = FALSE)
        "future::multisession"
      })

    # "future::multicore" can not be used on Windows.
    if (strategy == "future::multicore" && .Platform$OS.type == "windows") {
      warning(
        "`future::multicore` is not supported on Windows. ",
        "It was reset to `future::multisession`", call. = FALSE)
      strategy <- "future::multisession"
    }


    if (show_log) {
      IASDT.R::cat_time(
        paste(
          "Setting up", ifelse(n_cores > 1, "parallel", "sequential"),
          "processing with", n_cores, "core(s)"), ...)
    }

    withr::local_options(
      future.globals.maxSize = future_max_size * 1024^2,
      future.gc = TRUE, future.seed = TRUE,
      .local_envir = parent.frame())

    if (n_cores > 1) {
      future::plan(strategy = strategy, workers = n_cores)
    } else {
      future::plan("future::sequential", gc = TRUE)
    }
  }
  return(invisible(NULL))
}
