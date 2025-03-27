## |------------------------------------------------------------------------| #
# Set_parallel ----
## |------------------------------------------------------------------------| #

#' Set up or stop parallel processing plan
#'
#' Configures parallel processing with `future::plan()` or stops an existing
#' plan. When stopping, it resets to sequential mode.
#'
#' @param NCores Integer. Number of cores to use. If `NULL`, defaults to
#'   sequential mode. Default is `1`.
#' @param Strategy Character. The parallel processing strategy to use. Valid
#'   options are `future::sequential` (sequential), `future::multisession`
#'   (default), `future::multicore` (not supported on Windows), and
#'   `future::cluster`. If `Strategy` is not one of the valid options or if
#'   `future::multicore` on Windows PC, it defaults to `future::multisession`.
#'   See `future::plan()` for more details.
#' @param Stop Logical. If `TRUE`, stops any parallel cluster and resets to
#'   sequential mode. If `FALSE` (default), sets up a new plan.
#' @param Cat Logical. If `TRUE` (default), logs messages via
#'   [IASDT.R::CatTime()].
#' @param Future_maxSize Numeric. Maximum allowed total size (in megabytes) of
#'   global variables identified. See `future.globals.maxSize` argument of
#'   [future::future.options] for more details. Default is `8`.
#' @param ... Additional arguments to pass to [CatTime].
#' @export
#' @name Set_parallel
#' @author Ahmed El-Gabbas
#' @examples
#' # Prepare working on parallel
#' IASDT.R::Set_parallel(NCores = 2)
#' future::plan()
#'
#' # ---------------------------------------------
#'
#' # Stopping parallel processing
#' IASDT.R::Set_parallel(Stop = TRUE)
#' future::plan()
#'
#' # ---------------------------------------------
#'
#' # Prepare working on parallel using `future::cluster`
#' IASDT.R::Set_parallel(NCores = 2, Strategy = "future::cluster")
#' future::plan()
#'
#' # Stopping parallel processing
#' IASDT.R::Set_parallel(Stop = TRUE)
#' future::plan()

Set_parallel <- function(
    NCores = 1L, Strategy = "future::multisession", Stop = FALSE, Cat = TRUE,
    Future_maxSize = 8L, ...) {

  # Validate NCores input
  NCores <- ifelse((is.null(NCores) || NCores < 1), 1L, as.integer(NCores))

  # NCores can not be more than the available cores
  AvailableCores <- parallelly::availableCores()
  NCores <- ifelse(
    NCores > AvailableCores,
    {
      warning(
        "`NCores` > number of available cores. ",
        "It was reset to the number of available cores: ", AvailableCores,
        call. = FALSE)
      AvailableCores
    },
    NCores)

  if (Stop) {
    if (Cat) {
      IASDT.R::CatTime("Stopping parallel processing", ...)
    }

    # Stop any running future plan and reset to sequential
    future::plan("future::sequential", gc = TRUE)

  } else {

    # Strategy can not be NULL
    Strategy <- ifelse(
      is.null(Strategy),
      {
        message(
          "`Strategy` cannot be NULL. It was reset to `future::multisession`",
          call. = FALSE)
        "future::multisession"
      },
      Strategy)

    # Strategy should be a character vector of length 1
    if (length(Strategy) != 1) {
      stop("`Strategy` must be a character vector of length 1", call. = FALSE)
    }

    # Strategy can be only one of the following: "future::sequential",
    # "future::multisession", "future::multicore", "future::cluster".
    ValidStrategy <- c(
      "future::sequential", "future::multisession", "future::multicore",
      "future::cluster")
    Strategy <- ifelse(
      (Strategy %in% ValidStrategy),
      Strategy,
      {
        warning(
          "`Strategy` must be one of the following: `",
          paste(ValidStrategy, collapse = "`, `"),
          "` It was reset to `future::multisession`", call. = FALSE)
        "future::multisession"
      })

    # "future::multicore" can not be used on Windows.
    if (Strategy == "future::multicore" && .Platform$OS.type == "windows") {
      warning(
        "`future::multicore` is not supported on Windows. ",
        "It was reset to `future::multisession`", call. = FALSE)
      Strategy <- "future::multisession"
    }


    if (Cat) {
      IASDT.R::CatTime(
        paste(
          "Setting up", ifelse(NCores > 1, "parallel", "sequential"),
          "processing with", NCores, "core(s)"), ...)
    }

    withr::local_options(
      future.globals.maxSize = Future_maxSize * 1024^2,
      future.gc = TRUE, future.seed = TRUE,
      .local_envir = parent.frame())

    if (NCores > 1) {
      future::plan(strategy = Strategy, workers = NCores)
    } else {
      future::plan("future::sequential", gc = TRUE)
    }
  }
  return(invisible(NULL))
}
