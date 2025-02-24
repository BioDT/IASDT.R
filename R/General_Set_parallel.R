## |------------------------------------------------------------------------| #
# Set_parallel ----
## |------------------------------------------------------------------------| #

#' Set up or stop parallel processing plan
#'
#' Configures parallel processing with `future::plan()` or stops an existing
#' plan. When stopping, it resets to sequential mode.
#'
#' @param NCores Integer. Number of cores to use. If `NULL`, defaults to
#'   sequential mode.
#' @param Stop Logical. If `TRUE`, stops any parallel cluster and resets to
#'   sequential mode. If `FALSE` (default), sets up a new plan.
#' @param Cat Logical. If `TRUE` (default), logs messages via
#'   [IASDT.R::CatTime()].
#' @param Level Integer. The logging level for [CatTime]. Default is `0`.
#' @param Set_parallel Numeric. Maximum allowed total size (in megabytes) of
#'   global variables identified. See `` argument of [future::future.options]
#'   for more details
#' @return Invisible `NULL`. On Windows with `NCores > 1`, the cluster object
#'   (`c1`) is assigned to the parent environment and cleaned up automatically
#'   on exit.
#' @export
#' @name Set_parallel
#' @author Ahmed El-Gabbas
#' @examples
#' \dontrun{
#'   # Prepare working on parallel
#'   IASDT.R::Set_parallel(NCores = min(NCores, nrow(Beta_DF)), Level = 3)
#'
#'   if (.Platform$OS.type == "windows") {
#'     on.exit(try(snow::stopCluster("c1"), silent = TRUE), add = TRUE)
#'   }
#'   on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
#'
#'   # Stopping cluster
#'   IASDT.R::Set_parallel(Stop = TRUE, Cat = TRUE, Level = 3)
#' }

Set_parallel <- function(
    NCores = 1L, Stop = FALSE, Cat = TRUE, Level = 0L, Future_maxSize = 8) {

  # Validate NCores input
  NCores <- ifelse(
    (is.null(NCores) || NCores < 1),
    1L, as.integer(NCores))

  # NCores can not be more than the available cores
  NCores <- ifelse(
    NCores > parallelly::availableCores(),
    {
      warning(
        paste0(
          "`NCores` > number of available cores. ",
          "It was reset to the number of available cores",
          parallelly::availableCores()),
        call. = FALSE)
      parallelly::availableCores()
    },
    NCores)

  if (Stop) {
    if (Cat) {
      IASDT.R::CatTime("Stopping parallel processing", Level = Level)
    }

    # Stop any running future plan and reset to sequential
    future::plan("future::sequential", gc = TRUE)

  } else {

    if (Cat) {
      IASDT.R::CatTime(
        paste(
          "Setting up", ifelse(NCores > 1, "parallel", "sequential"),
          "processing with", NCores, "core(s)"),
        Level = Level)
    }

    withr::local_options(
      future.globals.maxSize = Future_maxSize * 1024^2,
      future.gc = TRUE, future.seed = TRUE,
      .local_envir = parent.frame())

    if (NCores > 1) {
      future::plan("future::multisession", workers = NCores)
    } else {
      future::plan("future::sequential", gc = TRUE)
    }
  }
  return(invisible(NULL))
}
