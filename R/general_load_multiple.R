## |------------------------------------------------------------------------| #
# load_multiple ----
## |------------------------------------------------------------------------| #

#' Load multiple `.RData` files together
#'
#' This function loads multiple `.RData` files either into a single list object
#' or directly into the global environment. It provides options for verbosity,
#' returning object names, and handling of non-existent files.
#' @param files Character vector. The paths of `.RData` files to be loaded.
#' @param verbose Logical. Whether to print progress messages. Default: `TRUE`.
#' @param single_object Logical. Whether to load all objects into a single list
#'   (`TRUE`) or directly into the global environment (`FALSE`). Defaults to
#'   `TRUE`.
#' @param return_names Logical. Whether to return the names of the loaded
#'   objects. Only effective when `single_object` is `FALSE`. Defaults to
#'   `TRUE`.
#' @author Ahmed El-Gabbas
#' @return If `single_object` is `TRUE`, returns a named list of all objects
#'   loaded from the specified files. If `single_object` is `FALSE` and
#'   `return_names` is `TRUE`, returns a character vector of the names of the
#'   objects loaded into the global environment. Otherwise, returns `NULL`.
#' @name load_multiple
#' @examples
#' (files <- system.file(
#'    "testdata", c("culcita_dat.RData", "gopherdat2.RData"), package = "lme4"))
#'
#' ls()
#'
#' # ---------------------------------------------------
#' # Load multiple *.RData files to one list object
#' # `single_object = TRUE`
#' # ---------------------------------------------------
#' MultiObj <- load_multiple(files = files, single_object = TRUE)
#' ls()
#'
#' str(MultiObj, 1)
#'
#' # ---------------------------------------------------
#' # Load multiple *.RData files separately to the current environment
#' # `single_object = FALSE`
#' # ---------------------------------------------------
#' ls()
#' rm("MultiObj")
#' ls()
#'
#' load_multiple(files = files, single_object = FALSE)
#'
#' str(Gdat, 1)
#'
#' str(culcita_dat, 1)
#'
#' # ---------------------------------------------------
#' # Load multiple *.RData files, one object already exists
#' # ---------------------------------------------------
#' ls()
#' rm("culcita_dat")
#' ls()
#'
#' try(load_multiple(files = files, single_object = FALSE))
#'
#' @export

load_multiple <- function(
    files = NULL, verbose = TRUE, single_object = TRUE, return_names = TRUE) {

  if (!inherits(files, "character")) {
    stop("`files` should be character object", call. = FALSE)
  }

  if (!all(file.exists(files))) {
    stop(
      "Some of these files do not exist. No objects were loaded!",
      call. = FALSE)
  }


  if (single_object) {

    ListNames <- stringr::str_replace_all(basename(files), "\\.RData$", "")

    Out <- purrr::map(.x = files, .f = IASDT.R::load_as) %>%
      stats::setNames(ListNames)
    return(Out)

  } else {

    Env1 <- new.env()

    purrr::walk(.x = files, .f = ~ load(file = .x, envir = Env1))

    if (any(ls(envir = Env1) %in% ls(envir = parent.frame()))) {
      stop(
        "Some of the new object names already exists in the current ",
        "environment. No files were loaded!", call. = FALSE)
    }

    ObjectsLoaded <- purrr::map(
      .x = files, .f = ~load(file = .x, envir = .GlobalEnv)) %>%
      unlist()

    if (verbose) {
      cat(
        paste0(crayon::blue(
          "Object:", crayon::red(ObjectsLoaded), "was loaded successfully.")),
        sep = "\n")
    }
    if (return_names) {
      return(ObjectsLoaded)
    } else {
      return(invisible(NULL))
    }
  }
}
