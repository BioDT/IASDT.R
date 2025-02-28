## |------------------------------------------------------------------------| #
# LoadMultiple ----
## |------------------------------------------------------------------------| #

#' Load multiple `.RData` files together
#'
#' This function loads multiple `.RData` files either into a single list object
#' or directly into the global environment. It provides options for verbosity,
#' returning object names, and handling of non-existent files.
#' @param Files Character vector. The paths of `.RData` files to be loaded.
#' @param Verbose Logical. Whether to print progress messages. Default: `TRUE`.
#' @param OneObject Logical. Whether to load all objects into a single list
#'   (`TRUE`) or directly into the global environment (`FALSE`). Defaults to
#'   `TRUE`.
#' @param ReturnNames Logical. Whether to return the names of the loaded
#'   objects. Only effective when `OneObject` is `FALSE`. Defaults to `TRUE`.
#' @author Ahmed El-Gabbas
#' @return If `OneObject` is `TRUE`, returns a named list of all objects loaded
#'   from the specified files. If `OneObject` is `FALSE` and `ReturnNames` is
#'   `TRUE`, returns a character vector of the names of the objects loaded into
#'   the global environment. Otherwise, returns `NULL`.
#' @name LoadMultiple
#' @examples
#' (Files <- system.file(
#'    "testdata", c("culcita_dat.RData", "gopherdat2.RData"), package = "lme4"))
#'
#' ls()
#'
#' # ---------------------------------------------------
#' # Load multiple *.RData files to one list object
#' # `OneObject = TRUE`
#' # ---------------------------------------------------
#' MultiObj <- LoadMultiple(Files = Files, OneObject = TRUE)
#' ls()
#'
#' str(MultiObj, 1)
#'
#' # ---------------------------------------------------
#' # Load multiple *.RData files separately to the current environment
#' # `OneObject = FALSE`
#' # ---------------------------------------------------
#' ls()
#' rm("MultiObj")
#' ls()
#'
#' LoadMultiple(Files = Files, OneObject = FALSE)
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
#' try(LoadMultiple(Files = Files, OneObject = FALSE))
#'
#' @export

LoadMultiple <- function(
    Files = NULL, Verbose = TRUE, OneObject = TRUE, ReturnNames = TRUE) {

  if (!inherits(Files, "character")) {
    stop("Files should be character object", call. = FALSE)
  }

  if (!all(file.exists(Files))) {
    stop("Some of these files do not exist. No objects were loaded!",
         call. = FALSE)
  }


  if (OneObject) {

    ListNames <- stringr::str_replace_all(basename(Files), "\\.RData$", "")

    Out <- purrr::map(.x = Files, .f = IASDT.R::LoadAs) %>%
      stats::setNames(ListNames)
    return(Out)

  } else {

    Env1 <- new.env()

    purrr::walk(.x = Files, .f = ~ load(file = .x, envir = Env1))

    if (any(ls(envir = Env1) %in% ls(envir = parent.frame()))) {
      stop(
        "Some of the new object names already exists in the current ",
        "environment. No files were loaded!", call. = FALSE)
    }

    ObjectsLoaded <- purrr::map(
      .x = Files, .f = ~load(file = .x, envir = .GlobalEnv)) %>%
      unlist()

    if (Verbose) {
      cat(
        paste0(crayon::blue(
          "Object:", crayon::red(ObjectsLoaded), "was loaded successfully.")),
        sep = "\n")
    }
    if (ReturnNames) {
      return(ObjectsLoaded)
    } else {
      return(invisible(NULL))
    }
  }
}
