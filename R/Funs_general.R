## quiets concerns of R CMD check re: the .'s that appear in pipelines
if (getRversion() >= "2.15.1")  utils::globalVariables(c("."))

# |---------------------------------------------------| #
# CatSep ----
# |---------------------------------------------------| #
#
#' Print separator(s) to the console
#'
#' Print separator(s) to the console
#' @param Rep number of separator lines; default 1 row
#' @param Extra1 number of extra empty lines before the separator; default: 0
#' @param Extra2 number of extra empty lines after the separator; default: 0
#' @param Char The character to be used as a separator; default "-"
#' @param CharReps Number of times the character is repeated; default: 50
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' CatSep()
#'
#' CatSep(2)
#'
#' CatSep(2,2,3)
#'
#' CatSep(2,2,3, Char = "*")
#'
#' CatSep(2,2,3, Char = "*", CharReps = 20)
#' @export

CatSep <- function(Rep = 1, Extra1 = 0, Extra2 = 0, Char = "-", CharReps = 50) {
  if (Extra1 > 0) {
    replicate(n = Extra1, expr = cat("\n"))
  }
  S <- c(rep(Char, CharReps)) %>%
    paste0(collapse = "")
  replicate(n = Rep, expr = cat(S, sep = "\n"))
  if (Extra2 > 0) {
    replicate(n = Extra2, expr = cat("\n"))
  }
  return(invisible(NULL))
}


# |---------------------------------------------------| #
# InfoChunk ----
# |---------------------------------------------------| #
#
#' Print Information chunk
#'
#' Print Information chunk
#' @param Message String passed to the `IASDT.R::CatTime` function
#' @param ... additional arguments for the `IASDT.R::CatSep` function
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' InfoChunk(Message = "Started")
#'
#' InfoChunk(Message = "finished", Char = "*", CharReps = 60)
#' @export

InfoChunk <- function(Message = "", ...) {
  IASDT.R::CatSep(..., Extra1 = 1)
  IASDT.R::CatTime(Message)
  IASDT.R::CatSep(..., Extra2 = 1)
}

# |---------------------------------------------------| #
# SaveSession ----
# |---------------------------------------------------| #
#
#' Save all objects (except functions) of the global environment as list items
#'
#' Save all objects (except functions) of the global environment as list items
#'
#' @param Path Path of where to save the output RData file
#' @param ExcludeObs objects not to save
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

SaveSession <- function(Path = getwd(), ExcludeObs = NULL) {

  IASDT.R::DirCreate(Path, Verbose = FALSE)
  ExcludeObs <- c(ExcludeObs, "Grid_10_sf_s", "Grid_10_Raster", "Bound_sf_Eur_s", "Bound_sf_Eur")

  AllObjs <- ls(envir = .GlobalEnv) %>%
    tibble::tibble(Object = .) %>%
    dplyr::mutate(
      Class = purrr::map_chr(
        .x = Object,
        .f = ~{
          get(.x, envir = .GlobalEnv) %>%
            class() %>%
            stringr::str_c(collapse = "_")
        }
      )) %>%
    dplyr::filter(
      Class != "function",
      magrittr::not(Object %in% ExcludeObs)) %>%
    dplyr::pull(Object)

  AllObjs <- AllObjs %>%
    purrr::map(
      .f = ~{
        Obj <- get(.x, envir = .GlobalEnv)
        if(class(Obj)[1] == "SpatRaster") {
          suppressWarnings(terra::wrap(Obj))
        } else {
          Obj
        }
      }) %>%
    setNames(AllObjs)

  FF2 <- lubridate::now(tzone = "CET") %>%
    purrr::map_chr(
      .f = ~{
        c(lubridate::year(.x), lubridate::month(.x),
          lubridate::day(.x), "__",
          lubridate::hour(.x), lubridate::minute(.x)) %>%
          sapply(stringr::str_pad, width = 2, pad = "0") %>%
          stringr::str_c(collapse = "") %>%
          stringr::str_replace_all("__", "_") %>%
          stringr::str_c("S_", ., collapse = "_")
      })
  IASDT.R::SaveAs(
    InObj = AllObjs, OutObj = FF2,
    OutPath = file.path(Path, paste0(FF2, ".RData")))

  AllObjs %>%
    lapply(pryr::object_size) %>%
    tibble::tibble(Obj = names(.), Size = as.numeric(.)) %>%
    dplyr::mutate(Size = Size / (1024 * 1024), Size = round(Size, 1)) %>%
    dplyr::select(Obj, Size) %>%
    dplyr::arrange(dplyr::desc(Size)) %>%
    return()
}

# |---------------------------------------------------| #
# SaveSessionInfo ----
# |---------------------------------------------------| #
#
#' Save all objects (except functions) of the global environment as list items
#'
#' Save all objects (except functions) of the global environment as list items
#'
#' @param Path Path of where to save the output RData file
#' @param SessionObj List of objects and their sizes (typically a a result of `IASDT::SaveSession` function)
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

SaveSessionInfo <- function(Path = getwd(), SessionObj = NULL) {
  FileName <- lubridate::now(tzone = "CET") %>%
    purrr::map_chr(
      .f = ~{
        c(lubridate::year(.x), lubridate::month(.x),
          lubridate::day(.x), "__",
          lubridate::hour(.x), lubridate::minute(.x)) %>%
          sapply(stringr::str_pad, width = 2, pad = "0") %>%
          stringr::str_c(collapse = "") %>%
          stringr::str_replace_all("__", "_") %>%
          stringr::str_c("S_", ., collapse = "_")
      }) %>%
    stringr::str_c(Path, "/", ., ".txt")

  capture.output(sessioninfo::session_info(), file = FileName)

  if (magrittr::not(is.null(SessionObj))) {
    capture.output(
      IASDT.R::InfoChunk("Objects in the current session"),
      file = FileName, append = TRUE)
    SessionObj %>%
      print(n = Inf) %>%
      capture.output(file = FileName, append = TRUE, type = "output")
  }

  return(invisible(NULL))
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# AssignIfNotExist ----
# |---------------------------------------------------| #
#
#' Assign a value to a variable, only if not existing in the global environment
#'
#' Assign a value to a variable, only if not existing in the global environment
#' @param Variable Variable name
#' @param Value Value to be assigned
#' @param Env environment to assign value to
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @examples
#' AssignIfNotExist(x, TRUE)
#' print(x)
#'
#' # --------------------------------------------------
#'
#' y <- 10
#'
#' # y exists and thus its value was not changed
#' AssignIfNotExist(y, TRUE)
#' print(y)

AssignIfNotExist <- function(Variable, Value, Env = globalenv()) {

  Variable <- as.character(rlang::ensyms(Variable))

  if (!Variable %in% ls(envir = Env)) {
    assign(x = Variable, value = Value, envir = Env)
  } else {
    "The `{Variable}` object already exists in the environment. Current value is:" %>%
      stringr::str_glue() %>%
      crayon::blue() %>%
      cat()

    print(rlang::env_get(Env, paste0(Variable)))
  }
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# LoadAs ----
# |---------------------------------------------------| #
#
#' Load RData file as specific object; i.e. rename loaded object
#'
#' Load RData file as specific object; i.e. rename loaded object
#' @param File path of file
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @examples
#' File <- system.file("testdata", "culcita_dat.RData", package = "lme4")
#'
#' # ---------------------------------------------------------
#' # loading RData using base library
#' # ---------------------------------------------------------
#' (load(File))
#'
#' ls()
#'
#' tibble::tibble(culcita_dat)
#'
#' rm(culcita_dat)
#'
#' # ---------------------------------------------------------
#' # Loading as custom object name
#' # ---------------------------------------------------------
#' NewObj <- LoadAs(File = File)
#'
#' ls()
#'
#' print(tibble::tibble(NewObj))

LoadAs <- function(File = NA) {
  InFile0 <- load(File)
  if (length(InFile0) == 1) {
    OutFile <- get(paste0(InFile0))
  } else {
    OutFile <- lapply(InFile0, function(x) {
      get(paste0(x))
    })
    names(OutFile) <- InFile0
  }
  return(OutFile)
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# LoadMultiple ----
# |---------------------------------------------------| #
#
#' Load multiple RData files together
#'
#' Load multiple RData files together
#' @param Files vector containing the list of `.RData` files to be loaded
#' @param Verbose Verbose? Only valid when `OneObject = FALSE`. Default = `TRUE`
#' @param OneObject Should the objects loaded as one object, each as a list item? If `FALSE`, then all original object names will be loaded to the global environment. Default = `TRUE`.
#' @param ReturnNames Should the name of the objects loaded be returned? Only valid when `OneObject = FALSE`. Default = `TRUE`.
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' (Files <- system.file("testdata", c("culcita_dat.RData", "gopherdat2.RData"), package = "lme4"))
#'
#' # ---------------------------------------------------
#' # Load multiple *.rd files to one list object
#' # ---------------------------------------------------
#' file.exists(Files)
#' ls()
#'
#' MultiObj <- LoadMultiple(Files = Files, OneObject = TRUE)
#' ls()
#'
#' str(MultiObj, 1)
#'
#' # ---------------------------------------------------
#' # Load multiple *.rd files current environment
#' # ---------------------------------------------------
#' rm("MultiObj")
#' ls()
#'
#' LoadMultiple(Files = Files, OneObject = FALSE)
#'
#' print(c("culcita_dat", "Gdat") %in% ls())
#'
#' str(Gdat, 1)
#'
#' str(culcita_dat, 1)
#'
#' # ---------------------------------------------------
#' # Load multiple *.rd files, one object already exists
#' # ---------------------------------------------------
#' rm("culcita_dat")
#' ls()
#'
#' print(c("culcita_dat", "Gdat") %in% ls())
#'
#' try(LoadMultiple(Files = Files, OneObject = FALSE))
#'
#' @export

LoadMultiple <- function(
    Files = NULL, Verbose = TRUE, OneObject = TRUE, ReturnNames = TRUE) {

  if (!inherits(Files, "character")) {
    cat(paste0(crayon::red("Files should be character object"), "\n"))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }

  if (all(file.exists(Files))) {

    if (OneObject) {
      ListNames <- basename(Files) %>%
        stringr::str_replace_all(".RData", "")

      ObjectsLoaded00 <- lapply(
        X = seq_along(Files),
        FUN = function(y) {
          assign(x = ListNames[y], value = LoadAs(Files[y]))
        }) %>%
        stats::setNames(ListNames)

      return(ObjectsLoaded00)

    } else {

      Env1 <- new.env()
      ObjectsLoaded <- lapply(
        X = seq_along(Files),
        FUN = function(y) {
          load(file = Files[y], envir = Env1)
        })

      if (any(ls(envir = Env1) %in% ls(envir = parent.frame()))) {
        cat(paste0(crayon::red("ERROR: Some of the new object names already exists in the current environment. No files were loaded!"), "\n"), sep = "")
        opt <- options(show.error.messages = FALSE)
        on.exit(options(opt))
        stop()
      } else {

        ObjectsLoaded <- lapply(
          X = seq_along(Files),
          FUN = function(y) {
            load(file = Files[y], envir = parent.frame(3))
          }) %>%
          do.call(what = c)

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
  } else {
    cat(paste0(crayon::red("Some of these files do not exist. No objects were loaded!"), "\n"))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }
  # return(invisible(NULL))
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# ReplaceSpace ----
# |---------------------------------------------------| #
#
#' Replace space with underscore
#'
#' Replace space with underscore
#' @param x string
#' @name ReplaceSpace
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' ReplaceSpace("Genus species")
#'
#' ReplaceSpace("Genus species subspecies")
#' @export

ReplaceSpace <- function(x) {
  stringr::str_replace_all(x, " ", "_")
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# ScrapLinks ----
# |---------------------------------------------------| #
#
#' Extract link texts and urls from a web page
#'
#' Extract link texts and urls from a web page
#' @param url the url
#' @name ScrapLinks
#' @return NULL
#' @importFrom rlang .data
#' @references https://gist.github.com/paulrougieux/e1ee769577b40cd9ed9db7f75e9a2cc2
#' @examples
#' ScrapLinks("https://github.com/")
#' @export

ScrapLinks <- function(url) {
  # Create an html document from the url
  webpage <- xml2::read_html(url) %>%
    rvest::html_nodes("a")

  # Extract the URLs
  url_ <- webpage %>%
    rvest::html_attr("href") %>%
    stringr::str_c(url, ., sep = "") %>%
    stringr::str_replace_all(paste0(url, url), url) %>%
    stringr::str_replace_all("//", "/")

  # Extract the link text
  link_ <- webpage %>%
    rvest::html_text() %>%
    stringr::str_replace_all("\n", "") %>%
    stringr::str_replace_all("\\s+", " ") %>%
    stringr::str_trim()

  tibble::tibble(link_text = link_, url = url_) %>%
    dplyr::arrange(.data$url, .data$link_text) %>%
    dplyr::distinct() %>%
    return()
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# DirCreate ----
# |---------------------------------------------------| #
#
#' Create directory if not existed
#'
#' Create directory if not existed
#' @param Path Path of the folder
#' @param Verbose logical; print a message whether the folder was created or already available. Default: `TRUE`
#' @name DirCreate
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' # create new folder (random name) in the temporary folder
#' Path2Create <- file.path(tempdir(), stringi::stri_rand_strings(1, 5))
#' file.exists(Path2Create)
#'
#' DirCreate(Path2Create)
#' DirCreate(Path2Create)
#' DirCreate(Path2Create, Verbose = FALSE)
#' file.exists(Path2Create)
#'
#' @export

DirCreate <- function(Path, Verbose = TRUE) {
  Path2 <- gsub("\\\\", "/", Path)
  if (dir.exists(Path) && Verbose) {
    CatTime(stringr::str_glue("Path: {crayon::bold(Path2)} - already exists"))
  } else {
    dir.create(Path, recursive = TRUE, showWarnings = FALSE)
    if (Verbose) {
      "Path: {crayon::bold(Path2)} created" %>%
        stringr::str_glue() %>%
        CatTime(Date = TRUE)
    }
  }
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# CatTime ----
# |---------------------------------------------------| #

#' Print text and a time stamp
#'
#' Print text and a time stamp
#' @param Text the text to print: default empty string (print time only)
#' @param NLines number of empty lines after the printing; default: 1
#' @param Date Also print date? Default value `FALSE`
#' @param TZ time zone (default: CET)
#' @param ... other arguments passed to `cat`
#' @name CatTime
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @examples
#' CatTime()
#' CatTime(Date = TRUE)
#'
#' CatTime("Time now")
#' CatTime("Time now", Date = TRUE)

CatTime <- function(Text = "", NLines = 1, Date = FALSE, TZ = "CET", ...) {
  DateFormat <- dplyr::if_else(Date, "%d/%m/%Y %X", "%X")
  Now <- lubridate::now(tzone = TZ)
  if (Text == "") {
    cat(format(Now, DateFormat), ...)
    cat(rep("\n", NLines))
  } else {
    Text <- rlang::quo_name(rlang::enquo(Text))
    cat(paste0(Text, " - ", format(Now, DateFormat)), ...)
    cat(rep("\n", NLines))
  }
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# Valid_URL ----
# |---------------------------------------------------| #
#
#' Check the validity of a URL
#'
#' Check the validity of a URL
#' @param url_in URL path
#' @param t timeout
#' @name Valid_URL
#' @references https://stackoverflow.com/questions/52911812/check-if-url-exists-in-r
#' @return NULL
#' @examples
#' urls <- c("http://www.amazon.com", "http://this.isafakelink.biz", "https://stackoverflow.com")
#' sapply(urls, Valid_URL)
#' @export

Valid_URL <- function(url_in, t = 2) {
  con <- url(url_in)
  check <- suppressWarnings(try(open.connection(con, open = "rt", timeout = t), silent = TRUE)[1])
  suppressWarnings(try(close.connection(con), silent = TRUE))
  ifelse(is.null(check), TRUE, FALSE)
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# GetMode ----
# |---------------------------------------------------| #

#' Calculate the mode of numbers
#'
#' Calculate the mode of numbers
#' @param v vector
#' @name GetMode
#' @references https://www.tutorialspoint.com/r/r_mean_median_mode.htm
#' @return NULL
#' @examples
#' GetMode(c(1:10,1,1,3,3,3,3))
#' @export

GetMode <- function(v) {
  unique_vals <- unique(v)
  unique_vals[which.max(tabulate(match(v, unique_vals)))]
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# List2RData ----
# |---------------------------------------------------| #

#' Split list items to separate RData files
#'
#' Split list items to separate RData files
#' @param List list object
#' @param Prefix prefix string
#' @param Dir output directory
#' @param Overwrite overwrite existing files?
#' @name List2RData
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @examples
#' # split iris data by species name
#' iris2 <- iris %>%
#'   tibble::tibble() %>%
#'   split(~Species)
#'
#' str(iris2, 1)
#'
#' (TMP_Folder <- file.path(tempdir(), stringi::stri_rand_strings(1, 5)))
#' list.files(TMP_Folder)
#'
#' List2RData(iris2, Dir = TMP_Folder)
#' list.files(TMP_Folder)

List2RData <- function(List, Prefix = "", Dir = getwd(), Overwrite = FALSE) {
  if (is.null(names(List))) {
    stop()
  } else {

    seq_along(List) %>%
      lapply(
        function(x) {

          if (Prefix != "") {
            Prefix <- paste0(Prefix, "_")
          }
          Name <- names(List[x])
          if (!dir.exists(Dir)) {
            DirCreate(Dir)
          }
          File <- paste0(Dir, "/", Prefix, Name, ".RData")

          if (file.exists(File)) {

            if (Overwrite) {
              assign(x = Name, value = List[[x]])
              save(list = Name, file = File)
            } else {
              "\n\nFile: {File} already exists" %>%
                stringr::str_glue() %>%
                cat()
            }
          } else {
            assign(x = Name, value = List[[x]])
            save(list = Name, file = File)
          }
        }
      ) %>%
      invisible()
  }
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# SplitVector ----
# |---------------------------------------------------| #

#' Split a vector into smaller chunks
#'
#' Split a vector into smaller chunks
#' @param Vector vector to split
#' @param NSplit number of splits
#' @param Prefix prefix string for the name of the resulted list
#' @name SplitVector
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' SplitVector(1:100, 3)
#'
#' # -------------------------------------------
#'
#' SplitVector(1:100, 2, "T")
#'
#' # -------------------------------------------
#'
#' \dontrun{
#' SplitVector(1:100)
#' }
#' @export

SplitVector <- function(Vector = NULL, NSplit = NULL, Prefix = "Chunk") {
  if (inherits(Vector, "NULL") || inherits(NSplit, "NULL")) {
    stop("Vector and NSplit parameters can not be NULL")
  }
  rep(1:NSplit, length.out = length(Vector)) %>%
    sort() %>%
    as.factor() %>%
    split(x = Vector) %>%
    stats::setNames(paste0(Prefix, "_", 1:NSplit))
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# SplitDF2Chunks ----
# |---------------------------------------------------| #

#' Split data.frame into smaller chunks
#'
#' Split data.frame into smaller chunks
#' @param DF dataframe to split
#' @param ChunkSize number of rows per chunk
#' @param NChunks number of chunks (only if `ChunkSize` not provided)
#' @param Prefix prefix
#' @name SplitDF2Chunks
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @examples
#' SplitDF2Chunks(mtcars, ChunkSize = 16)
#'
#' # -------------------------------------------
#'
#' SplitDF2Chunks(mtcars, NChunks = 3)
#'
#' # -------------------------------------------
#'
#' SplitDF2Chunks(mtcars, NChunks = 3, Prefix = "T")

SplitDF2Chunks <- function(
    DF = NULL, ChunkSize = NULL,
    NChunks = NULL, Prefix = "Chunk") {

  if (is.null(DF)) {
    cat(paste0(crayon::red("DF should be a loaded data frame"), "\n"))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }

  DefaultChunkSize <- DefaultNChunks <- FALSE
  if (is.null(NChunks)) {
    NChunks <- min(5, nrow(DF))
    DefaultNChunks <- TRUE
  }

  if (is.null(ChunkSize)) {
    ChunkSize <- ceiling(nrow(DF) / NChunks)
    DefaultChunkSize <- TRUE
  }

  if (any(!is.numeric(ChunkSize), !(ChunkSize > 1))) {
    cat(paste0(crayon::red("ChunkSize should be a numeric vector of length of 1"), "\n"))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }

  if (all(DefaultChunkSize, DefaultNChunks)) {
    cat(paste0(crayon::green("ChunkSize is not determined by user. The default split into 5 chunks is implemented"), "\n"))
  }


  if (nrow(DF) <= ChunkSize) {
    cat(paste0(crayon::red("ChunkSize is larger than the number of rows in the data frame!\nPlease use a smaller ChunkSize."), "\n"))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }

  DF <- tibble::as_tibble(DF)
  Out <- split(DF, (seq_len(nrow(DF)) - 1) %/% ChunkSize)
  names(Out) <- paste0(Prefix, "_", seq_along(Out))
  Out
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# SaveAs ----
# |---------------------------------------------------| #

#' Save an object with a different name
#'
#' Save an object with a different name
#' @param InObj input object
#' @param OutObj output object
#' @param OutPath save path
#' @name SaveAs
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @examples
#' TMP_Folder <- file.path(tempdir(), stringi::stri_rand_strings(1, 5))
#' DirCreate(TMP_Folder, Verbose = FALSE)
#'
#' # save iris data in `iris2.RData` with `iris2` object name
#' SaveAs(iris, "iris2", file.path(TMP_Folder, "iris2.RData"))
#'
#' list.files(TMP_Folder, pattern = "^.+.RData")
#'
#' (load(file.path(TMP_Folder, "iris2.RData")))
#'
#' tibble::tibble(iris2)

SaveAs <- function(InObj, OutObj, OutPath) {
  if (inherits(InObj, "character")) {
    InObj <- get(InObj)
  }
  OutObj <- eval(OutObj)
  assign(OutObj, InObj)
  save(list = OutObj, file = OutPath)
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# SaveMultiple ----
# |---------------------------------------------------| #

#' Save multiple objects to their respective `.RData` files
#'
#' Save multiple objects to their respective `.RData` files
#' @param Vars variables to save
#' @param OutFolder output path
#' @param Overwrite overwrite existing files? If file already exist, no files are saved unless `Overwrite` argument is set as `TRUE`
#' @param Prefix String prefix of the output file
#' @name SaveMultiple
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @examples
#' TMP_Folder <- file.path(tempdir(), stringi::stri_rand_strings(1, 5))
#' DirCreate(TMP_Folder)
#'
#' # ----------------------------------------------
#' # Save x1 and x2 to disk
#' # ----------------------------------------------
#' x1 = 10; x2 = 20
#'
#' SaveMultiple(Vars = cc(x1, x2), OutFolder = TMP_Folder)
#'
#' list.files(path = TMP_Folder, pattern = "^.+.RData")
#'
#' (x1Contents <- LoadAs(file.path(TMP_Folder, "x1.RData")))
#'
#' (x2Contents <- LoadAs(file.path(TMP_Folder, "x2.RData")))
#'
#' # ----------------------------------------------
#' # Use Prefix
#' # ----------------------------------------------
#'
#' SaveMultiple(Vars = cc(x1, x2), OutFolder = TMP_Folder, Prefix = "A_")
#'
#' list.files(path = TMP_Folder, pattern = "^.+.RData")
#'
#' # ----------------------------------------------
#' # File exists, no save
#' # ----------------------------------------------
#' try(SaveMultiple(Vars = c("x1", "x2"), OutFolder = TMP_Folder))
#'
#' # ----------------------------------------------
#' # overwrite existing file
#' # ----------------------------------------------
#' x1 = 100; x2 = 200; x3 = 300
#'
#' SaveMultiple(Vars = cc(x1, x2, x3), OutFolder = TMP_Folder, Overwrite = TRUE)
#'
#' (x1Contents <- LoadAs(file.path(TMP_Folder, "x1.RData")))
#'
#' (x2Contents <- LoadAs(file.path(TMP_Folder, "x2.RData")))
#'
#' (x3Contents <- LoadAs(file.path(TMP_Folder, "x3.RData")))

SaveMultiple <- function(
    Vars = NULL, OutFolder = getwd(), Overwrite = FALSE, Prefix = "", Verbose = FALSE) {

  SkipFun <- any(
    is.null(Vars), !inherits(Vars, "character"),
    any(!Vars %in% ls(envir = rlang::caller_env())))

  if (SkipFun) {
    cat(paste0(crayon::red("Vars should be a character vector for names of objects available in the global environment"), "\n"))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }

  GlobalVars <- ls(rlang::caller_env())
  VarsInGlobal <- Vars %>%
    purrr::map_lgl(~.x %in% GlobalVars) %>%
    all() %>%
    magrittr::not()

  if (VarsInGlobal) {
    "Please check that all variables to be saved exist at your global environment" %>%
      crayon::red() %>%
      cat(sep = "\n")

    MissedVars <- setdiff(Vars, GlobalVars)
    cat(crayon::red(paste0("Variable ", MissedVars, " does not exist.\n")))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  } else {
    if (is.null(OutFolder)) {
      OutFolder <- getwd()
    }
    if (!dir.exists(OutFolder)) {
      IASDT.R::DirCreate(OutFolder, Verbose = FALSE)
    }
  }

  FilesExist <- file.path(OutFolder, paste0(Prefix, Vars, ".RData")) %>%
    file.exists() %>%
    any()


  if (FilesExist && !Overwrite) {
    "Some files already exist.\nNo files are saved. Please use overwrite = TRUE" %>%
      crayon::red() %>%
      cat(sep = "\n")
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  }

  purrr::walk(
    .x = Vars,
    .f = ~{
      Val <- get(.x, envir = parent.frame(n = 4))
      assign(rlang::eval_tidy(.x), Val)
      save(
        list = paste0(.x),
        file = file.path(OutFolder, paste0(Prefix, .x, ".RData")))
    })


  if (all(file.exists(file.path(OutFolder, paste0(Prefix, Vars, ".RData"))))) {
    if (Verbose) {
      cat(paste0(crayon::blue("All files are saved to disk", crayon::red(OutFolder), "successfully."), "\n"))
    }
  } else {
    cat(paste0(crayon::red("Some files were not saved to disk!\nplease check again"), "\n"))
  }
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# cc ----
# |---------------------------------------------------| #

#' Concatenate without quotes
#'
#' Concatenate without quotes
#' @param ... one or more string to concatenate
#' @name cc
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @examples
#' cc(A, B, C)
#'
#' # -------------------------------------------
#'
#' cc(A, B, "A and B")

cc <- function(...) {
  rlang::ensyms(...) %>%
    as.character() %>%
    stringr::str_remove_all("`|`")
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# LoadPackages ----
# |---------------------------------------------------| #

#' Load package silently (+ install missing packages)
#'
#' Load package silently (+ install missing packages)
#' @param ... packages to load / install
#' @param List vector for the name of packages to be loaded
#' @param Verbose print a message of the package name and version
#' @name LoadPackages
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export
#' @examples
#' LoadPackages(raster)
#'
#' LoadPackages(terra, Verbose = TRUE)

LoadPackages <- function(..., List = NULL, Verbose = FALSE) {
  # Packages to load
  PG <- rlang::ensyms(...) %>%
    as.character() %>%
    c(List) %>%
    sort() %>%
    setdiff(as.character(.packages()))

  # packages to install
  InstPack <- utils::installed.packages() %>%
    as.data.frame() %>%
    "["("Package") %>%
    unlist() %>%
    setdiff(x = PG)

  if (length(InstPack) > 0) {
    cat("The following packages will be installed\n")
    InstPack %>%
      stringr::str_c("  >>>>>  ", .) %>%
      stringr::str_c("", collapse = "\n") %>%
      cat()

    InstPack %>%
      purrr::map(
        utils::install.packages, repos = "http://cran.us.r-project.org",
        dependencies = TRUE, quiet = TRUE) %>%
      utils::capture.output(file = nullfile()) %>%
      suppressMessages() %>%
      suppressWarnings()

    cat("\n")
  }

  # load packages
  PG %>%
    sapply(library, character.only = TRUE, quietly = TRUE, verbose = FALSE) %>%
    invisible() %>%
    suppressWarnings() %>%
    suppressMessages()

  if (Verbose & length(PG) > 0) {
    cat("\nThe following packages were loaded:\n")
    PG %>%
      purrr::map_chr(~{
        .x %>%
          utils::packageDescription() %>%
          "$"("Version") %>%
          as.character() %>%
          paste0("  >>>>>  ", .x, ": ", .)
      }) %>%
      stringr::str_c(collapse = "\n") %>%
      cat()
    cat("\n")
  }
  return(invisible(NULL))
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# NDecimals ----
# |---------------------------------------------------| #

#' Number of decimal places in a vector
#'
#' Number of decimal places in a vector
#' @param x vector
#' @name NDecimals
#' @author Ahmed El-Gabbas
#' @return integer; number of decimal places
#' @examples
#' NDecimals("13.45554545")
#'
#' # -------------------------------------------
#'
#' NDecimals(15.01500)
#'
#' NDecimals('15.01500')
#'
#' # -------------------------------------------
#'
#' NDecimals(13.45554545)
#' @export

NDecimals <- function(x) {
  Split <- x %>%
    format(scientific = FALSE) %>%
    stringr::str_split(pattern = "\\.", n = Inf, simplify = TRUE)

  if (length(Split) == 2) {
    Split[, 2] %>%
      nchar() %>%
      as.integer() %>%
      return()
  } else {
    return(0)
  }
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# System ----
# |---------------------------------------------------| #

#' Run bash script depending on the operating system
#'
#' Run bash script depending on the operating system
#' @param command bash command to implement
#' @param RObj whether to make the output of the command an R object
#' @param ... additional arguments to `shell` pr `system` functions
#' @name System
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' # print working directory
#' System("pwd")
#'
#' # first 5 files on the working directory
#' (A <- System("ls | head -n 5"))
#'
#' (A <- System("ls | head -n 5", RObj = FALSE))
#' @export

System <- function(command, RObj = TRUE, ...) {
  # Use shell() in windows OS; system() for linux OS
  # The running operating system (make also available in the global environment outside of the function)

  if (IASDT.R::CurrOS() == "Windows") {
    Out <- shell(cmd = command, intern = RObj, ...)
  }
  if (IASDT.R::CurrOS() == "Linux") {
    Out <- system(command = command, intern = RObj, ...)
  }
  return(Out)
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# CatPipe ----
# |---------------------------------------------------| #

#' print a message with time in the middle of the pipe
#'
#' print a message with time in the middle of the pipe
#' @param x input object to pipe
#' @param message printed message
#' @name CatPipe
#' @author Ahmed El-Gabbas
#' @return NULL
#' @keywords internal
#' @noRd

CatPipe <- function(x, message = NULL) {
  CatTime(message)
  return(x)
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# lapply_ ----
# |---------------------------------------------------| #

#' Apply a Function over a List or Vector (with no return value)
#'
#' Apply a Function over a List or Vector (with no return value)
#' @param X vector
#' @param FUN function
#' @param Silent should the output be silenced?
#' @param ... additional arguments to lapply
#' @name lapply_
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' par(mfrow = c(1,2), oma = c(0.25, 0.25, 0.25, 0.25), mar = c(3,3,3,1))
#' lapply(list(x = 100:110, y = 110:120), function(V) {
#'     plot(V, las = 1, main = "lapply")
#' })
#'
#' # -------------------------------------------
#'
#' par(mfrow = c(1,2), oma = c(0.25, 0.25, 0.25, 0.25), mar = c(3,3,3,1))
#' lapply_(list(x = 100:110, y = 110:120), function(V) {
#'     plot(V, las = 1, main = "lapply_")
#' })
#' @export

lapply_ <- function(X, FUN, Silent = TRUE, ...) {
  if (Silent) {
    invisible(lapply(X = X, FUN = FUN, ...))
  } else {
    lapply(X = X, FUN = FUN, ...)
  }
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# sapply_ ----
# |---------------------------------------------------| #

#' Apply a Function over a List or Vector (with no return value)
#'
#' Apply a Function over a List or Vector (with no return value)
#' @param X vector
#' @param FUN function
#' @param simplify should be the output be simplified
#' @param ... additional arguments to sapply
#' @name sapply_
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' par(mfrow = c(1,2), oma = c(0.25, 0.25, 0.25, 0.25), mar = c(3,3,3,1))
#' sapply(list(x = 100:110, y = 110:120), function(V){plot(V, las = 1, main = "sapply")})
#'
#' # -------------------------------------------
#'
#' # nothing returned or printed, only the plotting
#' par(mfrow = c(1,2), oma = c(0.25, 0.25, 0.25, 0.25), mar = c(3,3,3,1))
#' sapply_(list(x = 100:110, y = 110:120), function(V){plot(V, las = 1, main = "sapply_")})
#' @export

sapply_ <- function(X, FUN, simplify = TRUE, ...) {
  invisible(
    sapply(X = X, FUN = FUN, simplify = simplify, ...)
  )
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# sort_ ----
# |---------------------------------------------------| #

#' Sort alphanumeric strings
#'
#' Sort alphanumeric strings
#' @param x Vector to be sorted.
#' @param decreasing logical. Should the sort be increasing or decreasing? Note that descending=TRUE reverses the meanings of na.last and blanks.last. Default: `FALSE`
#' @param na.last for controlling the treatment of NA values. If TRUE, missing values in the data are put last; if FALSE, they are put first; if NA, they are removed. Default: `TRUE`
#' @param blank.last for controlling the treatment of blank values. If TRUE, blank values in the data are put last; if FALSE, they are put first; if NA, they are removed. Default: `FALSE`
#' @param numeric.type either "decimal" (default) or "roman". Are numeric values represented as decimal numbers (numeric.type="decimal") or as Roman numerals (numeric.type="roman")?
#' @param roman.case one of "upper", "lower", or "both". Are roman numerals represented using only capital letters ('IX') or lower-case letters ('ix') or both?
#' @name sort_
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' # example code
#' (AA <- paste0("V", 1:12))
#'
#' sort(AA)
#'
#' sort_(AA)
#'
#' @export

sort_ <- function(
    x, decreasing = FALSE, na.last = TRUE,
    blank.last = FALSE, numeric.type = c("decimal", "roman"),
    roman.case = c("upper", "lower", "both")) {

  gtools::mixedsort(
    x = x, decreasing = decreasing, na.last = na.last,
    blank.last = blank.last,
    numeric.type = numeric.type,
    roman.case = roman.case)
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# CurrOS ----
# |---------------------------------------------------| #

#' Current operating system
#'
#' Current operating system
#' @name CurrOS
#' @author Ahmed El-Gabbas
#' @return String for the current operating system
#' @examples
#' CurrOS()
#'
#' @export

CurrOS <- function() {
  as.character(Sys.info()["sysname"])
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# ht ----
# |---------------------------------------------------| #

#' Print head and tail of data frame
#'
#' Print head and tail of data frame
#' @name ht
#' @author Ahmed El-Gabbas
#' @return NULL
#' @param DF data frame to print
#' @param NRows N umber of rows to print at the top and bottom of data frame
#' @examples
#' ht(mtcars)
#'
#' # -------------------------------------------
#'
#' ht(mtcars, 2)
#'
#' # -------------------------------------------
#'
#' ht(mtcars, 6)
#' @export

ht <- function(DF, NRows = 5) {
  DF %>%
    data.table::data.table() %>%
    print(topn = NRows)
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# AddMissingCols ----
# |---------------------------------------------------| #

#' Add missing columns to data frame
#'
#' Add missing columns to data frame
#' @name AddMissingCols
#' @author Ahmed El-Gabbas
#' @return NULL
#' @param DT data frame
#' @param FillVal value to be used: default: NA_character
#' @param ... list of column names to add
#' @export
#' @examples
#' mtcars %>%
#'  dplyr::select(1:3) %>%
#'  AddMissingCols(FillVal = NA_character_, A, B, C) %>%
#'  AddMissingCols(FillVal = as.integer(10), D)
#'
#' # -------------------------------------------
#'
#' AddCols <- c("Add1", "Add2")
#' mtcars %>%
#'  dplyr::select(1:3) %>%
#'  AddMissingCols(FillVal = NA_real_, AddCols)

AddMissingCols <- function(DT, FillVal = NA_character_, ...) {
  Cols <- as.character(rlang::ensyms(...))

  if (any(Cols %in% ls(envir = parent.env(rlang::caller_env())))) {
    Cols <- get(Cols, envir = parent.env(rlang::caller_env()))
  }

  Cols2Add <- setdiff(Cols, names(DT))

  Add_DF <- rep(FillVal, length(Cols2Add)) %>%
    matrix(nrow = 1) %>%
    as.data.frame() %>%
    stats::setNames(Cols2Add) %>%
    tibble::as_tibble()

  if (length(Cols2Add) != 0) {
    DT <- tibble::add_column(DT, !!!Add_DF)
  }
  return(tibble::tibble(DT))
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# FileExt ----
# |---------------------------------------------------| #
#
#' Get file extension
#'
#' Get file extension
#' @name FileExt
#' @author Ahmed El-Gabbas
#' @return NULL
#' @param Path File path
#' @examples
#' FileExt("File.doc")
#' FileExt("D:/File.doc")
#'
#' @export

FileExt <- function(Path) {
  Ext <- Path %>%
    basename() %>%
    stringr::str_split("\\.", simplify = TRUE)
  Ext[length(Ext)]
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# FunctionsInPackage ----
# |---------------------------------------------------| #
#
#' List of functions in a package
#'
#' List of functions in a package
#' @name FunctionsInPackage
#' @author Ahmed El-Gabbas
#' @return NULL
#' @param Package name of the package
#' @export
#' @examples
#' FunctionsInPackage("raster")

FunctionsInPackage <- function(Package) {
  library(package = Package, character.only = TRUE)
  ls(paste0("package:", Package))
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# LoadedPackages ----
# |---------------------------------------------------| #
#
#' List of loaded packages
#'
#' List of loaded packages
#' @name LoadedPackages
#' @examples
#' LoadedPackages()
#'
#' LoadPackages(sf)
#'
#' LoadedPackages()
#'
#' @export

LoadedPackages <- function() {
  (.packages())
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# ReloadPackage ----
# |---------------------------------------------------| #
#
#' Reload an R package
#'
#' Reload an R package
#' @name ReloadPackage
#' @param Package The name of packages to reload
#' @examples
#' LoadPackages(sf)
#'
#' ReloadPackage("sf")
#'
#' @export

ReloadPackage <- function(Package = NULL) {
  if (is.null(Package)) {
    stop()
  }

  if (!Package %in% LoadedPackages()) {
    library(package = Package, character.only = TRUE)
  } else {
    PackagesFolders <- paste0(.libPaths(), "/")
    PackagesFolders <- c(outer(PackagesFolders, Package, FUN = "paste0"))
    PackagesFolders <- PackagesFolders[file.exists(PackagesFolders)]
    PackagesFolders <- gsub(pattern = "/", replacement = "//", x = PackagesFolders)
    lapply(PackagesFolders, function(x) {
      devtools::reload(x, quiet = FALSE)
    })
  }

  invisible(NULL)
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# AddImg2Plot ----
# |---------------------------------------------------| #
#
#' Add image to plot
#'
#' Add image to plot
#'
#' @name AddImg2Plot
#' @references [Click here](https://stackoverflow.com/questions/27800307/)
#' @export
#' @param obj an image file imported as an array (e.g. `png::readPNG`, `jpeg::readJPEG`)
#' @param x mid x coordinate for image
#' @param y mid y coordinate for image
#' @param width width of image (in x coordinate units)
#' @param interpolate (passed to `graphics::rasterImage`) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. Default = `TRUE`.
#' @examples
#' LoadPackages(png)
#' myurl <- "https://upload.wikimedia.org/wikipedia/commons/e/e1/Jupiter_%28transparent%29.png"
#' z <- tempfile()
#' download.file(myurl, z, mode="wb")
#' pic <- readPNG(z)
#' file.remove(z) # cleanup
#'
#' image(volcano)
#' AddImg2Plot(pic, x = 0.3, y = 0.5, width = 0.2)
#' AddImg2Plot(pic, x = 0.7, y = 0.7, width = 0.2)
#' AddImg2Plot(pic, x = 0.7, y = 0.2, width = 0.1)

AddImg2Plot <- function(
    obj, x = NULL, y = NULL, width = NULL, interpolate = TRUE) {

  if (any(is.null(x), is.null(y), is.null(width))) {
    stop("Must provide args 'x', 'y', and 'width'")
  }

  USR <- graphics::par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- graphics::par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1] / DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width / (USR[2] - USR[1]) * PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi / PIN[2] * (USR[4] - USR[3]) # height in units
  graphics::rasterImage(
    image = obj, xleft = x - (width / 2), xright = x + (width / 2),
    ybottom = y - (HEIu / 2), ytop = y + (HEIu / 2), interpolate = interpolate)
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# AddLine ----
# |---------------------------------------------------| #
#
#' Add a line to the current plot
#'
#' Add a line to the current plot
#'
#' @name AddLine
#' @references [Click here](https://stackoverflow.com/questions/27800307/)
#' @export
#' @param at relative location of where the line should be plotted
#' @param Outer plot the line out of plotting area. Default: `FALSE`
#' @param H Horizontal line? if H == FALSE, vertical line will be added. Default: `TRUE`
#' @param ... Other arguments
#' @examples
#' # Hotizontal line
#' par(oma = c(1, 1, 1, 1), mar = c(3, 3, 1, 1))
#' plot(1:100)
#' AddLine(at = 0.75)
#' AddLine(at = 0.25, Outer = TRUE, lwd = 2)
#' AddLine(at = 0.5, Outer = TRUE, lwd = 2, col = "red")
#'
#' # ---------------------------------------------
#'
#' # Vertical line
#' plot(1:100)
#' AddLine(H = FALSE, at = 0.75)
#' AddLine(H = FALSE, at = 0.25, Outer = TRUE, lwd = 2)
#' AddLine(H = FALSE, at = 0.5, Outer = TRUE, lwd = 2, col = "red")

AddLine <- function(at = NULL, Outer = FALSE, H = TRUE, ...) {
  if (Outer) graphics::par(xpd = TRUE)

  if (H) {
    graphics::abline(h = graphics::grconvertX(at, "npc"), ...)
  } else {
    graphics::abline(v = graphics::grconvertX(at, "npc"), ...)
  }

  if (Outer) {
    graphics::par(xpd = FALSE)
  }
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# ScriptLocation ----
# |---------------------------------------------------| #
#
#' The location of current script
#'
#' The location of current script
#'
#' @name ScriptLocation
#' @references [Click here](https://stackoverflow.com/questions/47044068/)
#' @importFrom rlang .data
#' @export
#' @examples
#'  \dontrun{
#' ScriptLocation()
#' }

ScriptLocation <-  function() {
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col = .data$value, into = c("key", "value"), sep = "=", fill = "right") %>%
    dplyr::filter(.data$key == "--file") %>%
    dplyr::pull(.data$value)
  if (length(this_file) == 0) {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(this_file)
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# Range2NewVal ----
# |---------------------------------------------------| #
#
#' Change values of `vector`/`data.frame`/`raster` within/outside a certain range to another value
#'
#' Change values of `vector`/`data.frame`/`raster` within/outside a certain range to another value
#'
#' @name Range2NewVal
#' @author Ahmed El-Gabbas
#' @param x vector / data.frame / raster object.
#' @param Between Range of values to be changed or kept.
#' @param MoreThan The value above which the original values will be changed. Only applied if `Between` is not used.
#' @param LessThan The value below which the original values will be changed. Only applied if `Between` and `MoreThan` are not used.
#' @param NewVal The new value to be assigned. Default: `NA`.
#' @param InvertSelection Should the selection be inverted? Only valid if `Between` used. If `InvertSelection = TRUE`, then the values not in the range of `Between` argument will be changed. Default: `FALSE`.
#'
#' @export
#' @examples
#' # Vector
#'
#' Range2NewVal(x =  1:10, Between = c(5, 8), NewVal = NA)
#'
#' Range2NewVal(x =  1:10, Between = c(5, 8), NewVal = NA, InvertSelection = TRUE)
#'
#' # ---------------------------------------------
#'
#' # tibble
#'
#' iris %>%
#'  tibble::as_tibble() %>%
#'  dplyr::slice_head(n = 50) %>%
#'  dplyr::select(-Sepal.Length, -Petal.Length, -Petal.Width) %>%
#'  dplyr::mutate(
#'    Sepal.Width.New = Range2NewVal(
#'        x = Sepal.Width, Between = c(3, 3.5), NewVal = NA, InvertSelection = FALSE),
#'    Sepal.Width.Rev = Range2NewVal(
#'        x = Sepal.Width, Between = c(3, 3.5), NewVal = NA, InvertSelection = TRUE)) %>%
#'  dplyr::arrange(-Sepal.Width) %>%
#'  print(n = 50)
#'
#' # ---------------------------------------------
#'
#' # raster
#'
#' LoadPackages(raster)
#'
#' RRR <- system.file("external/test.grd", package = "raster") %>%
#'     raster::raster()
#'
#' RRR2 <- Range2NewVal(x = RRR, LessThan = 500, NewVal = NA)
#' RRR3 <- Range2NewVal(x = RRR, MoreThan = 500, NewVal = NA)
#' par(mar = c(0.5, 0.5, 3, 3))
#' plot(raster::stack(RRR, RRR2, RRR3), nr = 1, main = c("Original", "<500 to NA", ">500 to NA"))
#'
#' RRR2 <- Range2NewVal(x = RRR, Between = c(1000, 1800), NewVal = 1800, InvertSelection = FALSE)
#' RRR3 <- Range2NewVal(x = RRR, Between = c(1000, 1800), NewVal = 1800, InvertSelection = TRUE)
#' plot(stack(RRR>=1000, RRR2, RRR3), nr = 1, main = c(">1000 ?", "<500 to NA", ">500 to NA"))

Range2NewVal <- function(
    x = NULL, Between = NULL,
    MoreThan = NULL, LessThan = NULL,
    NewVal = NA, InvertSelection = FALSE) {

  if (!is.null(Between)) {
    if (length(Between) != 2) stop()
    Min <- Between[1]
    Max <- Between[2]
    if (inherits(x, "RasterLayer")) {
      X1 <- X2 <- x
      X1[X1 >= Max] <- NA
      X1[!is.na(X1)] <- 1
      X2[X2 <= Min] <- NA
      X2[!is.na(X2)] <- 1
      X3 <- sum(X1, X2, na.rm = TRUE)

      if (InvertSelection) {
        x[X3 == 1] <- NewVal
      } else {
        x[X3 == 2] <- NewVal
      }
    } else {

      if (is.null(Max)) Max <- max(x)
      if (is.null(Min)) Min <- min(x)
      if (Max <= Min) stop()

      if (InvertSelection) {
        x[!(x >= Min & x <= Max)] <- NewVal
      } else {
        x[x >= Min & x <= Max] <- NewVal
      }
    }
  } else {

    if (!is.null(MoreThan)) {
      x[x > MoreThan] <- NewVal
    } else {
      if (!is.null(LessThan)) x[x < LessThan] <- NewVal
    }
  }
  x
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# ClearConsole ----
# |---------------------------------------------------| #
#
#' Clear the console
#'
#' Clear the console
#' This function is a lazy equivalent of `cat("\014")`
#'
#' @name ClearConsole
#' @export
#' @examples
#' \dontrun{
#' ClearConsole()
#' }

ClearConsole <- function() {
  cat("\014")
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# RequireMultiple ----
# |---------------------------------------------------| #

#' Load (or install) multiple packages at once
#'
#' @name RequireMultiple
#' @param ... Packages to load or install
#' @param Reload Should the already loaded packages re-loaded? Default: `FALSE`
#' @param InstallMissing Should the missing packages be installed and loaded? Default: `TRUE`
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' # Currently loaded packages
#' (P1 <- LoadedPackages())
#'
#' RequireMultiple(dismo, MASS, mgcv)
#'
#' # Loaded packages after implementing the function
#' (P2 <- LoadedPackages())
#'
#' # Which pacakages were loaded?
#' setdiff(P2, P1)

RequireMultiple <- function(..., Reload = FALSE, InstallMissing = TRUE) {
  options(warnings = -1)

  varnames <- sapply(substitute(list(...))[-1], deparse) %>%
    stringr::str_replace_all(pattern = '\"', "")

  sapply(varnames, function(x) {

    if (x %in% rownames(utils::installed.packages())) {

      if (!x %in% LoadedPackages()) {
        suppressPackageStartupMessages(suppressMessages(suppressWarnings(
          require(x, character.only = TRUE))))
        cat(crayon::red(">>"), "Library", crayon::blue(x), "was loaded successfully.\n")

      } else {

        if (Reload == TRUE) {
          suppressPackageStartupMessages(suppressMessages(suppressWarnings({
            ReloadPackage(Package = x)
            cat(crayon::red(">>"), "Library", crayon::blue(x), "was already loaded (re-loaded).\n")

          })))
        } else {
          cat(crayon::red(">>"), "Library", crayon::blue(x), "was already loaded (not re-loaded).\n")

        }
      }
    } else {

      if (InstallMissing) {
        suppressPackageStartupMessages(suppressMessages(suppressWarnings({
          utils::install.packages(
            x, quiet = TRUE, verbose = FALSE,
            repos = "https://cloud.r-project.org")
        })))

        if (x %in% rownames(utils::installed.packages())) {
          suppressPackageStartupMessages(suppressMessages(suppressWarnings({
            require(x, character.only = TRUE)
          })))
          cat(crayon::yellow(">>>>>>"), "Library", crayon::blue(x), "was installed and loaded. \n")
        } else {
          cat(crayon::yellow(">>>>>>"), "Library", crayon::blue(x), "can not be installed. \n")
        }
      } else {
        cat(crayon::red(">>"), "Library", crayon::blue(x), "is not installed. ")
      }
    }
  })
  return(invisible(NULL))
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# CheckRStudioVersion ----
# |---------------------------------------------------| #

#' Check if `RStudio` should be updated
#'
#' Check if `RStudio` should be updated
#'
#' @name CheckRStudioVersion
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' \dontrun{
#' CheckRStudioVersion()
#' }

CheckRStudioVersion <- function() {
  XPath <- ".flex-inhe:nth-child(8)"
  OnlineVersion <- "https://posit.co/download/rstudio-desktop/" %>%
    xml2::read_html() %>%
    rvest::html_node(XPath) %>%
    rvest::html_text2() %>%
    stringr::str_remove_all("RStudio-|.exe") %>%
    stringr::str_replace_all("-", ".")

  InstalledVersion <- rstudioapi::versionInfo() %>%
    "[["("long_version") %>%
    stringr::str_replace_all("\\+", "\\.")

  if (identical(OnlineVersion, InstalledVersion) == FALSE) {
    cat(
      crayon::blue(
        "R-Studio version:",
        crayon::red(crayon::bold(OnlineVersion)),
        "is available.\nInstalled R-studio version:",
        crayon::red(crayon::bold(InstalledVersion)),
        "\nPlease consider updating R-Studio.\n"))
  } else {
    cat(
      crayon::blue(
        "You are using the most recent version of R-Studio: v",
        crayon::red(crayon::bold(InstalledVersion)), ".",
        sep = ""))
  }
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# CheckQuartoVersion ----
# |---------------------------------------------------| #

#' Check if `Quarto` should be updated
#'
#' Check if `Quarto` should be updated
#' @name CheckQuartoVersion
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' CheckQuartoVersion()

CheckQuartoVersion <- function() {

  OnlineVersion <- "https://github.com/quarto-dev/quarto-cli/releases/" %>%
    xml2::read_html() %>%
    rvest::html_nodes(".Link--primary") %>%
    rvest::html_text2() %>%
    stringr::str_remove_all("v") %>%
    gtools::mixedsort(decreasing = TRUE) %>%
    "["(1)

  InstalledVersion <- system("quarto --version", intern = TRUE)

  if (identical(OnlineVersion, InstalledVersion) == FALSE) {
    cat(
      crayon::blue(
        "Quarto version:",
        crayon::red(crayon::bold(OnlineVersion)),
        "is available.\nInstalled Quarto version:",
        crayon::red(crayon::bold(InstalledVersion)),
        "\nPlease consider updating Quarto.\n"))
  } else {
    cat(
      crayon::blue(
        "You are using the most recent version of Quarto: v",
        crayon::red(crayon::bold(InstalledVersion)), ".",
        sep = ""))
  }
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# AllObjSizes ----
# |---------------------------------------------------| #

#' Size of objects in memory
#'
#' This function print the size allocated by any of loaded objects into `R`
#' @name AllObjSizes
#' @param GreaterThan Only show objects with size > specified value in MB. Default: 0, which means show all variables size
#' @author Ahmed El-Gabbas
#' @importFrom rlang .data
#' @export
#' @examples
#' AA1 <<- rep(1:1000, 10000)
#' AA2 <<- rep(1:1000, 100)
#'
#' AllObjSizes()
#'
#' AllObjSizes(GreaterThan = 1)
#'
#' AllObjSizes(GreaterThan = 50)

AllObjSizes <- function(GreaterThan = 0) {
  AllVars <- ls(envir = .GlobalEnv)

  if (length(AllVars) == 0) {
    cat("No Objects are available in the global environment!\n")
  } else {
    AllVarsSize <- AllVars %>%
      sapply(
        FUN = function(x) {
          pryr::object_size(get(x)) / (1024 * 1024)
        }) %>%
      dplyr::tibble() %>%
      stats::setNames("Size") %>%
      dplyr::mutate(
        Vars = AllVars,
        Size = as.numeric(.data$Size),
        Size = round(.data$Size, 4),
        Percent = 100 * (.data$Size / sum(.data$Size)),
        Percent = round(.data$Percent, 2),
        Percent = paste0(.data$Percent, " %")) %>%
      dplyr::arrange(dplyr::desc(.data$Size)) %>%
      dplyr::select(tidyselect::all_of(c("Vars", "Size", "Percent")))

    if (!is.na(GreaterThan)) {
      AllVarsSize <- dplyr::filter(AllVarsSize, .data$Size >= GreaterThan)
      if (nrow(AllVarsSize) > 0) {
        cat(crayon::blue(
          "---------------------------------------------\n",
          crayon::bold(nrow(AllVarsSize)),
          " Object(s) fulfill the criteria.\n",
          "---------------------------------------------\n",
          sep = ""),
          sep = "")
        print(AllVarsSize, row.names = FALSE, n = Inf)
        cat(crayon::blue(
          "Object sizes are in MB.\n",
          "---------------------------------------------\n", sep = ""),
          sep = "")
      } else {
        cat(crayon::red(
          paste0("No variables have Size > ",
                 GreaterThan, " MB\n")), sep = "")
      }
    } else {
      print(AllVarsSize, row.names = FALSE, n = Inf)
      cat(crayon::green("Object sizes are in MB"), sep = "")
    }
  }
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# KeepOnly ----
# |---------------------------------------------------| #

#' Keep only certain objects in memory, all other objects will be removed
#'
#' Keep only certain objects in memory, all other objects will be removed
#'
#' @name KeepOnly
#' @param Obj character vector for objects to be kept in memory
#' @param Verbose Should the names of kept and removed variables printed? Default: `TRUE`.
#' @author Ahmed El-Gabbas
#' @export
#' @examples
#' A <- B <- C <- 15
#' ls()
#'
#' KeepOnly("A")
#'
#' ls()
#' rm(list = ls())
#'
#'
#' A <- B <- C <- 15
#' KeepOnly(c("A","B"))
#' ls()

KeepOnly <- function(Obj = NULL, Verbose = TRUE) {
  if (is.null(Obj) || length(Obj) == 0) {
    stop()
  }
  AllObjects <- ls(pos = parent.frame())
  RemObjects <- setdiff(AllObjects, Obj)
  if (Verbose) {
    cat(crayon::red(paste0("Removed Variables (", length(RemObjects), "): ")), crayon::blue(paste0(seq_along(RemObjects), ":", RemObjects, collapse = " ||  ")), sep = "")
  }
  rm(list = RemObjects, pos = parent.frame())
  if (Verbose) {
    cat(crayon::red(paste0("\nKept Variables (", length(Obj), "): ")), crayon::blue(paste0(seq_along(Obj), ":", Obj, collapse = " ||  ")), "\n", sep = "")
  }
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#
# |---------------------------------------------------| #
# SourceSilent ----
# |---------------------------------------------------| #

#' Silently source R script
#'
#' Silently source R script
#'
#' @name SourceSilent
#' @param File path of the file to be sourced
#' @param ... additional arguments passed to `source` function
#' @param Messages Show messages; default: `TRUE`
#' @param Warnings Show warnings; default: `TRUE`
#' @author Ahmed El-Gabbas
#' @export

SourceSilent <- function(File, Messages = TRUE, Warnings = TRUE, ...) {

  if (Messages && Warnings) {
    File %>%
      source(...) %>%
      utils::capture.output(file = nullfile())
  }

  if (!Messages && !Warnings) {
    File %>%
      source(...) %>%
      utils::capture.output(file = nullfile()) %>%
      suppressMessages() %>%
      suppressWarnings()
  }

  if (Messages && !Warnings) {
    File %>%
      source(...) %>%
      utils::capture.output(file = nullfile()) %>%
      suppressWarnings()
  }

  if (!Messages && Warnings) {
    File %>%
      source(...) %>%
      utils::capture.output(file = nullfile()) %>%
      suppressMessages()
  }
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# NUnique ----
# |---------------------------------------------------| #

#' Number of unique values for all columns of a data frame
#'
#' Number of unique values for all columns of a data frame
#'
#' @name NUnique
#' @param Data input data
#' @references https://stackoverflow.com/questions/22196078/count-unique-values-for-every-column
#' @export
#' @examples
#' NUnique(mtcars)
#'
#' NUnique(iris)

NUnique <- function(Data) {
  Data %>%
    dplyr::summarise(
      dplyr::across(tidyselect::everything(), dplyr::n_distinct)) %>%
    tidyr::pivot_longer(tidyselect::everything()) %>%
    stats::setNames(c("Variable", "NUnique")) %>%
    dplyr::arrange(NUnique)
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# Package_RemoteSha ----
# |---------------------------------------------------| #

#' get remote sha of R packages
#'
#' get remote sha of R packages
#'
#' @name Package_RemoteSha
#' @param ... name of one or more R packages
#' @export
#' @examples
#' # Package_RemoteSha(IASDT.R, devtools)

Package_RemoteSha <- function(...) {
  Pk <- rlang::ensyms(...) %>%
    as.character()
  Pk %>%
    purrr::map_chr(~{
      pak::lib_status() %>%
        dplyr::filter(package == .x) %>%
        dplyr::pull("remotesha")
    }) %>%
    stats::setNames(Pk)
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# git_log ----
# |---------------------------------------------------| #

#' print detailed `git log` of the git repo located in the current working directory
#'
#' print detailed `git log` of the git repo located in the current working directory
#'
#' @name git_log
#' @export
#' @examples
#' git_log()

git_log <- function() {
  'git log --graph --pretty=format:"%Cred%h%Creset -%C(yellow)%d%Creset %s %Cgreen(%cr) %C(bold blue)<%an>%Creset" --abbrev-commit' %>%
    IASDT.R::System() %>%
    cat(sep = "\n")
}
