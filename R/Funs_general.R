# |---------------------------------------------------| #
# catSep ----
# |---------------------------------------------------| #

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
#' catSep()
#' catSep(2)
#' catSep(2,2,3)
#' @export

catSep <- function(Rep = 1, Extra1 = 0, Extra2 = 0, Char = "-", CharReps = 50) {
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

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# AssignIfNotExist ----
# |---------------------------------------------------| #

#' Assign a value to a variable, only if not existing in the global environment
#'
#' Assign a value to a variable, only if not existing in the global environment
#' @param Variable Variable name
#' @param Value Value to be assigned
#' @author Ahmed El-Gabbas
#' @return NULL
#' @examples
#' AssignIfNotExist(x, TRUE)
#' print(x)
#'
#' y <- 10
#' AssignIfNotExist(y, TRUE) # y exists and thus its value was not changed
#' print(y)
#' @export

AssignIfNotExist <- function(Variable, Value) {
  Variable <- rlang::ensyms(Variable) %>%
    as.character()
  if (!Variable %in% ls(envir = globalenv())) {
    assign(x = Variable, value = Value, envir = globalenv())
  }
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# LoadAs ----
# |---------------------------------------------------| #

#' Load RData file as specific object; i.e. rename loaded object
#'
#' Load RData file as specific object; i.e. rename loaded object
#' @param File path of file
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

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
#' \dontrun{
#' (Files <- system.file("testdata", c("culcita_dat.RData", "gopherdat2.RData"), package = "lme4"))
#'
#' # ---------------------------------------------------
#'
#' # Load multiple *.rd files to one list object
#' print(Files)
#' file.exists(Files)
#' MultiObj <- LoadMultiple(Files = Files, OneObject = TRUE)
#' MultiObj
#'
#' # ---------------------------------------------------
#'
#' # Load multiple *.rd files current environment
#' rm("MultiObj")
#' LoadMultiple(Files = Files, OneObject = FALSE)
#'
#' print(c("culcita_dat", "Gdat") %in% ls(globalenv()))
#' print(Gdat)
#' print(culcita_dat)
#'
#' # ---------------------------------------------------
#'
#' # Load multiple *.rd files, one object already exists
#' rm("culcita_dat")
#' ls()
#' print(c("culcita_dat", "Gdat") %in% ls(globalenv()))
#' LoadMultiple(Files = Files, OneObject = FALSE)
#'}
#' @export

LoadMultiple <- function(
    Files = NULL, Verbose = TRUE, OneObject = TRUE, ReturnNames = TRUE) {

  if (!methods::is(Files, "character")) {
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
            load(file = Files[y], envir = .GlobalEnv)
          }) %>%
          do.call(what = c)

        if (Verbose) {
          cat(
            paste0(crayon::blue(
              "Object:", crayon::red(ObjectsLoaded),
              "was loaded successfully.")),
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

#' Replace space with underscore
#'
#' Replace space with underscore
#' @param x string
#' @name ReplaceSpace
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

ReplaceSpace <- function(x) {
  stringr::str_replace(x, " ", "_")
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# scraplinks ----
# |---------------------------------------------------| #

#' Extract link texts and urls from a web page
#'
#' Extract link texts and urls from a web page
#' @param url the url
#' @name scraplinks
#' @author Ahmed El-Gabbas
#' @return NULL
#' @description
#' source: https://gist.github.com/paulrougieux/e1ee769577b40cd9ed9db7f75e9a2cc2
#' @examples
#' head(scraplinks("https://github.com/"))
#'
#' @export

scraplinks <- function(url) {
  # Create an html document from the url
  webpage <- xml2::read_html(url)
  # Extract the URLs
  url_ <- webpage %>%
    rvest::html_nodes("a") %>%
    rvest::html_attr("href")
  # Extract the link text
  link_ <- webpage %>%
    rvest::html_nodes("a") %>%
    rvest::html_text()
  return(tibble::tibble(link = link_, url = url_))
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# DirCreate ----
# |---------------------------------------------------| #

#' Create directory if not existed
#'
#' Create directory if not existed
#' @param Path Path of the folder
#' @name DirCreate
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

DirCreate <- function(Path) {
  if (dir.exists(Path)) {
    CatTime(stringr::str_glue("Path: {crayon::bold(Path)} - already exists"))
  } else {
    dir.create(Path, recursive = TRUE, showWarnings = FALSE)
    CatTime(stringr::str_glue("Path: {crayon::bold(Path)} created"))
  }
}


# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# CatTime ----
# |---------------------------------------------------| #

#' Print text and a time stamp
#'
#' Print text and a time stamp
#' @param Text the text to print
#' @param NLines number of empty lines after the printing; default: 1
#' @param ... other arguments passed to `cat`
#' @name CatTime
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

CatTime <- function(Text = NULL, NLines = 1, ...) {
  if (inherits(Text, "NULL")) {
    cat(format(Sys.time(), "%X"), ...)
    cat(rep("\n", NLines))
  } else {
    Text <- rlang::quo_name(rlang::enquo(Text))
    cat(paste0(Text, " - ", format(Sys.time(), "%X")), ...)
    cat(rep("\n", NLines))
  }
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# valid_url ----
# |---------------------------------------------------| #

#' Check the validity of a URL
#'
#' Check the validity of a URL
#' @param url_in URL path
#' @param t timeout
#' @name valid_url
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

valid_url <- function(url_in, t = 2) {
  con <- url(url_in)
  check <- suppressWarnings(try(open.connection(con, open = "rt", timeout = t), silent = TRUE)[1])
  suppressWarnings(try(close.connection(con), silent = TRUE))
  ifelse(is.null(check), TRUE, FALSE)
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# getmode ----
# |---------------------------------------------------| #

#' Calculate the mode of numbers
#'
#' Calculate the mode of numbers
#' @param v vector
#' @name getmode
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

getmode <- function(v) {
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

List2RData <- function(List, Prefix = "", Dir = getwd(), Overwrite = FALSE) {
  if (is.null(names(List))) {
    stop()
  } else {

    seq_len(List) %>%
      lapply(
        function(x) {

          if (Prefix != "") {
            Prefix <- paste0(Prefix, "_")
          }

          Name <- names(List[x])
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
#' @name SplitVector
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

SplitVector <- function(Vector = NULL, NSplit = NULL) {
  rep(1:NSplit, length.out = length(Vector)) %>%
    sort() %>%
    as.factor() %>%
    split(x = Vector) %>%
    stats::setNames(paste0("Chunk", "_", 1:NSplit))
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
#' @param NChunks number of chunks
#' @param Prefix prefix
#' @name SplitDF2Chunks
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

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
  names(Out) <- paste0(Prefix, "_", seq_len(Out))
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

SaveAs <- function(InObj, OutObj, OutPath) {
  if (methods::is(InObj, "character")) {
    InObj <- get(InObj)
  }

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
#' @param Overwrite overwrite existing files
#' @name SaveMultiple
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

SaveMultiple <- function(
    Vars = NULL, OutFolder = getwd(), Overwrite = FALSE) {

  if (any(is.null(Vars), !is.character(Vars),
          any(!Vars %in% ls(envir = globalenv())))) {
    cat(paste0(crayon::red("Vars should be a character vector for names of objects available in the global environment"), "\n"))
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop()
  } else {
    if (!all(sapply(Vars, function(x) {
      exists(x, envir = globalenv())
    }))) {
      # if (!all(exists(Vars, envir = globalenv()))) {
      "Please check that all variables to be saved exist at your global environment" %>%
        crayon::red() %>%
        cat(sep = "\n")
      MissedVars <- Vars[!exists(Vars, envir = globalenv())]
      cat(crayon::red(paste0("Variable ", MissedVars, " does not exist.\n")))
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      stop()
    } else {
      ID <- length(Vars)
      if (is.null(OutFolder)) {
        OutFolder <- getwd()
      }
    }

    FilesExist <- paste0(OutFolder, "/", Vars, ".RData") %>%
      file.exists()

    if (any(!Overwrite, FilesExist)) {
      "Some files already exist.\nNo files are saved. Please use overwrite = TRUE" %>%
        crayon::red() %>%
        cat(sep = "\n")
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt))
      stop()
    }

    invisible(sapply(1:ID, function(x) {
      save(
        list = paste0(Vars[x]),
        file = paste0(OutFolder, "/", Vars[x], ".RData"))
    }))

    if (all(file.exists(paste0(OutFolder, "/", Vars, ".RData")))) {
      cat(paste0(crayon::blue("All files are saved to disk", crayon::red(OutFolder), "successfully."), "\n"))
    } else {
      cat(paste0(crayon::red("Some files were not saved to disk!\nplease check again"), "\n"))
    }
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

cc <- function(...) {
  rlang::ensyms(...) %>% as.character()
}

# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
# |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

# |---------------------------------------------------| #
# LoadPackages ----
# |---------------------------------------------------| #

#' Load package silently and print version
#'
#' Load package silently and print version
#' @param Package package name
#' @name LoadPackages
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

LoadPackages <- function(Package) {
  PG <- rlang::quo_name(rlang::enquo(Package))
  suppressWarnings(suppressMessages(library(PG, character.only = TRUE)))

  Ver <- eval(parse(text = stringr::str_glue('packageVersion("{PG}")')))
  cat(stringr::str_glue(">>> {PG} - v {Ver}\n\n"))

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
#' @return NULL
#' @export

NDecimals <- function(x) {
  Split <- x %>%
    format(scientific = FALSE) %>%
    stringr::str_split(pattern = "\\.", n = Inf, simplify = TRUE)

  if (length(Split) == 2) {
    Split[, 2] %>%
      nchar() %>%
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
#' @param intern whether to make the output of the command an R object
#' @param ... additional arguments to `shell` pr `system` functions
#' @name System
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

System <- function(command, intern = TRUE, ...) {
  # Use shell() in windows OS; system() for linux OS

  # The running operating system (make also available in the global environment outside of the function)
  # CurrOS <- as.character(Sys.info()["sysname"])
  if (CurrOS() == "Windows") {
    return(shell(cmd = command, intern = intern, ...))
  }
  if (CurrOS() == "Linux") {
    return(system(command = command, intern = intern, ...))
  }
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
#' @export

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
#' @return NULL
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
#' ht(mtcars, 2)
#'
#' ht(mtcars, 6)
#' @export

ht <- function(DF, NRows = 5) {
  data.table::data.table(DF) %>%
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
#' @param DF data frame
#' @param FillVal value to be used: default: NA_character
#' @param ... list of column names to add
#' @examples
#' AddMissingCols(A, B, C, DT = mtcars, FillVal = NA_real_)
#'
#' AddCols <- c("Add1", "Add2")
#' AddMissingCols(AddCols, DT = mtcars, FillVal = NA_real_)
#' @export

AddMissingCols <- function(..., DT, FillVal = NA_character_) {
  Cols <- rlang::ensyms(...) %>%
    as.character()
  ArgInEnv <- (Cols %in% ls(envir = parent.env(environment())))

  if(any(ArgInEnv)) {
    Cols <- get(Cols, envir = parent.env(environment()))
  }

  Cols2Add <- setdiff(Cols, names(DT))
  if (length(Cols2Add) != 0) {
    DT[Cols2Add] <- FillVal
  }
  tibble::tibble(DT) %>%
    return()
}
