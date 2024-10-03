## |------------------------------------------------------------------------| #
# Mod_MergeChains ----
## |------------------------------------------------------------------------| #

#' Post-processing of model outputs
#'
#' This function performs post-processing of HMSC model outputs. It merges model
#' chains, checks for missing or incomplete model runs, and optionally prints
#' information about incomplete models. It can work with models run on HPC
#' environments and supports parallel processing.
#' @param Path_Model String. Path to the root directory of the fitted models
#'   without the trailing slash. Two folders will be created `Model_Fitted` and
#'   `Model_Coda` to store merged model and coda objects, respectively.
#' @param NCores Integer. The number of cores to use for parallel processing.
#'   This should be a positive integer.
#' @param ModInfoName String. The name of the file (without extension) where the
#'   processed model information will be saved. If `NULL`, it overwrites the
#'   `Model_Info.RData` file in the `Path_Model` directory. If provided, a new
#'   `.RData` file will be created with this name.
#' @param PrintIncomplete Logical. Indicates whether to print the names of
#'   models that were not successfully fitted to the console. Defaults to
#'   `TRUE`.
#' @param FromHPC Logical. Indicates whether the function is being run in an HPC
#'   environment. This affects how file paths are handled. Defaults to `TRUE`.
#' @param FromJSON Logical. Indicates whether to convert loaded models from JSON
#'   format before reading. Defaults to `FALSE`.
#' @name Mod_MergeChains
#' @author Ahmed El-Gabbas
#' @return The function does not return anything but saves the processed model
#'   information to disk.
#' @export

Mod_MergeChains <- function(
    Path_Model = NULL, NCores = NULL, ModInfoName = NULL,
    PrintIncomplete = TRUE, FromHPC = TRUE, FromJSON = FALSE) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Post_Path <- Post_Path <- M_Name_Fit <- Path_FittedMod <- Path_Coda <-
    NMissingChains <- MissingModels <- Model_Finished <- Path_ModProg <-
    Post_Aligned2 <- FittingTime <- FittingMemory <- NULL

  # # ..................................................................... ###

  # Checking arguments ----

  if (is.null(Path_Model) || is.null(NCores)) {
    stop("Path_Model, and NCores cannot be empty", call. = FALSE)
  }

  AllArgs <- ls(envir = environment())
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = "Path_Model", Type = "character")

  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "NCores", Type = "numeric")

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = c("PrintIncomplete", "FromHPC", "PlotVP"),
    Type = "logical")

  rm(AllArgs)

  if (!dir.exists(Path_Model)) {
    stop(paste0(
      "Path_Model directory (`", Path_Model, "`) does not exist"),
      call. = FALSE)
  }

  # # ..................................................................... ###

  # Creating paths -----

  Path_Fitted_Models <- file.path(Path_Model, "Model_Fitted")
  Path_Coda <- file.path(Path_Model, "Model_Coda")
  fs::dir_create(c(Path_Fitted_Models, Path_Coda))

  # # ..................................................................... ###

  # Loading model info ----

  Path_ModInfo <- file.path(Path_Model, "Model_Info.RData")

  if (!file.exists(Path_ModInfo)) {
    stop(paste0(
      "ModInfo file `", Path_ModInfo, "` does not exist"), call. = FALSE)
  }

  # # ..................................................................... ###

  # Remove temp files and incomplete RDs files ----

  Path_Model_Fit <- file.path(Path_Model, "Model_Fitting_HPC")
  tempFiles <- list.files(
    path = Path_Model_Fit, pattern = ".rds_temp$", full.names = TRUE)

  if (length(tempFiles) > 0) {
    IASDT.R::CatTime(
      paste0(
        "There are ", length(tempFiles),
        " unsuccessful model variants to be removed"))

    tempFilesRDs <- stringr::str_replace_all(tempFiles, ".rds_temp$", ".rds")
    tempFilesProgress <- stringr::str_replace_all(
      tempFiles, "_post.rds_temp", "_Progress.txt")

    purrr::walk(
      .x = c(tempFilesRDs, tempFiles, tempFilesProgress),
      .f = ~{
        if (file.exists(.x)) {
          file.remove(.x)
        }
      })
  }

  # # ..................................................................... ###

  # Prepare working on parallel -----
  withr::local_options(future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)

  if (NCores == 1) {
    future::plan("sequential", gc = TRUE)
  } else {
    c1 <- snow::makeSOCKcluster(NCores)
    on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
    future::plan("cluster", workers = c1, gc = TRUE)
    on.exit(future::plan("sequential", gc = TRUE), add = TRUE)
  }

  Model_Info2 <- IASDT.R::LoadAs(Path_ModInfo) %>%
    dplyr::mutate(
      # Check if any posterior files is missing
      Post_Missing = purrr::map_lgl(
        .x = Post_Path, .f = ~ all(!file.exists(.x))),
      # delete these columns if already exist from previous function execution
      Path_FittedMod = NULL, Path_Coda = NULL)

  # # ..................................................................... ###

  # Merge posteriors and save as Hmsc model / coda object

  Model_Info3 <- future.apply::future_lapply(
    X = seq_len(nrow(Model_Info2)),
    FUN = function(x) {

      if (Model_Info2$Post_Missing[[x]]) {
        return(list(
          Path_FittedMod = NA_character_,
          Path_Coda = NA_character_,
          Post_Aligned2 = NA)
        )
      }

      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # Merge fitted models
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      M_Name_Fit <- Model_Info2$M_Name_Fit[[x]]

      Path_FittedMod <- file.path(
        Path_Fitted_Models, paste0(M_Name_Fit, "_Model.RData"))

      # Check if model exists and RData is valid
      ModelFileOkay <- file.exists(Path_FittedMod) &&
        IASDT.R::CheckRData(Path_FittedMod)

      if (isFALSE(ModelFileOkay)) {

        # delete corrupted file
        if (file.exists(Path_FittedMod)) {
          fs::file_delete(Path_FittedMod)
        }

        # Get posteriors
        Posts <- purrr::map(
          .x = as.character(Model_Info2$Post_Path[[x]]),
          .f = IASDT.R::GetPosts, FromJSON = FromJSON)

        # Convert to Hmsc object
        # Try with `alignPost = TRUE`
        Model_Fit <- Hmsc::importPosteriorFromHPC(
          m = IASDT.R::LoadAs(Model_Info2$M_Init_Path[x]),
          postList = Posts, nSamples = Model_Info2$M_samples[x],
          thin = Model_Info2$M_thin[x],
          transient = Model_Info2$M_transient[x],
          alignPost = TRUE) %>%
          try(silent = TRUE)

        # If failed, use `alignPost = FALSE`
        if (inherits(Model_Fit, "try-error")) {
          Model_Fit <- try(
            Hmsc::importPosteriorFromHPC(
              m = IASDT.R::LoadAs(Model_Info2$M_Init_Path[x]),
              postList = Posts, nSamples = Model_Info2$M_samples[x],
              thin = Model_Info2$M_thin[x],
              transient = Model_Info2$M_transient[x], alignPost = FALSE),
            silent = TRUE)
          Post_Aligned2 <- FALSE
        } else {
          Post_Aligned2 <- TRUE
        }

        rm(Posts)
        invisible(gc())

        if (inherits(Model_Fit, "try-error")) {
          return(list(
            Path_FittedMod = NA_character_,
            Path_Coda = NA_character_,
            Post_Aligned2 = NA))
        } else {
          IASDT.R::SaveAs(
            InObj = Model_Fit, OutObj = paste0(M_Name_Fit, "_Model"),
            OutPath = Path_FittedMod)
        }
      } else {
        Post_Aligned2 <- Model_Info2$Post_Aligned[[x]]
      }

      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
      # Convert to Coda object
      # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

      Path_Coda <- file.path(Path_Coda, paste0(M_Name_Fit, "_Coda.RData"))

      if (!exists("Model_Fit") && !file.exists(Path_FittedMod)) {
        list(
          Path_FittedMod = NA_character_,
          Path_Coda = NA_character_,
          Post_Aligned2 = NA) %>%
          return()
      }

      if (!exists("Model_Fit")) {
        Model_Fit <- IASDT.R::LoadAs(Path_FittedMod)
      }

      # Check if coda file exists and RData is valid
      CodaFileOkay <- file.exists(Path_Coda) && IASDT.R::CheckRData(Path_Coda)

      if (isFALSE(CodaFileOkay)) {

        if (file.exists(Path_Coda)) {
          fs::file_delete(Path_Coda)
        }

        Mod_Coda <- Hmsc::convertToCodaObject(
          Model_Fit, spNamesNumbers = c(TRUE, FALSE),
          covNamesNumbers = c(TRUE, FALSE))

        IASDT.R::SaveAs(
          InObj = Mod_Coda, OutObj = paste0(M_Name_Fit, "_Coda"),
          OutPath = Path_Coda)

        rm(Mod_Coda)
      }

      rm(Model_Fit)
      invisible(gc())

      # Return list of objects
      list(
        Path_FittedMod = Path_FittedMod, Path_Coda = Path_Coda,
        Post_Aligned2 = Post_Aligned2) %>%
        return()
    },
    future.scheduling = Inf, future.seed = TRUE,
    future.globals = c(
      "Model_Info2", "Path_Fitted_Models", "FromJSON", "Path_Coda"),
    future.packages = c("Hmsc", "coda", "purrr", "IASDT.R", "dplyr"))

  if (NCores > 1) {
    snow::stopCluster(c1)
    future::plan("sequential", gc = TRUE)
  }

  # # ..................................................................... ###

  Model_Info2 <- dplyr::mutate(Model_Info2, ModelPosts = Model_Info3) %>%
    tidyr::unnest_wider("ModelPosts") %>%
    dplyr::mutate(Post_Aligned = dplyr::coalesce(Post_Aligned2)) %>%
    dplyr::select(-Post_Aligned2) %>%
    dplyr::mutate(

      # Check if both merged fitted model and coda file exist
      Model_Finished = purrr::map2_lgl(
        .x = Path_FittedMod, .y = Path_Coda,
        .f = ~all(file.exists(c(.x, .y)))),

      # Extract fitting time from the progress file
      FittingTime = purrr::map(
        .x = Path_ModProg,
        .f = ~{
          purrr::map(
            .x = .x,
            .f = function(File) {
              if (file.exists(File)) {
                readr::read_lines(file = File) %>%
                  stringr::str_subset("Whole Gibbs sampler elapsed") %>%
                  stringr::str_remove("Whole Gibbs sampler elapsed") %>%
                  stringr::str_trim() %>%
                  as.numeric() %>%
                  magrittr::divide_by(60) %>%
                  round(1)
              } else {
                NA_real_
              }
            }) %>%
            unlist()
        }),

      # Mean elapsed time
      FittingTimeMean = purrr::map2_dbl(
        .x =  Model_Finished, .y = FittingTime,
        .f = ~{
          if (.x) {
            mean(.y)
          } else {
            NA
          }
        }),

      # Maximum memory
      FittingMemory = purrr::map(
        .x = Path_ModProg,
        .f = ~{
          purrr::map(
            .x = .x,
            .f = function(File) {
              if (file.exists(File)) {
                readr::read_lines(file = File) %>%
                  stringr::str_subset("Maximum resident set size") %>%
                  stringr::str_remove_all(
                    "\t|:|Maximum resident set size \\(kbytes\\)") %>%
                  stringr::str_trim() %>%
                  as.numeric() %>%
                  magrittr::divide_by(1024 * 1024) %>%
                  round(2)
              } else {
                NA_real_
              }
            }) %>%
            unlist()
        }),


      # Mean memory
      FittingMemoryMean = purrr::map2_dbl(
        .x =  Model_Finished, .y = FittingMemory,
        .f = ~{
          if (.x) {
            mean(.y)
          } else {
            NA
          }
        }))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Print to the console the name of failed models and number of missing chain
  # files

  if (PrintIncomplete) {
    MissingModelVars <- Model_Info2 %>%
      dplyr::filter(!Model_Finished) %>%
      dplyr::mutate(
        NMissingChains = purrr::map_int(
          .x = Post_Path, .f = ~sum(!file.exists(.x))),
        MissingModels = paste0(M_Name_Fit, " (", NMissingChains, " chains)")
      ) %>%
      dplyr::pull(MissingModels) %>%
      IASDT.R::sort_()

    if (length(MissingModelVars) > 0) {
      IASDT.R::InfoChunk("Unsuccessful models", Extra2 = 0)
      purrr::walk(MissingModelVars, IASDT.R::CatTime, Level = 1)
    }
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save Model_Info to disk

  if (is.null(ModInfoName)) {
    IASDT.R::SaveAs(
      InObj = Model_Info2, OutObj = "Model_Info", OutPath = Path_ModInfo)
  } else {
    IASDT.R::SaveAs(
      InObj = Model_Info2, OutObj = ModInfoName,
      OutPath = file.path(Path_Model, paste0(ModInfoName, ".RData")))
  }

  IASDT.R::CatDiff(InitTime = .StartTime)

  return(invisible(NULL))
}
