## |------------------------------------------------------------------------| #
# Mod_MergeChains ----
## |------------------------------------------------------------------------| #

#' Post-processing of model outputs
#'
#' This function performs post-processing of HMSC model outputs. It merges model chains, checks for missing or incomplete model runs, and optionally prints information about incomplete models. It can work with models run on High-Performance Computing (HPC) environments and supports parallel processing.

#' @param Path_Model String. Path to the directory where model outputs are stored. This should include the path to the fitted models but should not have a trailing slash.
#' @param NCores Integer. The number of cores to use for parallel processing. This should be a positive integer.
#' @param ModInfoName String. The name of the file (without extension) where the processed model information will be saved. If `NULL`, it overwrites the `Model_Info.RData` file in the `Path_Model` directory. If provided, a new `.RData` file will be created with this name.
#' @param PrintIncomplete Logical. Indicates whether to print the names of models that were not successfully fitted to the console. Defaults to `TRUE`.
#' @param FromHPC Logical. Indicates whether the function is being run in an HPC environment. This affects how file paths are handled. Defaults to `TRUE`.
#' @param EnvFile String. The path to the file containing environment variables. Defaults to ".env". This file should contain at least the `Path_LUMI_Scratch` variable when running in an HPC environment.
#' @param FromJSON Logical. Indicates whether to convert loaded models from JSON format before reading. Defaults to `FALSE`.
#' @name Mod_MergeChains
#' @author Ahmed El-Gabbas
#' @return The function does not return anything but saves the processed model information to disk.
#' @export

Mod_MergeChains <- function(
    Path_Model = NULL, NCores = NULL, ModInfoName = NULL,
    PrintIncomplete = TRUE, FromHPC = TRUE, EnvFile = ".env", FromJSON = FALSE) {

  if (is.null(Path_Model) || is.null(NCores) || is.null(ModInfoName)) {
    stop("FilePath, NCores, and ModInfoName cannot be empty")
  }

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Post_Path <- Post_Missing <- Post_Path <- M_Init_Path <- M_samples <-
    M_thin <- M_transient <- M_Name_Fit <- Path_FittedMod <- Path_Coda <-
    NMissingChains <- MissingModels <- Model_Finished <- Path_ModProg <-
    Post_Aligned <- Post_Aligned2 <- FittingTime <- NULL

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Load .env file ----
  if (file.exists(EnvFile)) {
    readRenviron(EnvFile)
    Path_Scratch <- Sys.getenv("Path_LUMI_Scratch")
  } else {
    stop(paste0("Path for environment variables: ", EnvFile, " was not found"))
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Checking arguments ----
  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = c("Path_Model", "EnvFile"), Type = "character")

  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "NCores", Type = "numeric")

  IASDT.R::CheckArgs(
    AllArgs = AllArgs, Args = c("PrintIncomplete", "FromHPC", "PlotVP"),
    Type = "logical")

  rm(AllArgs)

  if (FromHPC) {
    Path_Model <- file.path(
      Path_Scratch, stringr::str_remove(Path_Model, paste0(Path_Scratch, "/")))
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Remove temp files and incomplete RDs files ----

  Path_Model_Fit <- file.path(Path_Model, "Model_Fitting_HPC")
  tempFiles <- list.files(
    path = Path_Model_Fit, pattern = ".rds_temp$", full.names = TRUE)

  if (length(tempFiles) > 0) {
    IASDT.R::CatTime(
      paste0("There are ", length(tempFiles),
             " unsuccessful model variants to be removed"))
    tempFilesRDs <- stringr::str_replace_all(tempFiles, ".rds_temp$", ".rds")
    tempFilesProgress <- stringr::str_replace_all(
      tempFiles, "_post.rds_temp", "_Progress.txt")
    purrr::walk(
      .x = c(tempFilesRDs, tempFiles, tempFilesProgress),
      .f = ~{
        if (file.exists(.x)) file.remove(.x)
      })
  }

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Prepare working on parallel

  c1 <- snow::makeSOCKcluster(NCores)
  on.exit(invisible(try(snow::stopCluster(c1), silent = TRUE)), add = TRUE)
  future::plan(future::cluster, workers = c1, gc = TRUE)

  Path_ModInfo <- file.path(Path_Model, "Model_Info.RData")
  Path_Fitted_Models <- file.path(Path_Model, "Model_Fitted")
  Path_Coda <- file.path(Path_Model, "Model_Coda")
  fs::dir_create(c(Path_Fitted_Models, Path_Coda))


  Model_Info2 <- IASDT.R::LoadAs(Path_ModInfo) %>%
    # Check if any posterior files is missing
    dplyr::mutate(
      Post_Missing = purrr::map_lgl(
        .x = Post_Path,
        .f = ~{
          if (FromHPC) {
            Post <- file.path(Path_Scratch, stringr::str_remove_all(.x, Path_Scratch))
          } else {
            Post <- .x
          }

          file.exists(Post) %>%
            all() %>%
            magrittr::not() %>%
            return()
        }),

      # delete these columns if already exist from previous function execution
      Path_FittedMod = NULL, Path_Coda = NULL)


  Model_Info2 <- Model_Info2 %>%
    dplyr::mutate(

      # Merge posteriors and save as coda object
      ModelPosts = furrr::future_pmap(
        .l = list(
          Post_Missing, Post_Path, M_Init_Path, M_samples,
          M_thin, M_transient, M_Name_Fit, Post_Aligned),
        .f = function(Post_Missing, Post_Path, M_Init_Path, M_samples,
                      M_thin, M_transient, M_Name_Fit, Post_Aligned) {

          if (Post_Missing) {
            list(
              Path_FittedMod = NA, Path_Coda = NA, Post_Aligned2 = NA) %>%
              return()

          } else {

            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            # Merge fitted models
            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

            Path_FittedMod <- file.path(
              Path_Fitted_Models, paste0(M_Name_Fit, "_Model.RData"))
            ModFitMissing <- magrittr::not(file.exists(Path_FittedMod))

            if (ModFitMissing) {
              if (FromHPC) {
                Post_Path <- file.path(Path_Scratch, Post_Path)
                M_Init_Path <- file.path(Path_Scratch, M_Init_Path)
              }

              # Get posteriors
              Posts <- purrr::map(
                as.character(Post_Path), IASDT.R::GetPosts, FromJSON = FromJSON)

              # Convert to Hmsc object - Try with `alignPost = TRUE`
              Model_Fit <- Hmsc::importPosteriorFromHPC(
                m = IASDT.R::LoadAs(M_Init_Path), postList = Posts,
                nSamples = M_samples, thin = M_thin, transient = M_transient,
                alignPost = TRUE) %>%
                try(silent = TRUE)

              # Convert to Hmsc object - If failed, use `alignPost = FALSE`
              if (inherits(Model_Fit, "try-error")) {
                Model_Fit <- try(
                  Hmsc::importPosteriorFromHPC(
                    m = IASDT.R::LoadAs(M_Init_Path), postList = Posts,
                    nSamples = M_samples, thin = M_thin,
                    transient = M_transient, alignPost = FALSE))
                Post_Aligned2 <- FALSE
              } else {
                Post_Aligned2 <- TRUE
              }

              rm(Posts)
              invisible(gc())

              if (inherits(Model_Fit, "try-error")) {
                paste0("Model ", M_Name_Fit, " failed to be merged!") %>%
                  IASDT.R::CatTime()
              } else {
                IASDT.R::SaveAs(
                  InObj = Model_Fit, OutObj = paste0(M_Name_Fit, "_Model"),
                  OutPath = Path_FittedMod)
              }
            } else {
              Post_Aligned2 <- Post_Aligned
            }

            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            # Convert to Coda object
            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

            Path_Coda <- file.path(Path_Coda, paste0(M_Name_Fit, "_Coda.RData"))
            CodaMissing <- magrittr::not(file.exists(Path_Coda))

            if (magrittr::not(ModFitMissing)) {
              Model_Fit <- IASDT.R::LoadAs(Path_FittedMod)
            }

            if (CodaMissing) {
              if (inherits(Model_Fit, "try-error")) {
                IASDT.R::CatTime("  >>>  No Coda object was exported")
              } else {
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
            }

            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            # Remove the prefix for scratch
            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

            if (FromHPC) {
              Path_FittedMod <- stringr::str_remove(
                Path_FittedMod, paste0(Path_Scratch, "/"))
              Path_Coda <- stringr::str_remove(
                Path_Coda, paste0(Path_Scratch, "/"))
            }
            invisible(gc())

            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            # Return list of objects
            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

            list(
              Path_FittedMod = Path_FittedMod, Path_Coda = Path_Coda,
              Post_Aligned2 = Post_Aligned2) %>%
              return()

          }},
        .progress = FALSE,
        .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)))

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  Model_Info2 <- Model_Info2 %>%
    tidyr::unnest_wider("ModelPosts") %>%
    dplyr::mutate(Post_Aligned = dplyr::coalesce(Post_Aligned2)) %>%
    dplyr::select(-Post_Aligned2) %>%
    dplyr::mutate(

      # Check if both merged fitted model and coda file exist
      Model_Finished = purrr::map2_lgl(
        .x = Path_FittedMod, .y = Path_Coda,
        .f = ~{
          ModelFiles <- c(.x, .y)
          if (FromHPC) {
            ModelFiles <- file.path(Path_Scratch, ModelFiles)
          }
          all(file.exists(ModelFiles))
        }),

      # Extract fitting time from the progress file
      FittingTime = purrr::map(
        .x = Path_ModProg,
        .f = ~{
          purrr::map(
            .x = .x,
            .f = function(File) {

              if (FromHPC) {
                File <- file.path(Path_Scratch, File)
              }

              if (file.exists(File)) {
                File %>%
                  readr::read_lines() %>%
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
        }))

  snow::stopCluster(c1)

  # # |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Print to the console the name of failed models and number of missing chain files

  if (PrintIncomplete) {
    MissingModelVars <- Model_Info2 %>%
      dplyr::filter(magrittr::not(Model_Finished)) %>%
      dplyr::mutate(
        NMissingChains = purrr::map_int(
          .x = Post_Path, .f = ~sum(magrittr::not(file.exists(.x)))),
        MissingModels = paste0(M_Name_Fit, " (", NMissingChains, " chains)")
      ) %>%
      dplyr::pull(MissingModels) %>%
      IASDT.R::sort_()

    if (length(MissingModelVars) > 0) {
      IASDT.R::InfoChunk("Unsuccessful models")
      MissingModelVars %>%
        paste0("   >>  ", ., collapse = "\n") %>%
        cat()
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
  return(invisible(NULL))
}
