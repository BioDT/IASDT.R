## |------------------------------------------------------------------------| #
# Mod_MergeChains ----
## |------------------------------------------------------------------------| #

#' Post-processing of model outputs
#'
#' Post-processing of model outputs
#'
#' @param Path_Model String. Path to save all the output, including the to be fitted models (without trailing slash)
#' @param NCores Integer. Number of parallel cores to use
#' @param ModInfoName String. Default: `NULL` which means overwrite the `Model_Info.RData` file. If `ModInfoName` is provided, a new `.RData` file will be created with this prefix for file name (excluding extension)
#' @param PrintIncomplete Logical. Print to the console the name of unfitted models
#' @param FromHPC Logical. Work from HPC? This is to adjust the file paths.#' @param EnvFile String. Path to read the environment variables. Default value: `.env`
#' @param EnvFile String. Path to read the environment variables. Default value: `.env`
#' @param PlotVP Logical. Plot variance partitioning using the `IASDT.R::PlotVarPar` function.
#' @name Mod_MergeChains
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Mod_MergeChains <- function(
    Path_Model = NULL, NCores = NULL, ModInfoName = NULL, PrintIncomplete = TRUE,
    FromHPC = TRUE, EnvFile = ".env", PlotVP = TRUE) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Post_Path <- Post_Missing <- Post_Path <- M_Init_Path <- M_samples <-
    M_thin <- M_transient <- M_Name_Fit <- Path_FittedMod <- Path_Coda <-
    NMissingChains <- MissingModels <- Model_Finished <- Path_ModProg <-
    Post_Aligned <- Post_Aligned2 <- FittingTime <- NULL

  # Load .env file
  if (file.exists(EnvFile)) {
    readRenviron(EnvFile)
    Path_Scratch <- Sys.getenv("Path_LUMI_Scratch")
  } else {
    MSG <- paste0(
      "Path for environment variables: ", EnvFile, " was not found")
    stop(MSG)
  }

  # Checking arguments
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


  if (FromHPC) {
    Path_Model <- file.path(Path_Scratch, Path_Model)
  }

  # remove temp files and incomplete RDs files
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

  # Prepare working on parallel
  c1 <- snow::makeSOCKcluster(NCores)
  future::plan(future::cluster, workers = c1, gc = TRUE)
  Path_ModInfo <- file.path(Path_Model, "Model_Info.RData")

  Path_Fitted_Models <- file.path(Path_Model, "Model_Fitted")
  Path_Coda <- file.path(Path_Model, "Model_Coda")
  Path_VarPar <- file.path(Path_Model, "Model_VarPar")
  Path_Preds <- file.path(Path_Model, "Model_Preds")
  Path_Eval <- file.path(Path_Model, "Model_Eval")
  c(Path_Fitted_Models, Path_Coda, Path_VarPar, Path_Preds, Path_Eval) %>%
    fs::dir_create()


  Model_Info2 <- Path_ModInfo %>%
    IASDT.R::LoadAs() %>%
    # Check if any posterior files is missing
    dplyr::mutate(
      Post_Missing = purrr::map_lgl(
        .x = Post_Path,
        .f = ~{
          if (FromHPC) PostP <- file.path(Path_Scratch, .x)
          magrittr::not(all(file.exists(PostP)))
        }),

      # delete these columns if already exist from previous function execution
      Path_FittedMod = NULL, Path_Coda = NULL,

      # Merge posteriors and save as coda object
      ModelPosts = furrr::future_pmap(
        .l = list(
          Post_Missing, Post_Path, M_Init_Path, M_samples,
          M_thin, M_transient, M_Name_Fit, Post_Aligned),
        .f = function(Post_Missing, Post_Path, M_Init_Path, M_samples,
                      M_thin, M_transient, M_Name_Fit, Post_Aligned) {

          if (Post_Missing) {
            list(Path_FittedMod = NA, Path_Coda = NA, Post_Aligned2 = NA) %>%
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
              Posts <- purrr::map(as.character(Post_Path), IASDT.R::GetPosts)

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
                    nSamples = M_samples, thin = M_thin, transient = M_transient,
                    alignPost = FALSE))
                Post_Aligned2 <- FALSE
              } else {
                Post_Aligned2 <- TRUE
              }

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
            # Compute predicted values
            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

            Path_PredVals <- file.path(
              Path_Preds, paste0(M_Name_Fit, "_Preds.RData"))

            if (magrittr::not(file.exists(Path_PredVals))) {
              Model_Preds <- Hmsc::computePredictedValues(
                hM = Model_Fit, nParallel = length(Model_Fit$postList)) %>%
                suppressWarnings()
              IASDT.R::SaveAs(
                InObj = Model_Preds, OutObj = paste0(M_Name_Fit, "_Preds"),
                OutPath = Path_PredVals)
            } else {
              Model_Preds <- IASDT.R::LoadAs(Path_PredVals)
            }

            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            # Compute measures of model fit
            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

            Path_ModEval <- file.path(
              Path_Eval, paste0(M_Name_Fit, "_ModEval.RData"))

            if (magrittr::not(file.exists(Path_ModEval))) {
              Model_Evals <- Hmsc::evaluateModelFit(
                hM = Model_Fit, predY = Model_Preds) %>%
                suppressWarnings()
              IASDT.R::SaveAs(
                InObj = Model_Evals, OutObj = paste0(M_Name_Fit, "_Evals"),
                OutPath = Path_ModEval)
            } else {
              Model_Evals <- IASDT.R::LoadAs(Path_ModEval)
            }

            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            # Compute variance partitioning
            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

            Path_ModVarPar <- file.path(
              Path_VarPar, paste0(M_Name_Fit, "_VP.RData"))

            if (magrittr::not(file.exists(Path_ModVarPar))) {
              Model_VP <- Hmsc::computeVariancePartitioning(Model_Fit)
              IASDT.R::SaveAs(
                InObj = Model_VP, OutObj = paste0(M_Name_Fit, "_VP"),
                OutPath = Path_ModVarPar)
            } else {
              Model_VP <- IASDT.R::LoadAs(Path_ModVarPar)
            }

            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            # Plot variance partitioning
            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

            if (PlotVP) {
              Path_VP_Plot <- file.path(
                Path_VarPar, paste0(M_Name_Fit, "_VP.jpeg"))

              IASDT.R::PlotVarPar(
                Model = NULL, PlotPath = Path_VP_Plot, ModelName = NULL,
                ModelEval = Model_Evals, VarPar = Model_VP,
                EnvFile = EnvFile)
            }

            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
            # Convert to Coda object
            # # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

            Path_Coda <- file.path(
              Path_Coda, paste0(M_Name_Fit, "_Coda.RData"))
            CodaMissing <- magrittr::not(file.exists(Path_Coda))

            if (magrittr::not(ModFitMissing)) {
              Model_Fit <- IASDT.R::LoadAs(Path_FittedMod)
            }
            if (CodaMissing) {
              if (inherits(Model_Fit, "try-error")) {
                IASDT.R::CatTime("  >>>  No Coad object was exported")
              } else {
                Mod_Coda <- Hmsc::convertToCodaObject(
                  Model_Fit, spNamesNumbers = c(TRUE, FALSE),
                  covNamesNumbers = c(TRUE, FALSE))
                IASDT.R::SaveAs(
                  InObj = Mod_Coda, OutObj = paste0(M_Name_Fit, "_Coda"),
                  OutPath = Path_Coda)
                rm(Mod_Coda)
              }
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
              Post_Aligned2 = Post_Aligned2, Path_PredVals = Path_PredVals,
              Path_ModEval = Path_ModEval, Path_ModVarPar = Path_ModVarPar) %>%
              return()
          }},
        .progress = FALSE,
        .options = furrr::furrr_options(seed = TRUE, scheduling = Inf)))


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
