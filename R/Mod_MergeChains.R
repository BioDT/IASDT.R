## |------------------------------------------------------------------------| #
# Mod_MergeChains ----
## |------------------------------------------------------------------------| #

#' Prepare habitat-specific data for the models
#'
#' Prepare habitat-specific data for the models
#'
#' @param Path_Model String. Path to save all the output, including the to be fitted models (without trailing slash)
#' @param NCores Integer. Number of parallel cores to use
#' @param ModInfoName String. Default: `NULL` which means overwrite the `Model_Info.RData` file. If `ModInfoName` is provided, a new `.RData` file will be created with this prefix for file name (excluding extension)
#' @param PrintIncomplete Logical. Print to the console the name of unfitted models
#' @name Mod_MergeChains
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Mod_MergeChains <- function(
    Path_Model = NULL, NCores = NULL, ModInfoName = NULL,
    PrintIncomplete = TRUE) {

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Post_Path <- Post_Missing <- Post_Path <- M_Init_Path <- M_samples <-
    M_thin <- M_transient <- M_Name_Fit <- Path_FittedMod <- Path_Coda <-
    NMissingChains <- MissingModels <- Model_Finished <- Path_ModPorg <-
    Post_Aligned <- Post_Aligned2 <- NULL

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "Path_Model", Type = "character")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "NCores", Type = "numeric")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "PrintIncomplete", Type = "logical")

  Path_ModInfo <- file.path(Path_Model, "Model_Info.RData")

  c1 <- snow::makeSOCKcluster(NCores)
  future::plan(future::cluster, workers = c1, gc = TRUE)


  Model_Info2 <- Path_ModInfo %>%
    IASDT.R::LoadAs() %>%
    dplyr::mutate(
      Post_Missing = purrr::map_lgl(
        .x = Post_Path,
        .f = ~magrittr::not(all(file.exists(.x)))),

      # delete these columns if already exist from previous function execution
      Path_FittedMod = NULL, Path_Coda = NULL,

      ModelPosts = furrr::future_pmap(
        .l = list(Post_Missing, Post_Path, M_Init_Path, M_samples,
                  M_thin, M_transient, M_Name_Fit, Post_Aligned),
        .f = function(Post_Missing, Post_Path, M_Init_Path, M_samples,
                      M_thin, M_transient, M_Name_Fit, Post_Aligned) {

          if (Post_Missing) {
            list(Path_FittedMod = NA, Path_Coda = NA, Post_Aligned2 = NA) %>%
              return()
          } else {

            Path_Fitted_Models <- file.path(Path_Model, "ModelFitted")
            fs::dir_create(Path_Fitted_Models)

            Path_FittedMod <- file.path(
              Path_Fitted_Models, paste0(M_Name_Fit, "_Fitted.RData"))
            Path_Coda <- file.path(
              Path_Fitted_Models, paste0(M_Name_Fit, "_Coda.RData"))

            ModFitMissing <- magrittr::not(file.exists(Path_FittedMod))
            CodaMissing <- magrittr::not(file.exists(Path_Coda))

            if (ModFitMissing) {
              Posts <- purrr::map(as.character(Post_Path), IASDT.R::GetPosts)
              Model_Fit <- try(
                Hmsc::importPosteriorFromHPC(
                  m = IASDT.R::LoadAs(M_Init_Path), postList = Posts,
                  nSamples = M_samples, thin = M_thin, transient = M_transient,
                  alignPost = TRUE)) %>%
                try(silent = TRUE)

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
                  InObj = Model_Fit, OutObj = M_Name_Fit,
                  OutPath = Path_FittedMod)
              }

            } else {
              Post_Aligned2 <- Post_Aligned
            }

            if (CodaMissing) {
              if (magrittr::not(ModFitMissing)) {
                Model_Fit <- IASDT.R::LoadAs(Path_FittedMod)
              }

              if (inherits(Model_Fit, "try-error")) {
                IASDT.R::CatTime("  >>>  No Coad object was exported")
              } else {
                Mod_Coda <- Hmsc::convertToCodaObject(
                  Model_Fit, spNamesNumbers = c(TRUE, FALSE),
                  covNamesNumbers = c(TRUE, FALSE))
                IASDT.R::SaveAs(
                  InObj = Mod_Coda, OutObj = paste0(M_Name_Fit, "_coda"),
                  OutPath = Path_Coda)
                rm(Mod_Coda)
              }
            }

            invisible(gc())
            list(
              Path_FittedMod = Path_FittedMod, Path_Coda = Path_Coda,
              Post_Aligned2 = Post_Aligned2) %>%
              return()
          }},
        .progress = FALSE,
        .options = furrr::furrr_options(seed = TRUE, scheduling = Inf))) %>%
    tidyr::unnest_wider("ModelPosts") %>%
    dplyr::mutate(Post_Aligned = dplyr::coalesce(Post_Aligned2)) %>%
    dplyr::select(-Post_Aligned2) %>%
    dplyr::mutate(
      Model_Finished = purrr::map2_lgl(
        .x = Path_FittedMod, .y = Path_Coda,
        .f = ~all(file.exists(c(.x, .y)))),
      FittingTime = purrr::map(
        .x = Path_ModPorg,
        .f = ~{
          purrr::map(
            .x = .x,
            .f = function(File) {
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
        }))

  snow::stopCluster(c1)

  if (PrintIncomplete) {
    IASDT.R::InfoChunk("Unsuccessful models")
    Model_Info2 %>%
      dplyr::filter(magrittr::not(Model_Finished)) %>%
      dplyr::mutate(
        NMissingChains = purrr::map_int(
          .x = Post_Path, .f = ~sum(magrittr::not(file.exists(.x)))),
        MissingModels = paste0(M_Name_Fit, " (", NMissingChains, " chains)")
      ) %>%
      dplyr::pull(MissingModels) %>%
      IASDT.R::sort_() %>%
      paste0("   >>  ", ., collapse = "\n") %>%
      cat()
  }

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
