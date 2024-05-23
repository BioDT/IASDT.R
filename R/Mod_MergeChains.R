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
#' @name Mod_MergeChains
#' @author Ahmed El-Gabbas
#' @return NULL
#' @export

Mod_MergeChains <- function(
    Path_Model = NULL, NCores = NULL, ModInfoName = NULL) {

  # Overwrite_ModInfo = TRUE,

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Post_Path <- Post_Missing <- Post_Path <- M_Init_Path <- M_samples <-
    M_thin <- M_transient <- M_Name_Fit <- NULL

  AllArgs <- ls()
  AllArgs <- purrr::map(
    AllArgs,
    function(x) get(x, envir = parent.env(env = environment()))) %>%
    stats::setNames(AllArgs)
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "Path_Model", Type = "character")
  IASDT.R::CheckArgs(AllArgs = AllArgs, Args = "NCores", Type = "numeric")

  Path_ModInfo <- file.path(Path_Model, "Model_Info.RData")

  c1 <- snow::makeSOCKcluster(NCores)
  future::plan(future::cluster, workers = c1, gc = TRUE)

  Model_Info2 <- Path_ModInfo %>%
    IASDT.R::LoadAs() %>%
    dplyr::mutate(
      Post_Missing = purrr::map_lgl(
        .x = Post_Path,
        .f = ~magrittr::not(all(file.exists(.x)))),
      ModelPosts = furrr::future_pmap(
        .l = list(Post_Missing, Post_Path, M_Init_Path, M_samples,
                  M_thin, M_transient, M_Name_Fit),
        .f = function(Post_Missing, Post_Path, M_Init_Path, M_samples,
                      M_thin, M_transient, M_Name_Fit) {

          if (Post_Missing) {
            return(list(Path_FittedMod = NULL, Path_Coda = NULL, Post_Aligned = NULL))
          } else {

            Path_FittedMod <- file.path(
              Path_Model, "ModelFitting", paste0(M_Name_Fit, "_Fitted.RData"))
            ModFitMissing <- magrittr::not(file.exists(Path_FittedMod))

            Path_Coda <- file.path(
              Path_Model, "ModelFitting", paste0(M_Name_Fit, "_Coda.RData"))
            CodaMissing <- magrittr::not(file.exists(Path_Coda))

            if (ModFitMissing) {
              Posts <- purrr::map(as.character(Post_Path), IASDT.R::GetPosts)
              Model_Fit <- try(Hmsc::importPosteriorFromHPC(
                m = IASDT.R::LoadAs(M_Init_Path), postList = Posts,
                nSamples = M_samples, thin = M_thin, transient = M_transient,
                alignPost = TRUE)) %>%
                try(silent = TRUE)

              if (inherits(Model_Fit, "try-error")) {
                Model_Fit <- try(Hmsc::importPosteriorFromHPC(
                  m = IASDT.R::LoadAs(M_Init_Path), postList = Posts,
                  nSamples = M_samples, thin = M_thin, transient = M_transient,
                  alignPost = FALSE))
                Post_Aligned <- FALSE
              } else {
                Post_Aligned <- TRUE
              }

              IASDT.R::SaveAs(
                InObj = Model_Fit, OutObj = M_Name_Fit, OutPath = Path_FittedMod)
            } else {
              Post_Aligned <- NULL
            }

            if (CodaMissing) {
              if (magrittr::not(ModFitMissing)) {
                Model_Fit <- IASDT.R::LoadAs(Path_FittedMod)
              }

              Mod_Coda <- Hmsc::convertToCodaObject(
                Model_Fit, spNamesNumbers = c(TRUE, FALSE),
                covNamesNumbers = c(TRUE, FALSE))

              IASDT.R::SaveAs(
                InObj = Mod_Coda, OutObj = paste0(M_Name_Fit, "_coda"),
                OutPath = Path_Coda)
              rm(Mod_Coda)
            }

            invisible(gc())
            list(Path_FittedMod = Path_FittedMod,
                 Path_Coda = Path_Coda,
                 Post_Aligned = Post_Aligned) %>%
              return()
          }},
        .progress = FALSE,
        .options = furrr::furrr_options(seed = TRUE, scheduling = Inf))) %>%
    tidyr::unnest_wider("ModelPosts")

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
