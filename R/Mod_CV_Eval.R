## |------------------------------------------------------------------------| #
# Mod_CV_Eval ----
## |------------------------------------------------------------------------| #

#' Evaluate cross-validated Hmsc Model
#'
#' This function evaluates Hmsc models trained via cross-validation using. It
#' merges posterior chains, makes predictions, and computes 4 evaluation
#' metrics: AUC, RMSE, Tjur's R2, and Boyce Index.
#'
#' @param Path_CV Character. The directory path where cross-validation models
#'   and outputs are stored.
#' @param predictEtaMean boolean flag indicating whether to use the estimated
#'   mean values of posterior predictive distribution for random effects
#'   corresponding for the new units. See `Hmsc:::predict.Hmsc` for more
#'   details.
#' @param NCores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @author Ahmed El-Gabbas
#' @export
#' @name Mod_CV_Eval


Mod_CV_Eval <- function(Path_CV = NULL, predictEtaMean = TRUE, NCores = 8L) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ModName <- Path_Pred <- Path_ModFitted <- val <- Path_ModFull <- IAS_ID <-
    valCoords <- CVName <- ModName <- Path_Eval <- Path_Pred_Eval_all <- NULL

  # # ..................................................................... ###

  # Merging posterior chains -----

  IASDT.R::CatTime("Merging posterior chains")

  CV_DT <- IASDT.R::LoadAs(IASDT.R::Path(Path_CV, "CV_DT.RData"))
  PredDir <- IASDT.R::Path(Path_CV, "Model_Prediction")
  EvalDir <- IASDT.R::Path(Path_CV, "Evaluation")
  TempDir <- IASDT.R::Path(PredDir, "Preds_Temp")
  fs::dir_create(c(PredDir, EvalDir, TempDir))

  purrr::walk(
    .x = seq_len(nrow(CV_DT)),
    .f = ~ {
      ObjName <- stringr::str_remove_all(
        basename(CV_DT$Path_ModFitted[.x]), "Model_Fitted_|.RData|.qs2")

      IASDT.R::CatTime(ObjName, Level = 1)
      Model_Out <- CV_DT$Path_ModFitted[.x]

      if (!IASDT.R::CheckData(Model_Out, warning = FALSE)) {
        Init <- IASDT.R::LoadAs(CV_DT$Path_ModInit[.x])
        Posts <- purrr::map(
          .x = CV_DT$Path_Post[[.x]],
          .f = function(Post) {
            readRDS(Post)$list[[1]]
          })

        Model <- Hmsc::importPosteriorFromHPC(
          m = Init, postList = Posts, nSamples = CV_DT$NSamples[.x],
          thin = CV_DT$Thin[.x], transient = TRUE)

        IASDT.R::SaveAs(InObj = Model, OutObj = ObjName, OutPath = Model_Out)
      }
    })

  # # ..................................................................... ###

  # Making predictions ----

  IASDT.R::CatTime("Making predictions")

  CV_DT <- CV_DT %>%
    dplyr::mutate(
      Path_Pred_Eval = purrr::pmap(
        .l = list(ModName, Path_ModFull, Path_ModFitted, val, valCoords),
        .f = function(ModName, Path_ModFull, Path_ModFitted, val, valCoords) {
          IASDT.R::CatTime(ModName, Level = 1)

          Path_Eval <- IASDT.R::Path(
            EvalDir, paste0("Eval_", ModName, ".RData"))
          Path_Pred <- IASDT.R::Path(PredDir, paste0("Preds_", ModName, ".qs2"))

          PredsOkay <- IASDT.R::CheckData(Path_Eval, warning = FALSE)

          if (isFALSE(PredsOkay)) {
            Model_Full <- IASDT.R::LoadAs(Path_ModFull)
            Model_CV <- IASDT.R::LoadAs(Path_ModFitted)
            XVal <- droplevels(Model_Full$XData[val, , drop = FALSE])
            # Pred_PA <- Model_Full$Y[val, ]
            # Pred_XY <- Model_Full$rL$sample$s[val, ]

            # prepareGradient
            # Gradient <- Hmsc::prepareGradient(
            #   hM = Model_CV, XDataNew = as.data.frame(XVal),
            #   sDataNew = list(sample = as.data.frame(valCoords)))

            rm(Model_Full, XVal, val, valCoords, envir = environment())
            invisible(gc())

            # predict / evaluate
            # Eval <- stats::predict(
            #   object = Model_CV, Gradient = Gradient, nParallel = NCores,
            #   expected = TRUE, predictEtaMean = predictEtaMean,
            #   TempDir = TempDir, PredDir = PredDir, Evaluate = TRUE,
            #   Evaluate_Name = paste0("Eval_", ModName), EvalDir = EvalDir,
            #   ModelName = ModName, RC = FALSE, Pred_PA = Pred_PA,
            #   Pred_XY = Pred_XY)
          }

          return(
            tibble::tibble(Path_Pred = Path_Pred, Path_Eval = Path_Eval))
        })) %>%
    tidyr::unnest_wider("Path_Pred_Eval")

  CV_Eval_Preds <- CV_DT %>%
    dplyr::select(CVName, Path_Pred, Path_Eval) %>%
    dplyr::summarise(
      .by = CVName,
      Path_Pred = list(Path_Pred), Path_Eval = list(Path_Eval)) %>%
    dplyr::mutate(
      Path_Pred_Eval_all = purrr::pmap(
        .l = list(CVName, Path_Pred, Path_Eval),
        .f = function(CVName, Path_Pred, Path_Eval) {
          Path_Preds <- IASDT.R::Path(
            PredDir, paste0("Preds_", CVName, "_all.RData"))
          Preds <- purrr::map_dfr(.x = Path_Pred, .f = IASDT.R::LoadAs) %>%
            IASDT.R::SaveAs(
              InObj = Preds, OutObj = paste0("Preds_", CVName, "_all"),
              OutPath = Path_Preds)

          Path_Evals <- IASDT.R::Path(
            EvalDir, paste0("Eval_", CVName, "_all.RData"))
          Evals <- purrr::map(
            .x = Path_Eval,
            .f = ~ {
              CV <- readr::parse_number(stringr::str_extract(.x, "CV..RData"))
              IASDT.R::LoadAs(.x) %>%
                dplyr::rename_with(
                  .fn = ~ paste0(.x, "_CV", CV),
                  .cols = c("RMSE", "AUC", "Boyce", "TjurR2"))
            }) %>%
            purrr::reduce(.f = dplyr::left_join, by = c("Sp", "IAS_ID")) %>%
            dplyr::select(
              IAS_ID,
              tidyselect::starts_with("AUC_"),
              tidyselect::starts_with("Boyce_"),
              tidyselect::starts_with("TjurR2_"),
              tidyselect::starts_with("RMSE_"))
          IASDT.R::SaveAs(
            InObj = Evals, OutObj = paste0("Eval_", CVName, "_all"),
            OutPath = Path_Evals)

          return(
            tibble::tibble(Path_Preds = Path_Preds, Path_Evals = Path_Evals))
        })) %>%
    tidyr::unnest_wider(Path_Pred_Eval_all)

  save(CV_DT, file = IASDT.R::Path(Path_CV, "CV_DT.RData"))
  save(CV_Eval_Preds, file = IASDT.R::Path(Path_CV, "CV_Eval_Preds.RData"))

  # # ..................................................................... ###

  IASDT.R::CatDiff(InitTime = .StartTime)

  return(invisible(NULL))
}
