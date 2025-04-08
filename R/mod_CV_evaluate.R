## |------------------------------------------------------------------------| #
# mod_CV_evaluate ----
## |------------------------------------------------------------------------| #

#' Evaluate cross-validated Hmsc Model
#'
#' This function evaluates Hmsc models trained via cross-validation using. It
#' merges posterior chains, makes predictions, and computes 4 evaluation
#' metrics: AUC, RMSE, Tjur's R2, and Boyce Index.
#'
#' @param path_CV Character. The directory path where cross-validation models
#'   and outputs are stored.
#' @param predict_eta_mean boolean flag indicating whether to use the estimated
#'   mean values of posterior predictive distribution for random effects
#'   corresponding for the new units. See `Hmsc:::predict.Hmsc` for more
#'   details.
#' @param n_cores Integer. Number of CPU cores to use for parallel processing.
#'   Default: 8.
#' @author Ahmed El-Gabbas
#' @export
#' @name mod_CV_evaluate

mod_CV_evaluate <- function(
  path_CV = NULL, predict_eta_mean = TRUE, n_cores = 8L) {

  # # ..................................................................... ###

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  ModName <- Path_Pred <- Path_ModFitted <- val <- Path_ModFull <- IAS_ID <-
    valCoords <- CV_name <- ModName <- Path_Eval <- Path_Pred_Eval_all <- NULL

  # # ..................................................................... ###

  # Merging posterior chains -----

  IASDT.R::cat_time("Merging posterior chains")

  CV_DT <- IASDT.R::load_as(IASDT.R::path(path_CV, "CV_DT.RData"))
  PredDir <- IASDT.R::path(path_CV, "Model_Prediction")
  EvalDir <- IASDT.R::path(path_CV, "Evaluation")
  TempDir <- IASDT.R::path(PredDir, "Preds_Temp")
  fs::dir_create(c(PredDir, EvalDir, TempDir))

  purrr::walk(
    .x = seq_len(nrow(CV_DT)),
    .f = ~ {
      ObjName <- stringr::str_remove_all(
        basename(CV_DT$Path_ModFitted[.x]), "Model_Fitted_|.RData|.qs2")

      IASDT.R::cat_time(ObjName, level = 1)
      Model_Out <- CV_DT$Path_ModFitted[.x]

      if (!IASDT.R::check_data(Model_Out, warning = FALSE)) {
        Init <- IASDT.R::load_as(CV_DT$Path_ModInit[.x])
        Posts <- purrr::map(
          .x = CV_DT$Path_Post[[.x]],
          .f = function(Post) {
            readRDS(Post)$list[[1]]
          })

        Model <- Hmsc::importPosteriorFromHPC(
          m = Init, postList = Posts, nSamples = CV_DT$NSamples[.x],
          thin = CV_DT$Thin[.x], transient = TRUE)

        IASDT.R::save_as(
          object = Model, object_name = ObjName, out_path = Model_Out)
      }
    })

  # # ..................................................................... ###

  # Making predictions ----

  IASDT.R::cat_time("Making predictions")

  CV_DT <- CV_DT %>%
    dplyr::mutate(
      Path_Pred_Eval = purrr::pmap(
        .l = list(ModName, Path_ModFull, Path_ModFitted, val, valCoords),
        .f = function(ModName, Path_ModFull, Path_ModFitted, val, valCoords) {
          IASDT.R::cat_time(ModName, level = 1)

          Path_Eval <- IASDT.R::path(
            EvalDir, paste0("Eval_", ModName, ".RData"))
          Path_Pred <- IASDT.R::path(PredDir, paste0("Preds_", ModName, ".qs2"))

          PredsOkay <- IASDT.R::check_data(Path_Eval, warning = FALSE)

          if (isFALSE(PredsOkay)) {
            Model_Full <- IASDT.R::load_as(Path_ModFull)
            Model_CV <- IASDT.R::load_as(Path_ModFitted)
            XVal <- droplevels(Model_Full$XData[val, , drop = FALSE])
            # pred_PA <- Model_Full$Y[val, ]
            # pred_XY <- Model_Full$rL$sample$s[val, ]

            # prepareGradient
            # Gradient <- Hmsc::prepareGradient(
            #   hM = Model_CV, XDataNew = as.data.frame(XVal),
            #   sDataNew = list(sample = as.data.frame(valCoords)))

            rm(Model_Full, XVal, val, valCoords, envir = environment())
            invisible(gc())

            # predict / evaluate
            # Eval <- stats::predict(
            #   object = Model_CV, Gradient = Gradient, nParallel = n_cores,
            #   expected = TRUE, predictEtaMean = predict_eta_mean,
            #   TempDir = TempDir, PredDir = PredDir, evaluate = TRUE,
            #   Evaluate_Name = paste0("Eval_", ModName), EvalDir = EvalDir,
            #   ModelName = ModName, RC = FALSE, pred_PA = pred_PA,
            #   pred_XY = pred_XY)
          }

          return(
            tibble::tibble(Path_Pred = Path_Pred, Path_Eval = Path_Eval))
        })) %>%
    tidyr::unnest_wider("Path_Pred_Eval")

  CV_Eval_Preds <- CV_DT %>%
    dplyr::select(CV_name, Path_Pred, Path_Eval) %>%
    dplyr::summarise(
      .by = CV_name,
      Path_Pred = list(Path_Pred), Path_Eval = list(Path_Eval)) %>%
    dplyr::mutate(
      Path_Pred_Eval_all = purrr::pmap(
        .l = list(CV_name, Path_Pred, Path_Eval),
        .f = function(CV_name, Path_Pred, Path_Eval) {
          Path_Preds <- IASDT.R::path(
            PredDir, paste0("Preds_", CV_name, "_all.RData"))
          Preds <- purrr::map_dfr(.x = Path_Pred, .f = IASDT.R::load_as) %>%
            IASDT.R::save_as(
              object = Preds, object_name = paste0("Preds_", CV_name, "_all"),
              out_path = Path_Preds)

          Path_Evals <- IASDT.R::path(
            EvalDir, paste0("Eval_", CV_name, "_all.RData"))
          Evals <- purrr::map(
            .x = Path_Eval,
            .f = ~ {
              CV <- readr::parse_number(stringr::str_extract(.x, "CV..RData"))
              IASDT.R::load_as(.x) %>%
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
          IASDT.R::save_as(
            object = Evals, object_name = paste0("Eval_", CV_name, "_all"),
            out_path = Path_Evals)

          return(
            tibble::tibble(Path_Preds = Path_Preds, Path_Evals = Path_Evals))
        })) %>%
    tidyr::unnest_wider(Path_Pred_Eval_all)

  save(CV_DT, file = IASDT.R::path(path_CV, "CV_DT.RData"))
  save(CV_Eval_Preds, file = IASDT.R::path(path_CV, "CV_Eval_Preds.RData"))

  # # ..................................................................... ###

  IASDT.R::cat_diff(init_time = .StartTime)

  return(invisible(NULL))
}
