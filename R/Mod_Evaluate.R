## |------------------------------------------------------------------------| #
# Mod_Evaluate ----
## |------------------------------------------------------------------------| #

#' Model Evaluation
#'
#' This function performs cross-validated model prediction and evaluation. It
#' merges posterior chains, and makes predictions for cross-validation folds. It
#' evaluates the model fit either for internal evaluation or on
#' cross-validation.
#'
#' @param Path_CV The path to the cross-validation folder where models and data
#'   are stored.
#' @param predictEtaMean Logical; whether to predict the mean value of the
#'   latent variable.
#' @param NCores Integer specifying the number of parallel cores for
#'   parallelization. Default: 8 cores.
#' @name Mod_Evaluate
#' @export

Mod_Evaluate <- function(Path_CV = NULL, predictEtaMean = TRUE, NCores = 8) {

  .StartTime <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  Model_Name <- Path_Pred <- Path_ModFitted <- val <- Path_Model_Full <-
    valCoords <- CurrName <- CVNames <- nfolds <- Model_Name <- NULL

  # # ..................................................................... ###

  # Merging posterior chains -----

  IASDT.R::CatTime("Merging posterior chains")

  CV_DT <- IASDT.R::LoadAs(file.path(Path_CV, "CV_DT.RData"))

  purrr::walk(
    .x = seq_len(nrow(CV_DT)),
    .f = ~ {

      ObjName <- stringr::str_remove_all(
        basename(CV_DT$Path_ModFitted[.x]), "Model_Fitted_|.RData")

      IASDT.R::CatTime(ObjName, Level = 1)

      Model_Out <- CV_DT$Path_ModFitted[.x]

      if (!file.exists(Model_Out) || !IASDT.R::CheckRData(Model_Out)) {
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

  # Making predictions - CV ----

  IASDT.R::CatTime("Making predictions - CV")

  CV_DT %>%
    dplyr::mutate(
      Pred = purrr::pmap(
        .l = list(
          Model_Name, Path_Pred, Path_Model_Full, 
          Path_ModFitted, val, valCoords),
        .f = function(Model_Name, Path_Pred, Path_Model_Full,
                      Path_ModFitted, val, valCoords) {

          IASDT.R::CatTime(Model_Name, Level = 1)
          PredsOkay <- file.exists(Path_Pred) && IASDT.R::CheckRData(Path_Pred)

          if (isFALSE(PredsOkay)) {
            Model_Full <- IASDT.R::LoadAs(Path_Model_Full)
            Model_CV <- IASDT.R::LoadAs(Path_ModFitted)
            XVal <- droplevels(Model_Full$XData[val, , drop = FALSE])

            # prepareGradient
            Gradient <- Hmsc::prepareGradient(
              hM = Model_CV, XDataNew = as.data.frame(XVal),
              sDataNew = list(sample = as.data.frame(valCoords)))

            rm(Model_Full, XVal, val, valCoords)

            # predict
            Preds <- stats::predict(
              object = Model_CV, Gradient = Gradient, nParallel = NCores,
              expected = TRUE, predictEtaMean = predictEtaMean) %>%
              simplify2array()

            IASDT.R::SaveAs(
              InObj = Preds, OutObj = paste0("Pred_", CurrName),
              OutPath = Path_Pred)

            rm(Preds)
            invisible(gc())
          }
          invisible(NULL)
        })) %>%
    invisible()

  # # ..................................................................... ###

  # Model evaluation - CV -----

  IASDT.R::CatTime("Model evaluation - CV")

  Eval_CV <- CV_DT %>%
    dplyr::select(CVNames, Path_Pred, Path_Model_Full, val, nfolds) %>%
    dplyr::group_by(CVNames) %>%
    dplyr::summarise_all(list) %>%
    dplyr::mutate(
      Path_Model_Full = purrr::map_chr(Path_Model_Full, unique),
      nfolds = purrr::map_int(nfolds, unique),
      Eval_CV = purrr::pmap(
        .l = list(CVNames, Path_Pred, Path_Model_Full, val, nfolds),
        .f = function(CVNames, Path_Pred, Path_Model_Full, val, nfolds) {

          IASDT.R::CatTime(CVNames, Level = 1)

          FullModel <- IASDT.R::LoadAs(Path_Model_Full)
          postN <- Reduce(sum, lapply(FullModel$postList, length))
          predArray <- array(NA, c(FullModel$ny, FullModel$ns, postN))

          for (Fold in seq_len(nfolds)) {
            predArray[val[[Fold]], , ] <- IASDT.R::LoadAs(Path_Pred[Fold])
            invisible(gc())
          }

          Eval_CV <- Hmsc::evaluateModelFit(FullModel, predY = predArray) %>%
            as.data.frame() %>%
            dplyr::tibble() %>%
            dplyr::mutate(
              Type = CVNames, Sp_ID = FullModel$spNames, .before = "RMSE")

          return(Eval_CV)
        })) %>%
    dplyr::pull(Eval_CV) %>%
    dplyr::bind_rows()

  # # ..................................................................... ###

  # Model evaluation - training ----

  IASDT.R::CatTime("Model evaluation - training")

  FullModel <- IASDT.R::LoadAs(unique(CV_DT$Path_Model_Full))

  preds <- Hmsc::computePredictedValues(
    hM = FullModel, nParallel = NCores, predictEtaMean = predictEtaMean,
    predictEtaMeanField = FALSE)

  Eval_Train <- Hmsc::evaluateModelFit(hM = FullModel, predY = preds) %>%
    as.data.frame() %>%
    dplyr::tibble() %>%
    dplyr::mutate(
      Type = "Training", Sp_ID = FullModel$spNames, .before = "RMSE")

  rm(preds, FullModel)
  invisible(gc())

  # # ..................................................................... ###

  Eval_DT <- dplyr::bind_rows(Eval_Train, Eval_CV)
  save(Eval_DT, file = file.path(dirname(Path_CV), "Eval_DT.RData"))

  # # ..................................................................... ###

  IASDT.R::CatDiff(InitTime = .StartTime)

  return(invisible(NULL))
}
