## |------------------------------------------------------------------------| #
# mod_cv_fit ----
## |------------------------------------------------------------------------| #

#' Prepare cross-validated Hmsc models for HPC fitting
#'
#' This function prepares cross-validated Hmsc models for fitting using HPC. It
#' handles data preparation, model initialisation, and generation of SLURM
#' commands.
#'
#' @param path_model Character. Path to a saved model file (`*.qs2`).
#' @param cv_name Character vector. Column name(s) in the model input data to be
#'   used to cross-validate the models (see [mod_prepare_data] and
#'   [mod_cv_prepare]). The function allows the possibility of using more than
#'   one way of assigning grid cells into cross-validation folders. If multiple
#'   names are provided, separate cross-validation models will be fitted for
#'   each cross-validation type. Currently, there are three cross-validation
#'   strategies: `cv_sac`, `cv_dist`, and `cv_large`. Defaults to `c("cv_dist",
#'   "cv_large")`.
#' @param partitions A vector for cross-validation created by
#'   [Hmsc::createPartition] or similar. Defaults to `NULL`, which means to use
#'   column name(s) provided in the `cv_name` argument. If the `partitions`
#'   vector is provided, the label used in the output files will be `cv_custom`.
#' @param init_par a named list of parameter values used for initialisation of
#'   MCMC states. See [Hmsc::computePredictedValues] for more information.
#'   Default: `NULL`.
#' @param job_name Character. Name of the submitted job(s) for SLURM. Default:
#'   `cv_models`.
#' @param updater named `list`. Which conditional updaters should be omitted?
#'   See [Hmsc::computePredictedValues] for more information. Defaults to
#'   `list(Gamma2 = FALSE, GammaEta = FALSE)` to disable the following warnings:
#'   `setting updater$Gamma2=FALSE due to specified phylogeny matrix` and
#'   `setting updater$GammaEta=FALSE: not implemented for spatial methods 'GPP'
#'   and 'NNGP'`.
#' @param align_posterior Logical. Whether the posterior of each chains should
#'   be aligned. See [Hmsc::computePredictedValues] for more information.
#'   Default: `TRUE`.
#' @param ... Additional arguments passed to the [IASDT.R::mod_slurm] function.
#' @details The function copies part of the [Hmsc::computePredictedValues]
#'   function, which currently does not support performing cross-validation
#'   using Hmsc-HPC.
#' @author Ahmed El-Gabbas
#' @importFrom rlang %||%
#' @inheritParams mod_slurm
#' @inheritParams mod_inputs
#' @export
#' @name mod_cv_fit

mod_cv_fit <- function(
    path_model = NULL, cv_name = c("cv_dist", "cv_large"), partitions = NULL,
    env_file = ".env", init_par = NULL, job_name = "cv_models",
    updater = list(Gamma2 = FALSE, GammaEta = FALSE),
    align_posterior = TRUE, to_json = FALSE, slurm_prepare = TRUE,
    memory_per_cpu = NULL, job_runtime = NULL, path_hmsc = NULL,
    precision = 64, ...) {

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Avoid "no visible binding for global variable" message
  # https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
  nfolds <- path_mod_init_rds <- cv <- model_name <- partition <- NULL

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Check input parameters -----
  ecokit::cat_time("Check input parameters")

  # character arguments
  ecokit::check_args(
    args_to_check = c("path_model", "job_name", "path_hmsc"),
    args_type = "character")
  # numeric arguments
  ecokit::check_args(
    args_to_check = c("to_json", "slurm_prepare"), args_type = "logical")

  if (!(precision %in% c(32, 64))) {
    ecokit::stop_ctx(
      "precision should be either of 32 or 64", precision = precision,
      include_backtrace = TRUE)
  }

  if (!file.exists(path_model)) {
    ecokit::stop_ctx(
      "Model path does not exist.", path_model = path_model,
      include_backtrace = TRUE)
  }

  if (!ecokit::check_env_file(env_file, warning = FALSE)) {
    ecokit::stop_ctx(
      "Environment file is not found or invalid.", env_file = env_file)
  }

  if (slurm_prepare) {
    # Validate memory_per_cpu
    memory_per_cpu <- .validate_slurm_ram(memory_per_cpu)
    job_runtime <- .validate_slurm_runtime(job_runtime)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Loading model -----
  ecokit::cat_time("Loading model")

  # Path of the cross-validation folder
  path_cv <- fs::path(dirname(dirname(path_model)), "model_fitting_cv")

  # Path of the model input data
  path_model_data <- fs::path(
    dirname(dirname(path_model)), "model_data_subset.RData")
  if (!ecokit::check_data(path_model_data)) {
    ecokit::stop_ctx(
      "Model data file not found",
      path_model_data = path_model_data, include_backtrace = TRUE)
  }

  # Loading full model object
  model_full <- ecokit::load_as(path_model)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Creating paths -----
  ecokit::cat_time("Creating paths")

  path_init <- fs::path(path_cv, "model_init")
  path_fitted <- fs::path(path_cv, "model_fitted")
  path_post <- fs::path(path_cv, "model_Posterior")
  path_log <- fs::path(path_post, "jobs_log")
  dir_pred <- fs::path(path_cv, "model_prediction")
  fs::dir_create(c(path_cv, path_init, path_fitted, path_post, dir_pred))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Cross-validation partitions ----
  ecokit::cat_time("Cross-validation partitions")

  if (is.null(partitions)) {
    # if custom partitions is not provided, extract the cv column(s) available
    # in the modelling data

    cv_data <- ecokit::load_as(path_model_data)$data_cv

    if (!(all(cv_name %in% names(cv_data)))) {
      # if any of the column names does not exist, stop the function
      missing_cv <- cv_name[isFALSE(cv_name %in% names(cv_data))]
      ecokit::stop_ctx(
        paste0(
          "`partitions` was not defined (NULL) and column(s) for cv folds ",
          paste(missing_cv, collapse = " + "),
          " can not be found in species data"),
        cv_name = cv_name, partitions = partitions,
        names_cv_data = names(cv_data), include_backtrace = TRUE)
    }

    # Extract cv folds from the cv column(s)
    partitions <- purrr::map(cv_name, ~ dplyr::pull(cv_data, .x))
    names(partitions) <- stringr::str_remove(cv_name, "^cv_")

  } else {
    # If partitions is provided directly to the function, use "cv_custom" as cv
    # name
    cv_name <- "cv_custom"
  }

  # Check the length of cv data equals the number of sampling units in the model
  if (any(purrr::map_int(partitions, length) != model_full$ny)) {
    ecokit::stop_ctx(
      paste0(
        "partitions parameter must be a vector of the same length of the ",
        "sampling  units of the the full model"),
      partitions = partitions, model_ny = model_full$ny,
      length_partitions = purrr::map_int(partitions, length),
      include_backtrace = TRUE)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Verbose ----

  # Verbose value should be extracted directly from the full model object. It
  # seems that models fitted with Hmsc-HPC have NULL verbose value (although not
  # tested). If the `verbose` value stored in the full model object is NULL, the
  # following assigns a value of 1000 to it (using the null coalescing operator
  # `%||%`).

  verbose <- model_full$verbose %||% 1000

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Number of chains ----
  n_chains <- length(model_full$postList)

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Prepare cross-validated initial models -----
  ecokit::cat_time("Prepare cross-validated initial models")

  # Load coordinates from the full model
  coords <- ecokit::load_as(path_model_data)$data_xy

  cv_data <- tibble::tibble(
    partition = partitions, cv_name = names(partitions)) %>%
    dplyr::mutate(

      # number of cross-validation folders for each cv type
      nfolds = purrr::map_int(.x = partition, .f = ~ length(unique(.x))),

      # prepare data for each cross-validation strategy
      cv_info = purrr::pmap(
        .l = list(partition, cv_name, nfolds),
        .f = function(partition, cv_name, nfolds) {

          ecokit::cat_time(
            paste0("Cross-validation Type: ", cv_name), level = 1L)

          cv_data0 <- purrr::map_dfr(
            .x = seq_len(nfolds),
            .f = function(k) {

              ecokit::cat_time(paste0("Fold ", k, "/", nfolds), level = 2L)

              # # |||||||||||||||||||||||
              # The following is adapted from Hmsc::computePredictedValues()
              # # |||||||||||||||||||||||

              train <- (partition != k)
              val <- (partition == k)

              switch(
                class(model_full$X)[1L],
                matrix = {
                  x_train <- model_full$X[train, , drop = FALSE]
                },
                list = {
                  x_train <- purrr::map(
                    .x = model_full$X, .f = ~ .x[train, , drop = FALSE])
                })

              if (model_full$ncRRR > 0) {
                XRRRTrain <- model_full$XRRR[train, , drop = FALSE]
                XRRRVal <- model_full$XRRR[val, , drop = FALSE]
              } else {
                XRRRTrain <- XRRRVal <- NULL
              }

              dfPi <- droplevels(model_full$dfPi[train, , drop = FALSE])

              y_train <- model_full$Y[train, , drop = FALSE]
              rownames(y_train) <- rownames(dfPi)

              # Initial models
              model_cv <- Hmsc::Hmsc(
                Y = y_train, Loff = model_full$Loff[train, , drop = FALSE],
                X = x_train, XRRR = XRRRTrain, ncRRR = model_full$ncRRR,
                XSelect = model_full$XSelect, distr = model_full$distr,
                studyDesign = dfPi, Tr = model_full$Tr, C = model_full$C,
                ranLevels = model_full$rL) %>%
                Hmsc::setPriors(
                  V0 = model_full$V0, f0 = model_full$f0,
                  mGamma = model_full$mGamma, UGamma = model_full$UGamma,
                  aSigma = model_full$aSigma, bSigma = model_full$bSigma,
                  rhopw = model_full$rhowp)

              model_cv$XFormula <- model_full$XFormula
              model_cv$XInterceptInd <- model_full$XInterceptInd

              # function to remove scaled attributes
              remove_attr <- function(x) {
                attr(x, "scaled:center") <- NULL
                attr(x, "scaled:scale") <- NULL
                x
              }

              model_cv$YScalePar <- model_full$YScalePar
              model_cv$YScaled <- scale(
                model_cv$Y, model_cv$YScalePar[1, ], model_cv$YScalePar[2, ])

              model_cv$XScalePar <- model_full$XScalePar
              model_cv$XScaled <- scale(
                model_cv$X, model_cv$XScalePar[1, ], model_cv$XScalePar[2, ])

              if (model_cv$ncRRR > 0) {
                model_cv$XRRRScalePar <- model_full$XRRRScalePar
                model_cv$XRRRScaled <- scale(
                  model_cv$XRRR,
                  model_cv$XRRRScalePar[1, ], model_cv$XRRRScalePar[2, ])
                model_cv$XRRRScaled <- remove_attr(model_cv$XRRRScaled)
              }

              model_cv$TrInterceptInd <- model_full$TrInterceptInd
              model_cv$TrScalePar <- model_full$TrScalePar
              model_cv$TrScaled <- scale(
                model_cv$Tr, model_cv$TrScalePar[1, ], model_cv$TrScalePar[2, ])
              model_cv$TrScaled <- remove_attr(model_cv$TrScaled)


              # Save unfitted model
              path_mod_init <- fs::path(
                path_init, paste0("init_mod_", cv_name, "_k", k, ".RData"))
              ecokit::save_as(
                object = model_cv,
                object_name = paste0("init_mod_", cv_name, "_k", k),
                out_path = path_mod_init)


              # initiate sampling and save initial models
              model_cv <- Hmsc::sampleMcmc(
                hM = model_cv, samples = model_full$samples,
                thin = model_full$thin, transient = model_full$transient,
                adaptNf = model_full$adaptNf, initPar = init_par,
                nChains = n_chains, updater = updater, verbose = verbose,
                alignPost = align_posterior, engine = "HPC")

              if (to_json) {
                model_cv <- jsonify::to_json(model_cv)
              }

              # Save model input as rds file
              path_mod_init_rds <- fs::path(
                path_init, paste0("init_mod_", cv_name, "_k", k, ".rds"))
              saveRDS(model_cv, file = path_mod_init_rds)

              path_mod_fitted <- fs::path(
                path_fitted, paste0("model_", cv_name, "_k", k, ".RData"))

              dfPi <- droplevels(model_full$dfPi[val, , drop = FALSE])
              val_coords <- coords[val, ]
              rownames(dfPi) <- rownames(val_coords) <- dfPi$sample

              tibble::tibble(
                path_mod_full = path_model,
                cv = k,
                path_mod_init = path_mod_init,
                path_mod_init_rds = path_mod_init_rds,
                path_mod_fitted = path_mod_fitted,
                val = list(val),
                val_coords = list(val_coords),
                XRRRVal = list(XRRRVal),
                dfPi = list(dfPi))
            })
          cv_data0
        })) %>%
    dplyr::select(-"partition") %>%
    tidyr::unnest("cv_info")

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Prepare model fitting commands -----
  ecokit::cat_time("Prepare model fitting commands")

  cv_data <- cv_data %>%
    dplyr::mutate(model_name = paste0(cv_name, cv), .after = cv) %>%
    dplyr::mutate(                # nolint: consecutive_mutate_linter
      model_prep = purrr::pmap(
        .l = list(path_mod_init_rds, cv, model_name),
        .f = function(path_mod_init_rds, cv, model_name) {

          cv_out <- purrr::map_dfr(
            .x = seq_len(n_chains),
            .f = function(chain) {

              # Turn off scientific notation
              withr::local_options(list(scipen = 999))

              # Path to save the posterior of the combination of cv and chain
              post_file <- fs::path(
                path_post,
                paste0("mod_", model_name, "_ch", chain, "_post.rds"))

              # Path to save the progress of model fitting
              mod_prog_file <- stringr::str_replace_all(
                post_file, "post.rds$", "Progress.txt")


              # Model fitting command
              command_hpc <- paste0(
                "/usr/bin/time -v python3 -m hmsc.run_gibbs_sampler", # nolint: absolute_paths_linter
                " --input ", ecokit::normalize_path(path_mod_init_rds),
                " --output ", ecokit::normalize_path(post_file),
                " --samples ", model_full$samples,
                " --transient ", model_full$transient,
                " --thin ", model_full$thin,
                " --verbose ", verbose,
                " --chain ", (chain - 1),
                " --fp ", precision,
                " >& ", ecokit::normalize_path(mod_prog_file))

              # data to be returned for each combination of cv and chain
              tibble::tibble(
                path_post = post_file, path_mod_progress = mod_prog_file,
                command_hpc = command_hpc, n_samples = model_full$samples,
                transient = model_full$transient, thin = model_full$thin)
            })

          dplyr::summarise_all(cv_out, list)

        })) %>%
    tidyr::unnest("model_prep") %>%
    dplyr::mutate_at(
      .vars = c("n_samples", "transient", "thin"),
      .funs = ~ as.integer(unique(unlist(.x))))

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save model fitting commands -----
  ecokit::cat_time("Save model fitting commands")

  purrr::walk(
    .x = cv_name,
    .f = ~{
      curr_cv_name <- stringr::str_remove_all(.x, "cv_")
      command_file <- fs::path(path_cv, paste0("commands_to_fit_", .x, ".txt"))
      f <- file(command_file, open = "wb")
      on.exit(invisible(try(close(f), silent = TRUE)), add = TRUE)

      cv_data %>%
        dplyr::filter(cv_name == curr_cv_name) %>%
        dplyr::pull(command_hpc) %>%
        unlist() %>%
        cat(sep = "\n", append = FALSE, file = f)
      close(f)
    })

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Prepare SLURM script -----

  if (slurm_prepare) {
    ecokit::cat_time("Prepare SLURM script")

    # Prepare SLURM file to submit ALL commands to HPC
    IASDT.R::mod_slurm(
      model_dir = path_cv, job_name = job_name, memory_per_cpu = memory_per_cpu,
      job_runtime = job_runtime, env_file = env_file, path_hmsc = path_hmsc,
      slurm_path_out = path_log, slurm_prefix = "cv_bash_fit", ...)
  }

  ## # ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  # Save summary data to disk -----
  ecokit::cat_time("Save summary data to disk")
  save(cv_data, file = fs::path(path_cv, "cv_data.RData"))

  return(invisible(NULL))
}
