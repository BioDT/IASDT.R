## |------------------------------------------------------------------------| #
# variance_partitioning_compute ----
## |------------------------------------------------------------------------| #

#' Computes and visualise variance partitioning of Hmsc models
#'
#' Computes and plots variance components with respect to given grouping of
#' fixed effects and levels of random effects. The The
#' `variance_partitioning_compute()` function inherits the main functionality
#' from the [Hmsc::computeVariancePartitioning] function, but with the added
#' functionality of parallel computation and using `TensorFlow`. The
#' `variance_partitioning_plot()` function generates plots for variance
#' partitioning as JPEG files. It allows for sorting the predictors and species;
#' e.g., by the mean value per predictor; and by original species order. It also
#' plots the raw variance partitioning (relative variance partitioning
#' multiplied by the training and testing (if supported) Tjur-R<sup>2</sup>
#' value).
#' @param path_model Character. Path to fitted `Hmsc` model object.
#' @param group vector of numeric values corresponding to group identifiers in
#'   groupnames. If the model was defined with `XData` and `XFormula`, the
#'   default is to use model terms.
#' @param group_names vector of names for each group of fixed effect. Should
#'   match `group`. If the model was defined with `XData` and `XFormula`, the
#'   default is to use the labels of model terms.
#' @param start index of first MCMC sample included. Default: `1L`.
#' @param na.ignore Logical. If `TRUE`, covariates are ignored for sites where
#'   the focal species is NA when computing variance-covariance matrices for
#'   each species.
#' @param n_cores Integer. Number of CPU cores to use for computing variance
#'   partitioning using `TensorFlow`. This is only effective when `use_tf` is
#'   `TRUE`. Default: `1`.
#' @param chunk_size Integer. Size of each chunk of samples to process in
#'   parallel. Only relevant for `TensorFlow`. Default: `50`.
#' @param verbose Logical. Whether to print progress messages. Default: `TRUE`.
#' @param temp_cleanup Logical. Whether to delete temporary files after
#'   processing. Default: `TRUE`.
#' @param vp_commands_only Logical. If `TRUE`, returns the commands to run the
#'   Python script. Default is `FALSE`. Only relevant when `use_tf` is `TRUE`.
#' @param vp_file Character. Name of the output file to save the results.
#'   Default: `varpar`.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param width,height Numeric. Width and height of the output plot in
#'   centimetres. Default: `30` and `15`, respectively.
#' @param is_cv_model Logical. Whether the model is a cross-validated model
#'   (`TRUE`) or fitted with the full dataset (`FALSE`; default). If `TRUE`, the
#'   explanatory and predictive power of the model will be used to estimate the
#'   raw variance partitioning.
#' @param temp_dir Character. Path to a temporary directory to store
#'   intermediate files. Default: `NULL`, which creates a temporary directory in
#'   the same parent directory as the model file.
#' @param axis_text Numeric. Size of the axis text. Default: `4`.
#' @author Ahmed El-Gabbas
#' @importFrom foreach %dopar%
#' @inheritParams predict_hmsc
#' @inheritParams predict_latent_factor
#' @name variance_partitioning
#' @rdname variance_partitioning
#' @order 1
#' @export

variance_partitioning_compute <- function(
    path_model, group = NULL, group_names = NULL, start = 1L, na.ignore = FALSE,
    n_cores = 8L, use_tf = TRUE, tf_environ = NULL, tf_use_single = FALSE,
    temp_cleanup = TRUE, chunk_size = 50L, verbose = TRUE, vp_file = "varpar",
    vp_commands_only = FALSE, temp_dir = NULL) {

  x <- NULL

  n_cores <- .validate_n_cores(n_cores)

  # # .................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # set up parallel processing
  doParallel::registerDoParallel(cores = n_cores)
  ecokit::load_packages(package_list = "foreach")
  withr::defer(doParallel::stopImplicitCluster())


  # # .................................................................... ###

  # Create folder for variance partitioning results
  path_var_par <- fs::path(
    dirname(dirname(path_model)),
    "model_postprocessing", "variance_partitioning")
  fs::dir_create(path_var_par)

  file_var_par <- fs::path(path_var_par, paste0(vp_file, ".RData"))
  if (ecokit::check_data(file_var_par, warning = FALSE)) {
    return(ecokit::load_as(file_var_par))
  }

  # # .................................................................... ###

  # Check if the virtual environment and Python scripts exist

  if (use_tf) {

    # Check python virtual environment
    if (isFALSE(vp_commands_only) && .Platform$OS.type == "windows" &&
        (is.null(tf_environ) || !dir.exists(tf_environ))) {
      ecokit::stop_ctx(
        paste0(
          "When running on Windows and `use_tf` is TRUE, `tf_environ` must ",
          "be specified and point to an existing directory with a ",
          "Python virtual environment"),
        vp_commands_only = vp_commands_only, OS = .Platform$OS.type,
        tf_environ = tf_environ, include_backtrace = TRUE)
    }

    # Determine the Python executable path

    # On Windows, the TF calculations has to be done through a valid virtual
    # environment; the path to the virtual environment must be specified in
    # `tf_environ`. On LUMI, this is not needed as the compatible python
    # installation is loaded automatically when loading `tensorflow` module.
    # When using another HPC system, the function needs to be adapted
    # accordingly to point to a valid python virtual environment.

    if (.Platform$OS.type == "windows") {
      python_executable <- fs::path(tf_environ, "Scripts", "python.exe")

      if (isFALSE(vp_commands_only) && !file.exists(python_executable)) {
        ecokit::stop_ctx(
          "Python executable not found in the virtual environment",
          python_executable = python_executable, include_backtrace = TRUE)
      }
    } else {
      python_executable <- "/usr/bin/time -v python3" # nolint: absolute_paths_linter
    }

    # Check GPU availability
    if (isFALSE(vp_commands_only) && .Platform$OS.type == "windows") {
      result <- system(
        paste0(
          python_executable,
          " -c \"import tensorflow as tf; print(len(",
          "tf.config.list_physical_devices('GPU')))\""),
        intern = TRUE)
      N_GPU <- result[length(result)]
      if (N_GPU == 0) {
        ecokit::cat_time(
          "No GPU found; Calculations will use CPU.",
          cat_timestamp = FALSE, cat_bold = TRUE, cat_red = TRUE,
          verbose = verbose)
      } else {
        ecokit::cat_time(
          paste0(N_GPU, " GPUs were found. Calculations will use GPU."),
          cat_timestamp = FALSE, cat_bold = TRUE, cat_red = TRUE,
          verbose = verbose)
      }
    }

    # Paths to the Python scripts
    script_geta <- system.file("VP_geta.py", package = "IASDT.R")
    script_getf <- system.file("VP_getf.py", package = "IASDT.R")
    script_gemu <- system.file("VP_gemu.py", package = "IASDT.R")
    if (!all(file.exists(c(script_geta, script_getf, script_gemu)))) {
      ecokit::stop_ctx(
        "Necessary python scripts do not exist",
        script_geta = script_geta, script_getf = script_getf,
        script_gemu = script_gemu, include_backtrace = TRUE)
    }
  }

  # # .................................................................... ###

  # Load model object ------
  ecokit::cat_time("Load model object", verbose = verbose)

  if (is.null(path_model) || !file.exists(path_model)) {
    ecokit::stop_ctx(
      "Model path is NULL or does not exist", path_model = path_model,
      include_backtrace = TRUE)
  }

  model_obj <- ecokit::load_as(path_model)

  # # .................................................................... ###

  ny <- model_obj$ny
  nc <- model_obj$nc
  ns <- model_obj$ns
  nr <- model_obj$nr

  if (is.null(group)) {
    # names of variables used in the model
    group_names <- names(model_obj$XData)

    # actual variables used in the model, including quadratic terms
    model_vars <- dimnames(model_obj$X)[[2]][-1]

    # group variable to combine variable and its quadratic terms together
    group <- purrr::map(
      model_vars, ~ which(stringr::str_detect(.x, group_names))) %>%
      unlist() %>%
      # add intercept to the first group
      c(.[1], .)
  }

  # If na.ignore=T, convert XData to a list
  if (na.ignore) {
    xl <- list()
    for (s in seq_len(ns)) {
      xl[[s]] <- model_obj$X
    }
    model_obj$X <- xl
  }

  switch(
    class(model_obj$X)[1L],
    matrix = {
      cMA <- stats::cov(model_obj$X)
    },
    list = {
      if (na.ignore) {
        cMA <- list()
        for (s in seq_len(ns)) {
          cMA[[s]] <- stats::cov(
            model_obj$X[[s]][which(model_obj$Y[, s] > -Inf), ])
        }
      } else {
        cMA <- lapply(model_obj$X, stats::cov)
      }
    }
  )

  # # .................................................................... ###

  # Prepare postList-----

  ecokit::cat_time("Prepare postList", verbose = verbose)
  postList <- Hmsc::poolMcmcChains(model_obj$postList, start = start)

  # Remove not-needed items

  Items2Delete <- c(
    "Eta", "Psi", "V", "sigma", "Delta", "Alpha",
    "rho", "wRRR", "PsiRRR", "DeltaRRR")
  postList <- purrr::map(
    .x = postList,
    .f = ~ {
      .x[Items2Delete] <- NULL
      return(.x)
    })
  invisible(gc())

  # Remove unnecessary elements from the model object

  Items2Delete <- c(
    "postList", "Y", "XScaled", "rL", "ranLevels", "XData", "dfPi",
    "studyDesign", "C", "Pi", "phyloTree", "XFormula", "XScalePar",
    "aSigma", "bSigma", "TrScaled", "YScalePar", "call", "rhopw",
    "distr", "V0", "UGamma", "YScaled")

  # This takes long time in some cases (when submitted via SLURM)
  # model_obj[names_to_remove] <- NULL

  # `trim_hmsc` is much faster
  model_obj <- IASDT.R::trim_hmsc(
    model = model_obj, names_to_remove = Items2Delete)
  invisible(gc())

  pool_n <- length(postList)
  ngroups <- max(group)

  # # .................................................................... ###
  # # .................................................................... ###

  # Prepare `la`/`lf`/`lmu` lists -----

  if (use_tf) {

    ecokit::cat_time(
      "Prepare/check VP files for `TensorFlow`", verbose = verbose)

    # Create the temporary directory
    if (is.null(temp_dir)) {
      temp_dir <- fs::path(dirname(dirname(path_model)), "temp_vp")
    }
    fs::dir_create(temp_dir)
    path_vp_input_files <- fs::path(temp_dir, "vp_input_files.txt")

    file_suffix <- stringr::str_pad(
      string = seq_len(pool_n), pad = "0", width = 4)

    # List of feather files resulted from `geta` function
    Files_la <- fs::path(temp_dir, paste0("vp_a_", file_suffix, ".feather"))
    files_la_exist <- all(fs::file_exists(Files_la))

    if (files_la_exist) {
      files_la_exist <- foreach::foreach(
        i = Files_la, .combine = c, .multicombine = TRUE,
        .packages = c("ecokit", "arrow")) %dopar% { # nolint: object_usage_linter
          ecokit::check_data(i, warning = FALSE)
        }
      files_la_exist <- all(files_la_exist)
    }
    if (files_la_exist) {
      ecokit::cat_time(
        "All `vp_a*.feather` files are available",
        verbose = verbose, level = 1L)
    }

    # List of feather files resulted from `getf` function
    Files_lf <- fs::path(temp_dir, paste0("vp_f_", file_suffix, ".feather"))
    files_lf_exist <- all(fs::file_exists(Files_lf))
    if (files_lf_exist) {
      files_lf_exist <- foreach::foreach(
        i = Files_lf, .combine = c, .multicombine = TRUE,
        .packages = c("ecokit", "arrow")) %dopar% { # nolint: object_usage_linter
          ecokit::check_data(i, warning = FALSE)
        }
      files_lf_exist <- all(files_lf_exist)
    }
    if (files_lf_exist) {
      ecokit::cat_time(
        "All `vp_f*.feather` files are available",
        verbose = verbose, level = 1L)
    }

    # List of feather files resulted from `gemu` function
    Files_lmu <- fs::path(temp_dir, paste0("vp_mu_", file_suffix, ".feather"))
    files_lmu_exist <- all(fs::file_exists(Files_lmu))
    if (files_lmu_exist) {
      files_lmu_exist <- foreach::foreach(
        i = Files_lmu, .combine = c, .multicombine = TRUE,
        .packages = c("ecokit", "arrow")) %dopar% { # nolint: object_usage_linter
          ecokit::check_data(i, warning = FALSE)
        }
      files_lmu_exist <- all(files_lmu_exist)
    }
    if (files_lmu_exist) {
      ecokit::cat_time(
        "All `vp_mu*.feather` files are available",
        verbose = verbose, level = 1L)
    }

    beta_files <- fs::path(
      temp_dir, paste0("vp_beta_", file_suffix, ".feather"))

    if (all(c(files_la_exist, files_lf_exist, files_lmu_exist))) {

      # All feather data are already processed and available on disk
      ecokit::cat_time(
        "Data for `la`/`lf`/`lmu` lists were already processed",
        verbose = verbose, level = 1L)

    } else {

      # Write the contents of Files_la, Files_lf, Files_lmu to a text file
      readr::write_lines(
        x = c(Files_la, Files_lf, Files_lmu), file = path_vp_input_files)

      ## Prepare la/lf/lmu lists using `TensorFlow` ----
      ecokit::cat_time(
        "Prepare la/lf/lmu lists using `TensorFlow`", verbose = verbose)

      ### Prepare data for `TensorFlow` ----
      ecokit::cat_time(
        "Prepare data for `TensorFlow`", level = 1L, verbose = verbose)

      #### X data -----
      # needed only to calculate `geta` and `getf` functions
      if (!all(c(files_la_exist, files_lf_exist))) {
        ecokit::cat_time("X", level = 2L, verbose = verbose)
        path_x <- fs::path(temp_dir, "vp_x.feather")
        if (!file.exists(path_x)) {
          arrow::write_feather(as.data.frame(model_obj$X), path_x)
        }
      }

      #### Tr / Gamma ------
      # needed only to calculate `geta` and `gemu` functions
      if (!all(c(files_la_exist, files_lmu_exist))) {
        ecokit::cat_time("Tr", level = 2L, verbose = verbose)
        path_tr <- fs::path(temp_dir, "vp_tr.feather")
        if (!file.exists(path_tr)) {
          arrow::write_feather(as.data.frame(model_obj$Tr), path_tr)
        }

        # Gamma - convert each list item into a column in a data frame
        ecokit::cat_time("Gamma", level = 2L, verbose = verbose)
        path_gamma <- fs::path(temp_dir, "vp_gamma.feather")
        if (!file.exists(path_gamma)) {
          gamma_data <- postList %>%
            purrr::map(~as.vector(.x[["Gamma"]])) %>%
            as.data.frame() %>%
            stats::setNames(paste0("Sample_", seq_len(ncol(.))))
          arrow::write_feather(gamma_data, path_gamma)
        }
      }

      #### Beta -----
      # only needed for `getf` function
      if (!files_lf_exist) {

        # Beta -- Each element of Beta is a matrix, so each list item is saved
        # to separate feather file
        ecokit::cat_time("Beta", level = 2L, verbose = verbose)

        beta_files_Exist <- all(fs::file_exists(beta_files))
        if (beta_files_Exist) {
          beta_files_Exist <- foreach::foreach(
            i = beta_files, .combine = c, .multicombine = TRUE,
            .packages = c("ecokit", "arrow")) %dopar% { # nolint: object_usage_linter
              ecokit::check_data(i, warning = FALSE)
            }
          beta_files_Exist <- all(beta_files_Exist)
        }

        if (!beta_files_Exist) {

          ecokit::cat_time(
            "Processing beta in parallel", level = 3L, verbose = verbose)

          beta_0 <- foreach::foreach(
            x = seq_along(postList), .combine = c, .multicombine = TRUE,
            .packages = c("ecokit", "fs", "purrr", "arrow", "magrittr"),
            .export = c("beta_files", "postList")) %dopar% { # nolint: object_usage_linter

              beta_file <- beta_files[x]

              if (ecokit::check_data(beta_file, warning = FALSE)) {
                return(NULL)
              }

              if (fs::file_exists(beta_file)) {
                try(fs::file_delete(beta_file), silent = TRUE)
              }

              Beta <- purrr::pluck(postList, x, "Beta") %>%
                as.data.frame()

              attempt <- 1
              repeat {

                arrow::write_feather(x = Beta, sink = beta_file)
                Sys.sleep(1)

                if (ecokit::check_data(beta_file, warning = FALSE)) {
                  break
                }

                if (attempt >= 5) {
                  ecokit::stop_ctx(
                    "Failed to create Beta file after multiple attempts",
                    beta_file = beta_file, attempt = attempt,
                    include_backtrace = TRUE)
                }

                attempt <- attempt + 1
              }
              return(NULL)
            }
          rm(beta_0, envir = environment())
        }
      }

      invisible(gc())

      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

      ### Processing geta -----

      if (files_la_exist) {

        ecokit::cat_time(
          "All `la` data were already available on disk",
          level = 1L, verbose = verbose)

      } else {

        ecokit::cat_time(
          "Processing `geta` function", level = 1L, verbose = verbose)
        path_out_a <- fs::path(temp_dir, "vp_a.feather") %>%
          ecokit::normalize_path()

        cmd_a <- paste(
          python_executable, script_geta,
          "--tr", ecokit::normalize_path(path_tr, must_work = TRUE),
          "--x", ecokit::normalize_path(path_x, must_work = TRUE),
          "--gamma", ecokit::normalize_path(path_gamma, must_work = TRUE),
          "--output", path_out_a,
          "--chunk_size", chunk_size)

        if (tf_use_single) {
          cmd_a <- paste0(cmd_a, " --use_single")
        }

        if (vp_commands_only) {

          # Save command to file
          # Redirect results of time to a log file
          path_log_a <- stringr::str_replace(path_out_a, ".feather", ".log")
          cmd_a <- paste0(cmd_a, paste0(" >> ", path_log_a, " 2>&1"))
          readr::write_lines(
            x = cmd_a,
            file = fs::path(temp_dir, "vp_a_command.txt"), append = FALSE)

        } else {

          # Run the command using system
          la <- system(cmd_a, wait = TRUE, intern =  TRUE)

          # Check for errors
          if (inherits(la, "error") || la[length(la)] != "Done") {
            ecokit::stop_ctx(
              "Error in computing geta", la = la, class_la = class(la),
              include_backtrace = TRUE)
          }

          if (length(la) != 1) {
            cat(la, sep = "\n")
          }
        }
      }

      invisible(gc())

      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

      ### Processing getf ----

      if (files_la_exist) {

        ecokit::cat_time(
          "All `lf` data were already available on disk",
          level = 1L, verbose = verbose)

      } else {

        ecokit::cat_time(
          "Processing `getf` function", level = 1L, verbose = verbose)
        path_out_f <- fs::path(temp_dir, "vp_f.feather") %>%
          ecokit::normalize_path()

        cmd_f <- paste(
          python_executable, script_getf,
          "--x", ecokit::normalize_path(path_x, must_work = TRUE),
          "--beta_dir", ecokit::normalize_path(temp_dir, must_work = TRUE),
          "--output", path_out_f,
          "--ncores", 1)

        if (tf_use_single) {
          cmd_f <- paste0(cmd_f, " --use_single")
        }

        if (vp_commands_only) {

          # Save command to file
          # Redirect results of time to a log file
          path_log_f <- stringr::str_replace(path_out_f, ".feather", ".log")
          cmd_f <- paste0(cmd_f, paste0(" >> ", path_log_f, " 2>&1"))
          readr::write_lines(
            x = cmd_f,
            file = fs::path(temp_dir, "vp_f_command.txt"),
            append = FALSE)

        } else {

          # Run the command using system
          lf <- system(cmd_f, wait = TRUE, intern =  TRUE)

          # Check for errors
          if (inherits(lf, "error") || lf[length(lf)] != "Done") {
            ecokit::stop_ctx(
              "Error in computing geta", lf = lf, class_lf = class(lf),
              include_backtrace = TRUE)
          }

          if (length(lf) != 1) {
            cat(lf, sep = "\n")
          }
        }
      }

      invisible(gc())

      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

      ### Processing gemu ----

      if (files_lmu_exist) {

        ecokit::cat_time(
          "All `lmu` data were already available on disk",
          level = 1L, verbose = verbose)

      } else {

        ecokit::cat_time(
          "Processing `gemu` function", level = 1L, verbose = verbose)
        path_out_mu <- fs::path(temp_dir, "vp_mu.feather") %>%
          ecokit::normalize_path()

        cmd_mu <- paste(
          python_executable, script_gemu,
          "--tr", ecokit::normalize_path(path_tr, must_work = TRUE),
          "--gamma", ecokit::normalize_path(path_gamma, must_work = TRUE),
          "--output", path_out_mu,
          "--ncores", 1,
          "--chunk_size", chunk_size)

        if (tf_use_single) {
          cmd_mu <- paste0(cmd_mu, " --use_single")
        }

        if (vp_commands_only) {

          # Save command to file
          # Redirect results of time to a log file
          path_log_mu <- stringr::str_replace(path_out_mu, ".feather", ".log")
          cmd_mu <- paste0(cmd_mu, paste0(" >> ", path_log_mu, " 2>&1"))
          readr::write_lines(
            x = cmd_mu,
            file = fs::path(temp_dir, "vp_mu_command.txt"),
            append = FALSE)

        } else {

          # Run the command using system
          lmu <- system(cmd_mu, wait = TRUE, intern =  TRUE)

          # Check for errors
          if (inherits(lmu, "error") || lmu[length(lmu)] != "Done") {
            ecokit::stop_ctx(
              "Error in computing geta", lmu = lmu, class_lmu = class(lmu),
              include_backtrace = TRUE)
          }

          if (length(lmu) != 1) {
            cat(lmu, sep = "\n")
          }
        }
      }
    }

    invisible(gc())

    if (vp_commands_only) {
      return(invisible(NULL))
    }

  } else {

    # Prepare la/lf/lmu lists using original R code -----

    ecokit::cat_time(
      "Prepare la/lf/lmu lists using original R code", verbose = verbose)

    ## geta -----
    ecokit::cat_time("Running geta", level = 1L, verbose = verbose)
    geta <- function(a) {
      switch(
        class(model_obj$X)[1L],
        matrix = {
          res <- model_obj$X %*% (t(model_obj$Tr %*% t(a$Gamma)))
        },
        list = {
          res <- matrix(NA, ny, ns)
          for (j in seq_len(ns)) {
            res[, j] <- model_obj$X[[j]] %*%
              (t(model_obj$Tr[j, ] %*% t(a$Gamma)))
          }
        })
      res
    }

    la <- lapply(postList, geta)

    invisible(gc())

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

    ## getf ------

    ecokit::cat_time("Running getf", level = 1L, verbose = verbose)
    getf <- function(a) {
      switch(
        class(model_obj$X)[1L],
        matrix = {
          res <- model_obj$X %*% (a$Beta)
        },
        list = {
          res <- matrix(NA, ny, ns)
          for (j in seq_len(ns)) res[, j] <- model_obj$X[[j]] %*% a$Beta[, j]
        })
      res
    }

    lf <- lapply(postList, getf)

    invisible(gc())

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

    ## gemu ------

    ecokit::cat_time("Running gemu", level = 1L, verbose = verbose)
    gemu <- function(a) {
      t(model_obj$Tr %*% t(a$Gamma))
    }

    lmu <- lapply(postList, gemu)

    invisible(gc())

  }

  # # .................................................................... ###

  # Running gebeta -------
  ecokit::cat_time("Running gebeta", verbose = verbose)
  gebeta <- function(a) {
    a$Beta
  }
  lbeta <- lapply(postList, gebeta)

  # # .................................................................... ###

  # Remove Gamma from postList ------
  ecokit::cat_time("Remove Gamma from postList", verbose = verbose)
  postList <- purrr::map(
    postList,
    ~ {
      .x["Gamma"] <- NULL
      .x
    })

  invisible(gc())

  # # .................................................................... ###
  # # .................................................................... ###

  # Computing variance partitioning -------

  if (use_tf) {

    # compute variance partitioning in parallel if `TensorFlow` was used

    ecokit::cat_time(
      "Computing variance partitioning in parallel in R", verbose = verbose)

    ecokit::cat_time(
      "Split `lbeta` list into small qs2 files", level = 1L, verbose = verbose)
    path_lbeta <- fs::path(temp_dir, "lbeta")
    purrr::walk(
      .x = seq_along(lbeta),
      .f = ~ {
        ecokit::save_as(
          object = lbeta[[.x]],
          out_path = fs::path(
            path_lbeta,
            paste0(
              "lbeta_", stringr::str_pad(.x, width = 4, pad = "0"), ".qs2")))
      })

    ecokit::cat_time(
      "Split `postList` list into small qs2 files",
      level = 1L, verbose = verbose)
    path_postList <- fs::path(temp_dir, "postList")
    purrr::walk(
      .x = seq_along(postList),
      .f = ~{
        ecokit::save_as(
          object = postList[[.x]],
          out_path = fs::path(
            path_postList,
            paste0(
              "post_", stringr::str_pad(.x, width = 4, pad = "0"), ".qs2")))
      })

    ecokit::cat_time(
      "removing `postList` and `lbeta` list objects",
      level = 1L, verbose = verbose)
    n_postList <- length(postList)
    rm(postList, lbeta, envir = environment())

    model_x <- model_obj$X
    model_obj <- IASDT.R::trim_hmsc(
      model = model_obj, names_to_remove = c("X", "covNames"))
    invisible(gc())

    ecokit::cat_time("Processing in parallel", level = 1L, verbose = verbose)

    vars_to_export <- c(
      "ngroups", "Files_la", "Files_lf", "Files_lmu", "path_lbeta", "nc",
      "model_x", "path_postList", "ns", "nr", "cMA", "group")
    packages_to_load <- c(
      "Matrix", "dplyr", "arrow", "purrr", "qs2", "methods",
      "stringr", "fs", "ecokit", "magrittr")

    Res <- foreach::foreach(
      i = seq_len(n_postList), .export = vars_to_export,
      .packages = packages_to_load) %dopar% { # nolint: object_usage_linter

        mm <- methods::getMethod("%*%", "Matrix")

        curr_postList <- fs::path(
          path_postList,
          paste0(
            "post_", stringr::str_pad(i, width = 4, pad = "0"), ".qs2")) %>%
          ecokit::load_as()
        Beta <- curr_postList$Beta
        Lambdas <- curr_postList$Lambda
        rm(curr_postList, envir = environment())
        invisible(gc())

        # Suppress warnings when no trait information is used in the models
        # cor(Beta[k, ], lmu[k, ]) : the standard deviation is zero
        data_lmu <- as.matrix(arrow::read_feather(Files_lmu[i]))      # nolint: object_name_linter
        curr_lbeta <- fs::path(      # nolint: object_name_linter
          path_lbeta,
          paste0(
            "lbeta_", stringr::str_pad(i, width = 4, pad = "0"), ".qs2")) %>%
          ecokit::load_as()
        r2t_beta <- purrr::map_dbl(
          .x = seq_len(nc),
          .f = ~ {
            suppressWarnings(stats::cor(curr_lbeta[.x, ], data_lmu[.x, ])^2)
          })

        rm(curr_lbeta, data_lmu, envir = environment())
        invisible(gc())

        fixed1 <- matrix(0, nrow = ns, ncol = 1)
        fixedsplit1 <- matrix(0, nrow = ns, ncol = ngroups)
        random1 <- matrix(0, nrow = ns, ncol = nr)

        data_la <- as.matrix(arrow::read_feather(Files_la[i]))
        data_lf <- as.matrix(arrow::read_feather(Files_lf[i]))
        # a <- data_la - matrix(rep(rowMeans(data_la), ns), ncol = ns)
        # f <- data_lf - matrix(rep(rowMeans(data_lf), ns), ncol = ns)
        a <- data_la - Matrix::rowMeans(data_la)
        f <- data_lf - Matrix::rowMeans(data_lf)

        res1 <- sum((Matrix::rowSums((a * f)) / (ns - 1))^2)
        res2 <- sum((Matrix::rowSums((a * a)) / (ns - 1)) *
                      (rowSums((f * f)) / (ns - 1)))
        r2t_y <- res1 / res2

        for (j in seq_len(ns)) {
          switch(
            class(model_x)[1L],
            matrix = {
              cM <- cMA
            },
            list = {
              cM <- cMA[[j]]
            })

          ftotal <- Matrix::crossprod(Beta[, j], mm(cM, Beta[, j]))
          fixed1[j] <- fixed1[j] + ftotal

          for (k in seq_len(ngroups)) {
            sel <- (group == k)
            fpart <- Matrix::crossprod(
              Beta[sel, j], mm(cM[sel, sel], Beta[sel, j]))
            fixedsplit1[j, k] <- fixedsplit1[j, k] + fpart
          }
        }
        rm(ftotal, fpart, envir = environment())

        for (level in seq_len(nr)) {
          Lambda <- Lambdas[[level]]
          nf <- dim(Lambda)[[1]]
          for (factor in seq_len(nf)) {
            random1[, level] <- random1[, level] +
              t(Lambda[factor, ]) * Lambda[factor, ]
          }
        }
        rm(Lambda, envir = environment())

        if (nr > 0) {
          tot <- fixed1 + Matrix::rowSums(random1)
          fixed <- fixed1 / tot
          for (level in seq_len(nr)) {
            random <- random1[, level] / tot
          }
        } else {
          fixed <- matrix(1, nrow = ns, ncol = 1)
          random <- random1
        }

        # fixedsplit <- matrix(0, nrow = ns, ncol = ngroups)
        # for (k in seq_len(ngroups)) {
        #   fixedsplit[, k] <- fixedsplit1[, k] / Matrix::rowSums(fixedsplit1)
        # }

        fixedsplit <- fixedsplit1 / Matrix::rowSums(fixedsplit1)

        rm(
          fixedsplit1, fixed1, random1, a, f, res1, res2, tot,
          envir = environment())

        return(
          list(
            fixed = fixed, random = random, fixedsplit = fixedsplit,
            r2t_y = r2t_y, r2t_beta = r2t_beta))
      }

    # Check if foreach returned a flattened list
    if (length(Res) == 5 * n_postList) {
      # Reshape into list of n_postList elements, each with 5 items
      Res <- lapply(
        X = seq_len(n_postList),
        FUN = function(i) {
          start_idx <- (i - 1) * 5 + 1
          list(
            fixed = Res[[start_idx]],
            random = Res[[start_idx + 1]],
            fixedsplit = Res[[start_idx + 2]],
            r2t_y = Res[[start_idx + 3]],
            r2t_beta = Res[[start_idx + 4]]
          )
        })
    }

    # Summarise the results
    ecokit::cat_time("Summarise the results", level = 1L, verbose = verbose)
    fixed <- Reduce("+", purrr::map(Res, ~ .x$fixed)) / pool_n
    random <- Reduce("+", purrr::map(Res, ~ .x$random)) / pool_n
    fixedsplit <- Reduce("+", purrr::map(Res, ~ .x$fixedsplit)) / pool_n
    r2t_y <- Reduce("+", purrr::map(Res, ~ .x$r2t_y)) / pool_n
    r2t_beta <- Reduce("+", purrr::map(Res, ~ .x$r2t_beta)) / pool_n

    rm(Res, model_x, envir = environment())
    invisible(gc())

  } else {

    ecokit::cat_time(
      "Computing variance partitioning sequentially", verbose = verbose)

    mm <- methods::getMethod("%*%", "Matrix")

    fixed <- matrix(0, nrow = ns, ncol = 1)
    fixedsplit <- matrix(0, nrow = ns, ncol = ngroups)
    random <- matrix(0, nrow = ns, ncol = nr)
    r2t_y <- 0
    r2t_beta <- rep(0, nc)

    for (i in seq_len(pool_n)) {

      if (i %% 200 == 0) {
        ecokit::cat_time(
          sprintf(
            "Processing iteration %d of %d", i, pool_n),
          verbose = verbose)
      }

      data_la <- la[[i]]
      data_lf <- lf[[i]]
      data_lmu <- lmu[[i]]

      # Suppress warnings when no trait information is used in the models
      # cor(Beta[k, ], lmu[k, ]) : the standard deviation is zero
      suppressWarnings({
        for (k in seq_len(nc)) {
          r2t_beta[k] <- r2t_beta[k] +
            stats::cor(lbeta[[i]][k, ], data_lmu[k, ])^2
        }
      })

      fixed1 <- matrix(0, nrow = ns, ncol = 1)
      fixedsplit1 <- matrix(0, nrow = ns, ncol = ngroups)
      random1 <- matrix(0, nrow = ns, ncol = nr)
      Beta <- postList[[i]]$Beta
      Lambdas <- postList[[i]]$Lambda

      a <- data_la - matrix(rep(rowMeans(data_la), ns), ncol = ns)
      f <- data_lf - matrix(rep(rowMeans(data_lf), ns), ncol = ns)

      res1 <- sum((rowSums((a * f)) / (ns - 1))^2)
      res2 <- sum((rowSums((a * a)) / (ns - 1)) *
                    (rowSums((f * f)) / (ns - 1)))
      r2t_y <- r2t_y + res1 / res2

      for (j in seq_len(ns)) {
        switch(
          class(model_obj$X)[1L],
          matrix = {
            cM <- cMA
          },
          list = {
            cM <- cMA[[j]]
          })

        ftotal <- Matrix::crossprod(Beta[, j], mm(cM, Beta[, j]))
        fixed1[j] <- fixed1[j] + ftotal

        for (k in seq_len(ngroups)) {
          sel <- (group == k)
          fpart <- Matrix::crossprod(
            Beta[sel, j], mm(cM[sel, sel], Beta[sel, j]))
          fixedsplit1[j, k] <- fixedsplit1[j, k] + fpart
        }
      }

      for (level in seq_len(nr)) {
        Lambda <- Lambdas[[level]]
        nf <- dim(Lambda)[[1]]
        for (factor in seq_len(nf)) {
          random1[, level] <- random1[, level] +
            t(Lambda[factor, ]) * Lambda[factor, ]
        }
      }

      if (nr > 0) {
        tot <- fixed1 + rowSums(random1)
        fixed <- fixed + fixed1 / tot
        for (level in seq_len(nr)) {
          random[, level] <- random[, level] + random1[, level] / tot
        }
      } else {
        fixed <- fixed + matrix(1, nrow = ns, ncol = 1)
      }

      for (k in seq_len(ngroups)) {
        fixedsplit[, k] <- fixedsplit[, k] +
          fixedsplit1[, k] / rowSums(fixedsplit1)
      }
    }

    fixed <- fixed / pool_n
    random <- random / pool_n
    fixedsplit <- fixedsplit / pool_n
    r2t_y <- r2t_y / pool_n
    r2t_beta <- r2t_beta / pool_n
  }

  # # .................................................................... ###

  vals <- matrix(0, nrow = ngroups + nr, ncol = ns)
  for (i in seq_len(ngroups)) {
    vals[i, ] <- fixed * fixedsplit[, i]
  }
  for (i in seq_len(nr)) {
    vals[ngroups + i, ] <- random[, i]
  }

  names(r2t_beta) <- model_obj$covNames
  leg <- group_names
  for (r in seq_len(nr)) {
    leg <- c(leg, paste0("Random: ", model_obj$rLNames[r]))
  }

  vals <- data.frame(vals) %>%
    tibble::tibble() %>%
    stats::setNames(model_obj$spNames) %>%
    dplyr::mutate(variable = leg, .before = 1)

  vp <- list(
    vals = vals,
    species_names = model_obj$spNames,
    R2T = list(Beta = r2t_beta, Y = r2t_y),
    group = group,
    groupnames = group_names)

  # # .................................................................... ###

  # Save the results
  ecokit::cat_time("Save the variance partitioning results", verbose = verbose)
  ecokit::save_as(object = vp, object_name = vp_file, out_path = file_var_par)

  vp$file <- file_var_par

  # # .................................................................... ###

  if (temp_cleanup) {

    ecokit::cat_time("Clean up temporary files", verbose = verbose)

    if (use_tf) {

      if (fs::dir_exists(path_lbeta)) {
        try({
          ecokit::system_command(
            paste0("rm -rf ", ecokit::normalize_path(path_lbeta)),
            ignore.stderr = TRUE, ignore.stdout = TRUE)
        },
        silent = TRUE)
      }
      if (fs::dir_exists(path_postList)) {
        try({
          ecokit::system_command(
            paste0("rm -rf ", ecokit::normalize_path(path_postList)),
            ignore.stderr = TRUE, ignore.stdout = TRUE)
        },
        silent = TRUE)
      }
    }

    temp_dir <- ecokit::normalize_path(temp_dir)
    try(
      expr = {
        file_paths <- list.files(
          path = temp_dir,
          pattern = "(vp_).+(feather|log)$", full.names = TRUE)

        fs::file_delete(file_paths)
      },
      silent = TRUE)

    if (fs::dir_exists(temp_dir) && length(fs::dir_ls(temp_dir)) == 0) {
      try({
        ecokit::system_command(
          paste0("rm -rf ", temp_dir),
          ignore.stderr = TRUE, ignore.stdout = TRUE)
      },
      silent = TRUE)
    }
  }

  # # .................................................................... ###

  ecokit::cat_diff(init_time = .start_time, verbose = verbose)

  doParallel::stopImplicitCluster()

  return(vp)
}
