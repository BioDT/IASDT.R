## |------------------------------------------------------------------------| #
# variance_partitioning_compute ----
## |------------------------------------------------------------------------| #

#' Computes and visualise variance partitioning of Hmsc models
#'
#' The **`variance_partitioning_compute()`** function computes variance
#' components with respect to given grouping of fixed effects and levels of
#' random effects. This function inherits the main functionality from the
#' `Hmsc::computeVariancePartitioning` function, but with the added
#' functionality of parallel computation and using `TensorFlow`.<br/>The
#' **`variance_partitioning_plot()`** function generates plots for variance
#' partitioning as JPEG files. It allows for sorting the predictors and species;
#' e.g., by the mean value per predictor; and by original species order. It also
#' plots the raw variance partitioning (relative variance partitioning
#' multiplied by the Tjur-R<sup>2</sup> value).
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
#'   partitioning using `TensorFlow`. This is only effective when `use_TF` is
#'   `TRUE`. Default: `1`.
#' @param strategy Character. The parallel processing strategy to use. Valid
#'   options are "sequential", "multisession" (default), "multicore", and
#'   "cluster". See [future::plan()] and [ecokit::set_parallel()] for details.
#' @param chunk_size Integer. Size of each chunk of samples to process in
#'   parallel. Only relevant for `TensorFlow`. Default: `50`.
#' @param verbose Logical. Whether to print progress messages. Default: `TRUE`.
#' @param temp_cleanup Logical. Whether to delete temporary files after
#'   processing. Default: `TRUE`.
#' @param VP_commands_only Logical. If `TRUE`, returns the commands to run the
#'   Python script. Default is `FALSE`. Only relevant when `use_TF` is `TRUE`.
#' @param VP_file Character. Name of the output file to save the results.
#'   Default: `VarPar`.
#' @param env_file Character. Path to the environment file containing paths to
#'   data sources. Defaults to `.env`.
#' @param width,height Numeric. Width and height of the output plot in
#'   centimetres. Default: `30` and `15`, respectively.
#' @param Axis_text Numeric. Size of the axis text. Default: `4`.
#' @author Ahmed El-Gabbas
#' @inheritParams predict_hmsc
#' @inheritParams predict_latent_factor
#' @name variance_partitioning
#' @rdname variance_partitioning
#' @order 1
#' @export

variance_partitioning_compute <- function(
    path_model, group = NULL, group_names = NULL, start = 1L, na.ignore = FALSE,
    n_cores = 8L, strategy = "multisession", use_TF = TRUE, TF_environ = NULL,
    TF_use_single = FALSE, temp_cleanup = TRUE, chunk_size = 50L,
    verbose = TRUE, VP_file = "VarPar", VP_commands_only = FALSE) {

  if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores <= 0) {
    ecokit::stop_ctx(
      "n_cores must be a single positive integer.", n_cores = n_cores,
      include_backtrace = TRUE)
  }

  if (!is.character(strategy)) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector",
      strategy = strategy, class_strategy = class(strategy))
  }
  if (strategy == "sequential") {
    n_cores <- 1L
  }
  if (length(strategy) != 1L) {
    ecokit::stop_ctx(
      "`strategy` must be a character vector of length 1",
      strategy = strategy, length_strategy = length(strategy))
  }
  valid_strategy <- c("sequential", "multisession", "multicore", "cluster")
  if (!strategy %in% valid_strategy) {
    ecokit::stop_ctx("Invalid `strategy` value", strategy = strategy)
  }

  # # .................................................................... ###

  .start_time <- lubridate::now(tzone = "CET")

  # # ..................................................................... ###

  # packages to be loaded in parallel
  pkg_to_export <- ecokit::load_packages_future(
    packages = c(
      "Matrix", "dplyr", "arrow", "purrr", "IASDT.R", "qs2", "methods",
      "stringr", "fs", "ecokit", "magrittr", "Hmsc"),
    strategy = strategy)

  # # .................................................................... ###

  # Check if the virtual environment and Python scripts exist

  if (use_TF) {

    # Check python virtual environment
    if (isFALSE(VP_commands_only) && .Platform$OS.type == "windows" &&
        (is.null(TF_environ) || !dir.exists(TF_environ))) {
      ecokit::stop_ctx(
        paste0(
          "When running on Windows and `use_TF` is TRUE, `TF_environ` must ",
          "be specified and point to an existing directory with a ",
          "Python virtual environment"),
        VP_commands_only = VP_commands_only, OS = .Platform$OS.type,
        TF_environ = TF_environ, include_backtrace = TRUE)
    }

    # Determine the Python executable path

    # On Windows, the TF calculations has to be done through a valid virtual
    # environment; the path to the virtual environment must be specified in
    # `TF_environ`. On LUMI, this is not needed as the compatible python
    # installation is loaded automatically when loading `tensorflow` module.
    # When using another HPC system, the function needs to be adapted
    # accordingly to point to a valid python virtual environment.

    if (.Platform$OS.type == "windows") {
      python_executable <- fs::path(TF_environ, "Scripts", "python.exe")

      if (isFALSE(VP_commands_only) && !file.exists(python_executable)) {
        ecokit::stop_ctx(
          "Python executable not found in the virtual environment",
          python_executable = python_executable, include_backtrace = TRUE)
      }
    } else {
      python_executable <- "/usr/bin/time -v python3"
    }

    # Check GPU availability
    if (isFALSE(VP_commands_only) && .Platform$OS.type == "windows") {
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
    Script_geta <- system.file("VP_geta.py", package = "IASDT.R")
    Script_getf <- system.file("VP_getf.py", package = "IASDT.R")
    Script_gemu <- system.file("VP_gemu.py", package = "IASDT.R")
    if (!all(file.exists(c(Script_geta, Script_getf, Script_gemu)))) {
      ecokit::stop_ctx(
        "Necessary python scripts do not exist",
        Script_geta = Script_geta, Script_getf = Script_getf,
        Script_gemu = Script_gemu, include_backtrace = TRUE)
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

  Model <- ecokit::load_as(path_model)

  # # .................................................................... ###

  # Create folder for variance partitioning results
  Path_VarPar <- fs::path(
    dirname(dirname(path_model)),
    "Model_Postprocessing", "Variance_Partitioning")
  fs::dir_create(Path_VarPar)

  # # .................................................................... ###

  ny <- Model$ny
  nc <- Model$nc
  ns <- Model$ns
  nr <- Model$nr

  if (is.null(group)) {
    # names of variables used in the model
    group_names <- names(Model$XData)

    # actual variables used in the model, including quadratic terms
    ModelVars <- dimnames(Model$X)[[2]][-1]

    # group variable to combine variable and its quadratic terms together
    group <- purrr::map(
      ModelVars, ~ which(stringr::str_detect(.x, group_names))) %>%
      unlist() %>%
      # add intercept to the first group
      c(.[1], .)
  }

  # If na.ignore=T, convert XData to a list
  if (na.ignore) {
    xl <- list()
    for (s in seq_len(ns)) {
      xl[[s]] <- Model$X
    }
    Model$X <- xl
  }

  switch(
    class(Model$X)[1L],
    matrix = {
      cMA <- stats::cov(Model$X)
    },
    list = {
      if (na.ignore) {
        cMA <- list()
        for (s in seq_len(ns)) {
          cMA[[s]] <- stats::cov(Model$X[[s]][which(Model$Y[, s] > -Inf), ])
        }
      } else {
        cMA <- lapply(Model$X, stats::cov)
      }
    }
  )

  # # .................................................................... ###

  # Prepare postList-----

  ecokit::cat_time("Prepare postList", verbose = verbose)
  postList <- Hmsc::poolMcmcChains(Model$postList, start = start)

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
  # Model[names_to_remove] <- NULL

  # `trim_hmsc` is much faster
  Model <- IASDT.R::trim_hmsc(model = Model, names_to_remove = Items2Delete)
  invisible(gc())

  poolN <- length(postList)
  ngroups <- max(group)

  # # .................................................................... ###
  # # .................................................................... ###

  # Prepare `la`/`lf`/`lmu` lists -----

  if (use_TF) {

    ecokit::cat_time("Prepare VP files for `TensorFlow`", verbose = verbose)

    # Create the temporary directory
    Path_Temp <- fs::path(dirname(dirname(path_model)), "TEMP_VP")
    fs::dir_create(Path_Temp)
    path_VP_input_files <- fs::path(Path_Temp, "VP_input_files.txt")

    FileSuffix <- stringr::str_pad(
      string = seq_len(poolN), pad = "0", width = 4)

    if (n_cores == 1) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = n_cores, strategy = strategy, show_log = FALSE)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    # List of feather files resulted from `geta` function
    Files_la <- fs::path(Path_Temp, paste0("VP_A_", FileSuffix, ".feather"))
    Files_la_Exist <- all(fs::file_exists(Files_la))
    if (Files_la_Exist) {
      Files_la_Exist <- future.apply::future_lapply(
        X = Files_la,
        FUN = ecokit::check_data, warning = FALSE,
        future.seed = TRUE, future.packages = c("arrow", "ecokit")) %>%
        unlist() %>%
        all()
    }

    # List of feather files resulted from `getf` function
    Files_lf <- fs::path(Path_Temp, paste0("VP_F_", FileSuffix, ".feather"))
    Files_lf_Exist <- all(fs::file_exists(Files_lf))
    if (Files_lf_Exist) {
      Files_lf_Exist <- future.apply::future_lapply(
        X = Files_lf, FUN = ecokit::check_data, warning = FALSE,
        future.seed = TRUE, future.packages = c("arrow", "ecokit")) %>%
        unlist() %>%
        all()
    }

    # List of feather files resulted from `gemu` function
    Files_lmu <- fs::path(Path_Temp, paste0("VP_Mu_", FileSuffix, ".feather"))
    Files_lmu_Exist <- all(fs::file_exists(Files_lmu))
    if (Files_lmu_Exist) {
      Files_lmu_Exist <- future.apply::future_lapply(
        X = Files_lmu, FUN = ecokit::check_data, warning = FALSE,
        future.seed = TRUE, future.packages = c("arrow", "ecokit")) %>%
        unlist() %>%
        all()
    }

    Beta_Files <- fs::path(
      Path_Temp, paste0("VP_Beta_", FileSuffix, ".feather"))

    if (all(c(Files_la_Exist, Files_lf_Exist, Files_lmu_Exist))) {

      # All feather data are already processed and available on disk
      ecokit::cat_time(
        "Data for `la`/`lf`/`lmu` lists were already processed",
        verbose = verbose)

    } else {

      # Write the contents of Files_la, Files_lf, Files_lmu to a text file
      readr::write_lines(
        x = c(Files_la, Files_lf, Files_lmu), file = path_VP_input_files)

      ## Prepare la/lf/lmu lists using `TensorFlow` ----
      ecokit::cat_time(
        "Prepare la/lf/lmu lists using `TensorFlow`", verbose = verbose)

      ### Prepare data for `TensorFlow` ----
      ecokit::cat_time(
        "Prepare data for `TensorFlow`", level = 1L, verbose = verbose)

      #### X data -----
      # needed only to calculate `geta` and `getf` functions
      if (!all(c(Files_la_Exist, Files_lf_Exist))) {
        ecokit::cat_time("X", level = 2L, verbose = verbose)
        Path_X <- fs::path(Path_Temp, "VP_X.feather")
        if (!file.exists(Path_X)) {
          arrow::write_feather(as.data.frame(Model$X), Path_X)
        }
      }

      #### Tr / Gamma ------
      # needed only to calculate `geta` and `gemu` functions
      if (!all(c(Files_la_Exist, Files_lmu_Exist))) {
        ecokit::cat_time("Tr", level = 2L, verbose = verbose)
        Path_Tr <- fs::path(Path_Temp, "VP_Tr.feather")
        if (!file.exists(Path_Tr)) {
          arrow::write_feather(as.data.frame(Model$Tr), Path_Tr)
        }

        # Gamma - convert each list item into a column in a data frame
        ecokit::cat_time("Gamma", level = 2L, verbose = verbose)
        Path_Gamma <- fs::path(Path_Temp, "VP_Gamma.feather")
        if (!file.exists(Path_Gamma)) {
          Gamma_data <- postList %>%
            purrr::map(~as.vector(.x[["Gamma"]])) %>%
            as.data.frame() %>%
            stats::setNames(paste0("Sample_", seq_len(ncol(.))))
          arrow::write_feather(Gamma_data, Path_Gamma)
        }
      }

      #### Beta -----
      # only needed for `getf` function
      if (!Files_lf_Exist) {

        # Beta -- Each element of Beta is a matrix, so each list item is saved
        # to separate feather file
        ecokit::cat_time("Beta", level = 2L, verbose = verbose)

        Beta_Files_Exist <- all(fs::file_exists(Beta_Files))
        if (Beta_Files_Exist) {
          Beta_Files_Exist <- future.apply::future_lapply(
            X = Beta_Files, FUN = ecokit::check_data, warning = FALSE,
            future.seed = TRUE, future.packages = c("arrow", "ecokit")) %>%
            unlist() %>%
            all()
        }

        if (!Beta_Files_Exist) {

          ecokit::cat_time(
            "Processing beta in parallel", level = 3L, verbose = verbose)

          Beta0 <- future.apply::future_lapply(
            X = seq_along(postList),
            FUN = function(x) {

              Beta_File <- Beta_Files[x]

              if (ecokit::check_data(Beta_File, warning = FALSE)) {
                return(NULL)
              }

              if (fs::file_exists(Beta_File)) {
                try(fs::file_delete(Beta_File), silent = TRUE)
              }

              Beta <- purrr::pluck(postList, x, "Beta") %>%
                as.data.frame()

              attempt <- 1
              repeat {

                arrow::write_feather(x = Beta, sink = Beta_File)
                Sys.sleep(1)

                if (ecokit::check_data(Beta_File, warning = FALSE)) {
                  break
                }

                if (attempt >= 5) {
                  ecokit::stop_ctx(
                    "Failed to create Beta file after multiple attempts",
                    Beta_File = Beta_File, attempt = attempt,
                    include_backtrace = TRUE)
                }

                attempt <- attempt + 1
              }
              return(NULL)
            },
            # Setting future.scheduling = Inf makes calculations sequential!
            # future.scheduling = Inf,
            future.seed = TRUE, future.packages = pkg_to_export,
            future.globals = c("Beta_Files", "postList"))

          rm(Beta0)
        }

      }

      invisible(gc())

      # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

      ### Processing geta -----

      if (Files_la_Exist) {

        ecokit::cat_time(
          "All `la` data were already available on disk",
          level = 1L, verbose = verbose)

      } else {

        ecokit::cat_time(
          "Processing `geta` function", level = 1L, verbose = verbose)
        Path_Out_a <- fs::path(Path_Temp, "VP_A.feather") %>%
          ecokit::normalize_path()

        cmd_a <- paste(
          python_executable, Script_geta,
          "--tr", ecokit::normalize_path(Path_Tr, must_work = TRUE),
          "--x", ecokit::normalize_path(Path_X, must_work = TRUE),
          "--gamma", ecokit::normalize_path(Path_Gamma, must_work = TRUE),
          "--output", Path_Out_a,
          "--chunk_size", chunk_size)

        if (TF_use_single) {
          cmd_a <- paste0(cmd_a, " --use_single")
        }

        if (VP_commands_only) {

          # Save command to file
          # Redirect results of time to a log file
          path_log_a <- stringr::str_replace(Path_Out_a, ".feather", ".log")
          cmd_a <- paste0(cmd_a, paste0(" >> ", path_log_a, " 2>&1"))
          readr::write_lines(
            x = cmd_a,
            file = fs::path(Path_Temp, "VP_A_Command.txt"), append = FALSE)

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

      if (Files_la_Exist) {

        ecokit::cat_time(
          "All `lf` data were already available on disk",
          level = 1L, verbose = verbose)

      } else {

        ecokit::cat_time(
          "Processing `getf` function", level = 1L, verbose = verbose)
        Path_Out_f <- fs::path(Path_Temp, "VP_F.feather") %>%
          ecokit::normalize_path()

        cmd_f <- paste(
          python_executable, Script_getf,
          "--x", ecokit::normalize_path(Path_X, must_work = TRUE),
          "--beta_dir", ecokit::normalize_path(Path_Temp, must_work = TRUE),
          "--output", Path_Out_f,
          "--ncores", 1)

        if (TF_use_single) {
          cmd_f <- paste0(cmd_f, " --use_single")
        }

        if (VP_commands_only) {

          # Save command to file
          # Redirect results of time to a log file
          path_log_f <- stringr::str_replace(Path_Out_f, ".feather", ".log")
          cmd_f <- paste0(cmd_f, paste0(" >> ", path_log_f, " 2>&1"))
          readr::write_lines(
            x = cmd_f,
            file = fs::path(Path_Temp, "VP_F_Command.txt"),
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

      if (Files_lmu_Exist) {

        ecokit::cat_time(
          "All `lmu` data were already available on disk",
          level = 1L, verbose = verbose)

      } else {

        ecokit::cat_time(
          "Processing `gemu` function", level = 1L, verbose = verbose)
        Path_Out_mu <- fs::path(Path_Temp, "VP_Mu.feather") %>%
          ecokit::normalize_path()

        cmd_mu <- paste(
          python_executable, Script_gemu,
          "--tr", ecokit::normalize_path(Path_Tr, must_work = TRUE),
          "--gamma", ecokit::normalize_path(Path_Gamma, must_work = TRUE),
          "--output", Path_Out_mu,
          "--ncores", 1,
          "--chunk_size", chunk_size)

        if (TF_use_single) {
          cmd_mu <- paste0(cmd_mu, " --use_single")
        }

        if (VP_commands_only) {

          # Save command to file
          # Redirect results of time to a log file
          path_log_mu <- stringr::str_replace(Path_Out_mu, ".feather", ".log")
          cmd_mu <- paste0(cmd_mu, paste0(" >> ", path_log_mu, " 2>&1"))
          readr::write_lines(
            x = cmd_mu,
            file = fs::path(Path_Temp, "VP_mu_Command.txt"),
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

    # stopping the cluster
    if (n_cores > 1) {
      ecokit::set_parallel(stop_cluster = TRUE, show_log = FALSE)
      future::plan("sequential", gc = TRUE)
    }

    invisible(gc())

    if (VP_commands_only) {
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
        class(Model$X)[1L],
        matrix = {
          res <- Model$X %*% (t(Model$Tr %*% t(a$Gamma)))
        },
        list = {
          res <- matrix(NA, ny, ns)
          for (j in seq_len(ns)) {
            res[, j] <- Model$X[[j]] %*% (t(Model$Tr[j, ] %*% t(a$Gamma)))
          }
        })
      return(res)
    }

    la <- lapply(postList, geta)

    invisible(gc())

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

    ## getf ------

    ecokit::cat_time("Running getf", level = 1L, verbose = verbose)
    getf <- function(a) {
      switch(
        class(Model$X)[1L],
        matrix = {
          res <- Model$X %*% (a$Beta)
        },
        list = {
          res <- matrix(NA, ny, ns)
          for (j in seq_len(ns)) res[, j] <- Model$X[[j]] %*% a$Beta[, j]
        })
      return(res)
    }

    lf <- lapply(postList, getf)

    invisible(gc())

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

    ## gemu ------

    ecokit::cat_time("Running gemu", level = 1L, verbose = verbose)
    gemu <- function(a) {
      res <- t(Model$Tr %*% t(a$Gamma))
      return(res)
    }

    lmu <- lapply(postList, gemu)

    invisible(gc())

  }

  # # .................................................................... ###

  # Running gebeta -------
  ecokit::cat_time("Running gebeta", verbose = verbose)
  gebeta <- function(a) {
    res <- a$Beta
    return(res)
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

  if (use_TF) {

    ecokit::cat_time(
      "Computing variance partitioning in parallel", verbose = verbose)

    ecokit::cat_time(
      "Split `lbeta` list into small qs2 files", level = 1L, verbose = verbose)
    path_lbeta <- fs::path(Path_Temp, "lbeta")
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
    path_postList <- fs::path(Path_Temp, "postList")
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

    Model <- IASDT.R::trim_hmsc(
      model = Model, names_to_remove = setdiff(names(Model), "X"))
    invisible(gc())

    # temporary for debugging
    ecokit::all_objects_sizes(greater_than = 1L, in_function = TRUE)

    if (n_cores == 1) {
      future::plan("sequential", gc = TRUE)
    } else {
      ecokit::set_parallel(
        n_cores = n_cores, level = 1L, future_max_size = 800L,
        strategy = strategy, cat_timestamp = FALSE)
      withr::defer(future::plan("sequential", gc = TRUE))
    }

    ecokit::cat_time("Processing in parallel", level = 1L, verbose = verbose)

    Res <- future.apply::future_lapply(
      X = seq_len(n_postList),
      FUN = function(i) {

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
        DT_lmu <- as.matrix(arrow::read_feather(Files_lmu[i]))      # nolint: object_name_linter
        curr_lbeta <- fs::path(      # nolint: object_name_linter
          path_lbeta,
          paste0(
            "lbeta_", stringr::str_pad(i, width = 4, pad = "0"), ".qs2")) %>%
          ecokit::load_as()
        R2T.Beta <- purrr::map_dbl(
          .x = seq_len(nc),
          .f = ~ {
            suppressWarnings(stats::cor(curr_lbeta[.x, ], DT_lmu[.x, ])^2)
          })

        rm(curr_lbeta, DT_lmu, envir = environment())
        invisible(gc())

        fixed1 <- matrix(0, nrow = ns, ncol = 1)
        fixedsplit1 <- matrix(0, nrow = ns, ncol = ngroups)
        random1 <- matrix(0, nrow = ns, ncol = nr)

        DT_la <- as.matrix(arrow::read_feather(Files_la[i]))
        DT_lf <- as.matrix(arrow::read_feather(Files_lf[i]))
        # a <- DT_la - matrix(rep(rowMeans(DT_la), ns), ncol = ns)
        # f <- DT_lf - matrix(rep(rowMeans(DT_lf), ns), ncol = ns)
        a <- DT_la - rowMeans(DT_la)
        f <- DT_lf - rowMeans(DT_lf)

        res1 <- sum((rowSums((a * f)) / (ns - 1))^2)
        res2 <- sum((rowSums((a * a)) / (ns - 1)) *
                      (rowSums((f * f)) / (ns - 1)))
        R2T.Y <- res1 / res2

        for (j in seq_len(ns)) {
          switch(
            class(Model$X)[1L],
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
          tot <- fixed1 + rowSums(random1)
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
        #   fixedsplit[, k] <- fixedsplit1[, k] / rowSums(fixedsplit1)
        # }

        fixedsplit <- fixedsplit1 / rowSums(fixedsplit1)

        return(
          list(
            fixed = fixed, random = random, fixedsplit = fixedsplit,
            R2T.Y = R2T.Y, R2T.Beta = R2T.Beta))

      },
      future.scheduling = Inf, future.seed = TRUE,
      future.packages = pkg_to_export,
      future.globals = c(
        "ngroups", "Files_la", "Files_lf", "Files_lmu", "path_lbeta", "nc",
        "Model", "path_postList", "ns", "nr", "cMA", "group"))

    # stopping the cluster
    if (n_cores > 1) {
      ecokit::set_parallel(
        stop_cluster = TRUE, level = 2L, cat_timestamp = FALSE)
      future::plan("sequential", gc = TRUE)
    }

    # Summarise the results
    ecokit::cat_time("Summarise the results", level = 1L, verbose = verbose)
    fixed <- Reduce("+", purrr::map(Res, ~ .x$fixed)) / poolN
    random <- Reduce("+", purrr::map(Res, ~ .x$random)) / poolN
    fixedsplit <- Reduce("+", purrr::map(Res, ~ .x$fixedsplit)) / poolN
    R2T.Y <- Reduce("+", purrr::map(Res, ~ .x$R2T.Y)) / poolN
    R2T.Beta <- Reduce("+", purrr::map(Res, ~ .x$R2T.Beta)) / poolN

  } else {

    ecokit::cat_time(
      "Computing variance partitioning sequentially", verbose = verbose)

    mm <- methods::getMethod("%*%", "Matrix")

    fixed <- matrix(0, nrow = ns, ncol = 1)
    fixedsplit <- matrix(0, nrow = ns, ncol = ngroups)
    random <- matrix(0, nrow = ns, ncol = nr)
    R2T.Y <- 0
    R2T.Beta <- rep(0, nc)

    for (i in seq_len(poolN)) {

      if (i %% 200 == 0) {
        ecokit::cat_time(
          sprintf("Processing iteration %d of %d", i, poolN), verbose = verbose)
      }

      DT_la <- la[[i]]
      DT_lf <- lf[[i]]
      DT_lmu <- lmu[[i]]

      # Suppress warnings when no trait information is used in the models
      # cor(Beta[k, ], lmu[k, ]) : the standard deviation is zero
      suppressWarnings({
        for (k in seq_len(nc)) {
          R2T.Beta[k] <- R2T.Beta[k] +
            stats::cor(lbeta[[i]][k, ], DT_lmu[k, ])^2
        }
      })

      fixed1 <- matrix(0, nrow = ns, ncol = 1)
      fixedsplit1 <- matrix(0, nrow = ns, ncol = ngroups)
      random1 <- matrix(0, nrow = ns, ncol = nr)
      Beta <- postList[[i]]$Beta
      Lambdas <- postList[[i]]$Lambda

      a <- DT_la - matrix(rep(rowMeans(DT_la), ns), ncol = ns)
      f <- DT_lf - matrix(rep(rowMeans(DT_lf), ns), ncol = ns)

      res1 <- sum((rowSums((a * f)) / (ns - 1))^2)
      res2 <- sum((rowSums((a * a)) / (ns - 1)) *
                    (rowSums((f * f)) / (ns - 1)))
      R2T.Y <- R2T.Y + res1 / res2

      for (j in seq_len(ns)) {
        switch(
          class(Model$X)[1L],
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

    fixed <- fixed / poolN
    random <- random / poolN
    fixedsplit <- fixedsplit / poolN
    R2T.Y <- R2T.Y / poolN
    R2T.Beta <- R2T.Beta / poolN
  }

  # # .................................................................... ###

  vals <- matrix(0, nrow = ngroups + nr, ncol = ns)
  for (i in seq_len(ngroups)) {
    vals[i, ] <- fixed * fixedsplit[, i]
  }
  for (i in seq_len(nr)) {
    vals[ngroups + i, ] <- random[, i]
  }

  names(R2T.Beta) <- Model$covNames
  colnames(vals) <- Model$spNames
  leg <- group_names
  for (r in seq_len(nr)) {
    leg <- c(leg, paste0("Random: ", Model$rLNames[r]))
  }
  rownames(vals) <- leg

  VP <- list()
  VP$vals <- vals
  VP$R2T <- list(Beta = R2T.Beta, Y = R2T.Y)
  VP$group <- group
  VP$groupnames <- group_names

  # # .................................................................... ###

  # Save the results
  ecokit::cat_time("Save the variance partitioning results", verbose = verbose)

  File_VarPar <- fs::path(Path_VarPar, paste0(VP_file, ".RData"))
  ecokit::save_as(object = VP, object_name = VP_file, out_path = File_VarPar)

  VP$File <- File_VarPar

  # # .................................................................... ###

  if (temp_cleanup) {

    ecokit::cat_time("Clean up temporary files", verbose = verbose)

    if (use_TF) {

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

    Path_Temp <- ecokit::normalize_path(Path_Temp)
    try(
      expr = {
        file_paths <- list.files(
          path = Path_Temp,
          pattern = "(VP_).+(feather|log)$", full.names = TRUE)

        fs::file_delete(file_paths)
      },
      silent = TRUE)

    if (fs::dir_exists(Path_Temp) && length(fs::dir_ls(Path_Temp)) == 0) {
      try({
        ecokit::system_command(
          paste0("rm -rf ", Path_Temp),
          ignore.stderr = TRUE, ignore.stdout = TRUE)
      },
      silent = TRUE)
    }
  }

  # # .................................................................... ###

  ecokit::cat_diff(init_time = .start_time, verbose = verbose)

  return(VP)
}
