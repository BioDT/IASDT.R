## |------------------------------------------------------------------------| #
# VarPar_Compute ----
## |------------------------------------------------------------------------| #

#' Computes variance partitioning of Hmsc models
#'
#' Computes variance components with respect to given grouping of fixed effects
#' and levels of random effects. This function is a wrapper around the
#' `Hmsc::computeVariancePartitioning`, but with the added functionality of
#' parallel computation using TensorFlow.
#'
#' @param Path_Model Character. Path to fitted `Hmsc` model object.
#' @param group vector of numeric values corresponding to group identifiers in
#'   groupnames. If the model was defined with `XData` and `XFormula`, the
#'   default is to use model terms.
#' @param groupnames vector of names for each group of fixed effect. Should
#'   match `group`. If the model was defined with `XData` and `XFormula`, the
#'   default is to use the labels of model terms.
#' @param start index of first MCMC sample included
#' @param na.ignore logical. If TRUE, covariates are ignored for sites where the
#'   focal species is NA when computing variance-covariance matrices for each
#'   species
#' @param NCores Integer specifying the number of parallel cores for
#'   parallelization. This is only used when `UseTF` is TRUE.
#' @param Chunk_size Integer. Number of samples to process in each chunk. This
#'   is only used when `UseTF` is TRUE.
#' @param Verbose Logical. Indicates whether progress messages should be
#'   displayed. Defaults to `TRUE`.
#' @param OutFileName Character. Name of the output file to save the results.
#' @export
#' @name VarPar_Compute
#' @author Ahmed El-Gabbas
#' @inheritParams Predict_Hmsc
#' @inheritParams Predict_LF
#' @export

VarPar_Compute <- function(
    Path_Model, group = NULL, groupnames = NULL, start = 1, na.ignore = FALSE,
    NCores = 6, UseTF = TRUE, TF_Environ = NULL,
    TF_use_single = FALSE, Chunk_size = 50, Verbose = TRUE,
    OutFileName = "VarPar") {

  # # .................................................................... ###

  if (isFALSE(Verbose)) {
    sink(file = nullfile())
    on.exit(sink(), add = TRUE)
  }

  .StartTime <- lubridate::now(tzone = "CET")

  # # .................................................................... ###

  # Check if the virtual environment and Python scripts exist

  if (UseTF) {
    if (is.null(TF_Environ) || !dir.exists(TF_Environ)) {
      stop(
        paste0(
          "When `UseTF` is TRUE, `TF_Environ` must be specified and should ",
          "point to an existing directory with a Python environment"),
        call. = FALSE)
    }

    # Determine the Python executable path
    python_executable <- if (.Platform$OS.type == "windows") {
      file.path(TF_Environ, "Scripts", "python.exe")
    } else {
      file.path(TF_Environ, "bin", "python")
    }

    if (!file.exists(python_executable)) {
      stop(
        "Python executable not found in the virtual environment.",
        call. = FALSE)
    }

    # Paths to the Python scripts
    Script_geta <- system.file("VP_geta.py", package = "IASDT.R")
    Script_getf <- system.file("VP_getf.py", package = "IASDT.R")
    Script_gemu <- system.file("VP_gemu.py", package = "IASDT.R")
    if (!all(file.exists(c(Script_geta, Script_getf, Script_gemu)))) {
      stop(
        "Necessary python scripts do not exist",
        call. = FALSE)
    }
  }

  # # .................................................................... ###

  IASDT.R::CatTime("Load the model object")

  # Load the model object
  if (is.null(Path_Model) || !file.exists(Path_Model)) {
    stop(
      paste0("Model path is NULL or does not exist: ", Path_Model),
      call. = FALSE)
  }

  hM <- IASDT.R::LoadAs(Path_Model)

  # # .................................................................... ###

  Path_VarPar <- file.path(
    dirname(dirname(Path_Model)),
    "Model_Postprocessing/Variance_Partitioning")
  fs::dir_create(Path_VarPar)

  # # .................................................................... ###

  ny <- hM$ny
  nc <- hM$nc
  ns <- hM$ns
  nr <- hM$nr

  if (is.null(group)) {
    ## default: use terms
    if (nc > 1) {
      if (is.null(hM$XFormula)) {
        stop("no XFormula: you must give 'group' and 'groupnames'")
      }
      group <- attr(hM$X, "assign")
      if (group[1] == 0) { # assign (Intercept) to group
        group[1] <- 1
      }
      groupnames <- attr(
        stats::terms(hM$XFormula, data = hM$XData), "term.labels")
    } else {
      group <- c(1)
      groupnames <- hM$covNames[1]
    }
  }

  ngroups <- max(group)
  fixed <- matrix(0, nrow = ns, ncol = 1)
  fixedsplit <- matrix(0, nrow = ns, ncol = ngroups)
  random <- matrix(0, nrow = ns, ncol = nr)

  R2T.Y <- 0
  R2T.Beta <- rep(0, nc)

  # If na.ignore=T, convert XData to a list
  if (na.ignore) {
    xl <- list()
    for (s in 1:ns) {
      xl[[s]] <- hM$X
    }
    hM$X <- xl
  }

  switch(
    class(hM$X)[1L],
    matrix = {
      cMA <- stats::cov(hM$X)
    },
    list = {
      if (na.ignore) {
        cMA <- list()
        for (s in 1:ns) {
          cMA[[s]] <- stats::cov(hM$X[[s]][which(hM$Y[, s] > -Inf), ])
        }
      } else {
        cMA <- lapply(hM$X, stats::cov)
      }
    }
  )

  IASDT.R::CatTime("Prepare postList")
  postList <- Hmsc::poolMcmcChains(hM$postList, start = start) %>%
    purrr::map(
      ~ {
        Items2Delete <- c(
          "Eta", "Psi", "V", "sigma", "Delta", "Alpha",
          "rho", "wRRR", "PsiRRR", "DeltaRRR")
        .x[Items2Delete] <- NULL
        return(.x)
      }
    )

  # Remove unnecessary elements from the model object
  IASDT.R::CatTime("Remove model object elements")
  names_to_remove <- c(
    "postList", "Y", "XScaled", "rL", "ranLevels", "XData", "dfPi",
    "studyDesign", "C", "Pi", "phyloTree", "XFormula", "XScalePar",
    "aSigma", "bSigma", "TrScaled", "YScalePar", "call", "rhopw",
    "distr", "V0", "UGamma", "YScaled")
  hM[names_to_remove] <- NULL

  invisible(gc())

  # # .................................................................... ###

  if (UseTF) {

    IASDT.R::CatTime("Prepare la/lf/lmu lists using TensorFlow")

    # Create the temporary directory
    Path_Temp <- file.path(dirname(dirname(Path_Model)), "TEMP_VP")
    fs::dir_create(Path_Temp)

    # Prepare data
    IASDT.R::CatTime("Prepare data for TensorFlow", Level = 1)

    # X
    IASDT.R::CatTime("X", Level = 2)
    Path_X <- file.path(Path_Temp, "VP_X.feather")
    arrow::write_feather(as.data.frame(hM$X), Path_X)

    # Tr
    IASDT.R::CatTime("Tr", Level = 2)
    Path_Tr <- file.path(Path_Temp, "VP_Tr.feather")
    arrow::write_feather(as.data.frame(hM$Tr), Path_Tr)

    # Gamma - convert each list item into a column in a data frame
    IASDT.R::CatTime("Gamma", Level = 2)
    Path_Gamma <- file.path(Path_Temp, "VP_Gamma.feather")
    Gamma_data <- postList %>%
      purrr::map(~as.vector(.x[["Gamma"]])) %>%
      as.data.frame() %>%
      stats::setNames(paste0("Sample_", seq_len(ncol(.))))
    arrow::write_feather(Gamma_data, Path_Gamma)

    # Beta -- Each element of Beta is a matrix, so each list item is saved to
    # separate feather file
    IASDT.R::CatTime("Beta", Level = 2)
    purrr::walk(
      .x = seq_along(postList),
      .f = ~ {
        Beta <- postList[[.x]][["Beta"]] %>%
          stats::setNames(paste0("Sample_", .x)) %>%
          as.data.frame()

        arrow::write_feather(
          x = Beta,
          sink = file.path(
            Path_Temp, paste0("VP_Beta_", sprintf("%04d", .x), ".feather")))
        return(NULL)
      }
    )

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Running geta
    IASDT.R::CatTime("Running geta", Level = 1)

    Path_Out_a <- file.path(Path_Temp, "VP_A.feather")
    cmd_a <- paste(
      python_executable, Script_geta,
      "--tr", Path_Tr,
      "--x", Path_X,
      "--gamma", Path_Gamma,
      "--output", Path_Out_a,
      "--ncores", NCores,
      "--chunk_size", Chunk_size)

    if (TF_use_single) {
      cmd_a <- c(cmd_a, "--use_single")
    }

    # Run the command using system or system2
    la <- system(cmd_a, wait = TRUE, intern  = TRUE)

    # Check for errors
    if (!inherits(la, "error") || length(la) != 0 || la[length(la)] == "Done") {
      stop(paste0("Error in computing geta: ", la), call. = FALSE)
    }

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Running getf
    IASDT.R::CatTime("Running getf", Level = 1)

    Path_Out_f <- file.path(Path_Temp, "VP_F.feather")

    cmd_f <- paste(
      python_executable, Script_getf,
      "--x", shQuote(Path_X),
      "--beta_dir", shQuote(Path_Temp),
      "--output", shQuote(Path_Out_f),
      "--ncores", NCores)

    if (TF_use_single) {
      cmd_f <- c(cmd_f, "--use_single")
    }

    # Run the command using system or system2
    lf <- system(cmd_f, wait = TRUE, intern  = TRUE)
    # Check for errors
    if (!inherits(lf, "error") || length(lf) != 0 || lf[length(lf)] == "Done") {
      stop(paste0("Error in computing getf: ", lf), call. = FALSE)
    }

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Running gemu
    IASDT.R::CatTime("Running gemu", Level = 1)

    Path_Out_mu <- file.path(Path_Temp, "VP_Mu.feather")

    cmd_mu <- paste(
      python_executable, Script_gemu,
      "--tr", Path_Tr,
      "--gamma", Path_Gamma,
      "--output", Path_Out_mu,
      "--ncores", NCores,
      "--chunk_size", Chunk_size)

    if (TF_use_single) {
      cmd_mu <- c(cmd_mu, "--use_single")
    }

    # Run the command using system or system2
    lmu <- system(cmd_mu, wait = TRUE, intern  = TRUE)
    if (!inherits(lmu, "error") || length(lmu) != 0 ||
        lmu[length(lmu)] == "Done") {
      stop(paste0("Error in computing gemu: ", lmu), call. = FALSE)
    }

    # |||||||||||||||||||||||||||||||||||||||||||||||||||||||

    # Reading the results frpm feather files
    IASDT.R::CatTime("Reading the results frpm feather files", Level = 1)

    IASDT.R::CatTime("Prepare working on parallel", Level = 2)
    withr::local_options(
      future.globals.maxSize = 8000 * 1024^2, future.gc = TRUE)

    if (NCores == 1) {
      future::plan("future::sequential", gc = TRUE)
    } else {
      c1 <- snow::makeSOCKcluster(NCores)
      on.exit(try(snow::stopCluster(c1), silent = TRUE), add = TRUE)
      future::plan("future::cluster", workers = c1, gc = TRUE)
      on.exit(future::plan("future::sequential", gc = TRUE), add = TRUE)
    }

    IASDT.R::CatTime("la data", Level = 2)
    la <- furrr::future_map(
      .x = seq_along(postList),
      .f = ~ {
        stringr::str_pad(string = .x, width = 4, pad = "0") %>%
          paste0("_", ., ".feather") %>%
          stringr::str_replace(Path_Out_a, ".feather", .) %>%
          arrow::read_feather() %>%
          as.matrix(.x) %>%
          unname() %>%
          return()
      },
      .options = furrr::furrr_options(seed = TRUE, scheduling = Inf))

    IASDT.R::CatTime("lf data", Level = 2)
    lf <- furrr::future_map(
      .x = seq_along(postList),
      .f = ~ {
        stringr::str_pad(string = .x, width = 4, pad = "0") %>%
          paste0("_", ., ".feather") %>%
          stringr::str_replace(Path_Out_f, ".feather", .) %>%
          arrow::read_feather() %>%
          as.matrix(.x) %>%
          unname() %>%
          return()
      },
      .options = furrr::furrr_options(seed = TRUE, scheduling = Inf))

    IASDT.R::CatTime("lmu data", Level = 2)
    lmu <- furrr::future_map(
      .x = seq_along(postList),
      .f = ~ {
        stringr::str_pad(string = .x, width = 4, pad = "0") %>%
          paste0("_", ., ".feather") %>%
          stringr::str_replace(Path_Out_mu, ".feather", .) %>%
          arrow::read_feather() %>%
          as.matrix(.x) %>%
          unname() %>%
          return()
      },
      .options = furrr::furrr_options(seed = TRUE, scheduling = Inf))

    IASDT.R::CatTime("stop cluster", Level = 2)
    snow::stopCluster(c1)
    future::plan("future::sequential", gc = TRUE)
    invisible(gc())

    fs::dir_delete(Path_Temp)

  } else {

    IASDT.R::CatTime("Prepare la/lf/lmu lists using original R code")

    IASDT.R::CatTime("Running geta", Level = 1)
    geta <- function(a) {
      switch(
        class(hM$X)[1L],
        matrix = {
          res <- hM$X %*% (t(hM$Tr %*% t(a$Gamma)))
        },
        list = {
          res <- matrix(NA, ny, ns)
          for (j in 1:hM$ns) {
            res[, j] <- hM$X[[j]] %*% (t(hM$Tr[j, ] %*% t(a$Gamma)))
          }
        })
      return(res)
    }
    la <- lapply(postList, geta)

    IASDT.R::CatTime("Running getf", Level = 1)
    getf <- function(a) {
      switch(
        class(hM$X)[1L],
        matrix = {
          res <- hM$X %*% (a$Beta)
        },
        list = {
          res <- matrix(NA, ny, ns)
          for (j in 1:ns) res[, j] <- hM$X[[j]] %*% a$Beta[, j]
        })
      return(res)
    }
    lf <- lapply(postList, getf)

    IASDT.R::CatTime("Running gemu", Level = 1)
    gemu <- function(a) {
      res <- t(hM$Tr %*% t(a$Gamma))
      return(res)
    }
    lmu <- lapply(postList, gemu)

  }

  # # .................................................................... ###

  IASDT.R::CatTime("Running gebeta", Level = 1)
  gebeta <- function(a) {
    res <- a$Beta
    return(res)
  }
  lbeta <- lapply(postList, gebeta)

  # # .................................................................... ###

  # Remove Gamma from postList
  IASDT.R::CatTime("Remove Gamma from postList")
  postList <- purrr::map(
    postList,
    ~ {
      .x["Gamma"] <- NULL
      .x
    })

  # # .................................................................... ###

  IASDT.R::CatTime("Compute variance partitioning")

  poolN <- length(postList) # pooled chains

  for (i in seq_len(poolN)) {

    # Suppress warnings when no trait information is used in the models
    # cor(Beta[k, ], lmu[k, ]) : the standard deviation is zero
    suppressWarnings({
      for (k in 1:nc) {
        R2T.Beta[k] <- R2T.Beta[k] +
          stats::cor(lbeta[[i]][k, ], lmu[[i]][k, ])^2
      }
    })

    fixed1 <- matrix(0, nrow = ns, ncol = 1)
    fixedsplit1 <- matrix(0, nrow = ns, ncol = ngroups)
    random1 <- matrix(0, nrow = ns, ncol = nr)
    Beta <- postList[[i]]$Beta
    Lambdas <- postList[[i]]$Lambda

    f <- lf[[i]]
    a <- la[[i]]

    a <- a - matrix(rep(rowMeans(a), ns), ncol = ns)
    f <- f - matrix(rep(rowMeans(f), ns), ncol = ns)
    res1 <- sum((rowSums((a * f)) / (ns - 1))^2)
    res2 <- sum((rowSums((a * a)) / (ns - 1)) *
                  (rowSums((f * f)) / (ns - 1)))
    R2T.Y <- R2T.Y + res1 / res2

    for (j in 1:ns) {
      switch(
        class(hM$X)[1L],
        matrix = {
          cM <- cMA
        },
        list = {
          cM <- cMA[[j]]
        }
      )
      ftotal <- t(Beta[, j]) %*% cM %*% Beta[, j]
      fixed1[j] <- fixed1[j] + ftotal
      for (k in 1:ngroups) {
        sel <- (group == k)
        fpart <- t(Beta[sel, j]) %*% cM[sel, sel] %*% Beta[sel, j]
        fixedsplit1[j, k] <- fixedsplit1[j, k] + fpart
      }
    }

    for (level in seq_len(nr)) {
      Lambda <- Lambdas[[level]]
      nf <- dim(Lambda)[[1]]
      for (factor in 1:nf) {
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
    for (k in 1:ngroups) {
      fixedsplit[, k] <- fixedsplit[, k] +
        fixedsplit1[, k] / rowSums(fixedsplit1)
    }
  }

  fixed <- fixed / poolN
  random <- random / poolN
  fixedsplit <- fixedsplit / poolN
  R2T.Y <- R2T.Y / poolN
  R2T.Beta <- R2T.Beta / poolN

  vals <- matrix(0, nrow = ngroups + nr, ncol = ns)
  for (i in 1:ngroups) {
    vals[i, ] <- fixed * fixedsplit[, i]
  }
  for (i in seq_len(nr)) {
    vals[ngroups + i, ] <- random[, i]
  }

  names(R2T.Beta) <- hM$covNames
  colnames(vals) <- hM$spNames
  leg <- groupnames
  for (r in seq_len(nr)) {
    leg <- c(leg, paste("Random: ", hM$rLNames[r], sep = ""))
  }
  rownames(vals) <- leg

  VP <- list()
  VP$vals <- vals
  VP$R2T <- list(Beta = R2T.Beta, Y = R2T.Y)
  VP$group <- group
  VP$groupnames <- groupnames

  # # .................................................................... ###

  # Save the results
  IASDT.R::CatTime("Save the variance partitioning results")

  File_VarPar <- file.path(Path_VarPar, paste0(OutFileName, ".RData"))
  IASDT.R::SaveAs(InObj = VP, OutObj = OutFileName, OutPath = File_VarPar)

  VP$File <- File_VarPar

  # # .................................................................... ###

  IASDT.R::CatDiff(InitTime = .StartTime)

  return(VP)
}
