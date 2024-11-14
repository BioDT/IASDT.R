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
    NCores = 6, UseTF = TRUE, TF_Environ = NULL, TF_use_single = TRUE,
    Verbose = TRUE, OutFileName = "VarPar") {

  # # .................................................................... ###

  if (isFALSE(Verbose)) {
    sink(file = nullfile())
    on.exit(sink(), add = TRUE)
  }

  .StartTime <- lubridate::now(tzone = "CET")

  VarPar_Parallel <- VarPar_Sequential <- is_gpu_available <-
    check_modules <- NULL

  # # .................................................................... ###

  # Load the model object
  if (is.null(Path_Model) || !file.exists(Path_Model)) {
    stop(
      paste0("Model path is NULL or does not exist: ", Path_Model),
      call. = FALSE)
  }

  hM <- IASDT.R::LoadAs(Path_Model)

  # # .................................................................... ###

  Path_VarPar <- file.path(
    dirname(dirname(Path_Model)), "Model_Postprocessing/Variance_Partitioning")
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
  names_to_remove <- c(
    "postList", "Y", "XScaled", "rL", "ranLevels", "XData", "dfPi",
    "studyDesign", "C", "Pi", "phyloTree", "XFormula", "XScalePar",
    "aSigma", "bSigma", "TrScaled", "YScalePar", "call", "rhopw",
    "distr", "V0", "UGamma", "YScaled")
  hM[names_to_remove] <- NULL

  invisible(gc())

  # # .................................................................... ###

  if (UseTF) {

    IASDT.R::CatTime("Check TensorFlow settings", Level = 1)

    # Check if TF_Environ directory exists
    if (is.null(TF_Environ) || !dir.exists(TF_Environ)) {
      stop(
        paste0(
          "When `UseTF` is TRUE, `TF_Environ` must be specified and should ",
          "point to an existing directory with a Python environment"),
        call. = FALSE)
    }

    PythonScripts <- c(
      system.file("VarPar.py", package = "IASDT.R"),
      system.file("Utilities.py", package = "IASDT.R"))

    # Check if PythonScripts exists
    if (!all(file.exists(PythonScripts))) {
      stop(
        "Necessary Python scripts do not exist in the package files",
        call. = FALSE)
    }

    # Suppress TensorFlow warnings and disable optimizations
    Sys.setenv(TF_CPP_MIN_LOG_LEVEL = "3", TF_ENABLE_ONEDNN_OPTS = "0")

    # Activate the python environment
    reticulate::use_virtualenv(TF_Environ, required = TRUE)
    purrr::walk(PythonScripts, reticulate::source_python)

    # Check if all necessary python modules exist
    MissingModules <- check_modules(
      module_list = c(
        "os", "tensorflow", "numpy", "logging",
        "concurrent.futures", "functools"),
      print_status = TRUE)

    if (length(MissingModules) > 0) {
      stop(
        paste0(
          "The following module(s) are missing: ",
          paste0(MissingModules, collapse = "; ")),
        call. = FALSE)
    } else {
      IASDT.R::CatTime(
        "All necessary python modules are available", Level = 2)
    }

    # Check GPU availability
    GPU <- is_gpu_available(print_status = FALSE)

    if (GPU) {
      IASDT.R::CatTime("GPU is available", Level = 2)
    } else {
      IASDT.R::CatTime("No GPU is available", Level = 2)
    }


    if (NCores > 1) {
      results <- VarPar_Parallel(
        hM_X = hM$X, hM_Tr = hM$Tr,
        postList = postList, use_single = TF_use_single,
        num_threads = NCores)
    } else {
      results <- VarPar_Sequential(
        hM_X = hM$X, hM_Tr = hM$Tr,
        postList = postList, use_single = TF_use_single)
    }

    la <- purrr::map(results, ~ .x$la)
    lf <- purrr::map(results, ~ .x$lf)
    lmu <- purrr::map(results, ~ .x$lmu)
    lbeta <- purrr::map(results, ~ .x$lbeta)

    rm(results)
    invisible(gc())

  } else {

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

    gemu <- function(a) {
      res <- t(hM$Tr %*% t(a$Gamma))
      return(res)
    }
    lmu <- lapply(postList, gemu)

    gebeta <- function(a) {
      res <- a$Beta
      return(res)
    }
    lbeta <- lapply(postList, gebeta)
  }

  # # .................................................................... ###

  # Remove Gamma from postList
  postList <- purrr::map(
    postList,
    ~ {
      .x["Gamma"] <- NULL
      .x
    })

  # # .................................................................... ###

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
  File_VarPar <- file.path(Path_VarPar, paste0(OutFileName, ".RData"))
  IASDT.R::SaveAs(InObj = VP, OutObj = OutFileName, OutPath = File_VarPar)

  VP$File <- File_VarPar

  # # .................................................................... ###

  IASDT.R::CatDiff(InitTime = .StartTime)
  return(VP)
}
