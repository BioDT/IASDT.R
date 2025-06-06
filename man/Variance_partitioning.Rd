% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod_variance_partitioning_compute.R,
%   R/mod_variance_partitioning_plot.R
\name{variance_partitioning}
\alias{variance_partitioning}
\alias{variance_partitioning_compute}
\alias{variance_partitioning_plot}
\title{Computes and visualise variance partitioning of Hmsc models}
\usage{
variance_partitioning_compute(
  path_model,
  group = NULL,
  group_names = NULL,
  start = 1L,
  na.ignore = FALSE,
  n_cores = 8L,
  use_TF = TRUE,
  TF_environ = NULL,
  TF_use_single = FALSE,
  temp_cleanup = TRUE,
  chunk_size = 50L,
  verbose = TRUE,
  VP_file = "VarPar",
  VP_commands_only = FALSE
)

variance_partitioning_plot(
  path_model,
  env_file = ".env",
  VP_file = "VarPar",
  use_TF = TRUE,
  TF_environ = NULL,
  n_cores = 1L,
  strategy = "multisession",
  width = 30,
  height = 15,
  Axis_text = 4
)
}
\arguments{
\item{path_model}{Character. Path to fitted \code{Hmsc} model object.}

\item{group}{vector of numeric values corresponding to group identifiers in
groupnames. If the model was defined with \code{XData} and \code{XFormula}, the
default is to use model terms.}

\item{group_names}{vector of names for each group of fixed effect. Should
match \code{group}. If the model was defined with \code{XData} and \code{XFormula}, the
default is to use the labels of model terms.}

\item{start}{index of first MCMC sample included. Default: \code{1L}.}

\item{na.ignore}{Logical. If \code{TRUE}, covariates are ignored for sites where
the focal species is NA when computing variance-covariance matrices for
each species.}

\item{n_cores}{Integer. Number of CPU cores to use for computing variance
partitioning using \code{TensorFlow}. This is only effective when \code{use_TF} is
\code{TRUE}. Default: \code{1}.}

\item{use_TF}{Logical. Whether to use \code{TensorFlow} for calculations. Defaults
to \code{TRUE}.}

\item{TF_environ}{Character. Path to the Python environment. This argument is
required if \code{use_TF} is \code{TRUE} under Windows. Defaults to \code{NULL}.}

\item{TF_use_single}{Logical. Whether to use single precision for the
\code{TensorFlow} calculations. Defaults to \code{FALSE}.}

\item{temp_cleanup}{Logical. Whether to delete temporary files after
processing. Default: \code{TRUE}.}

\item{chunk_size}{Integer. Size of each chunk of samples to process in
parallel. Only relevant for \code{TensorFlow}. Default: \code{50}.}

\item{verbose}{Logical. Whether to print progress messages. Default: \code{TRUE}.}

\item{VP_file}{Character. Name of the output file to save the results.
Default: \code{VarPar}.}

\item{VP_commands_only}{Logical. If \code{TRUE}, returns the commands to run the
Python script. Default is \code{FALSE}. Only relevant when \code{use_TF} is \code{TRUE}.}

\item{env_file}{Character. Path to the environment file containing paths to
data sources. Defaults to \code{.env}.}

\item{strategy}{Character. The parallel processing strategy to use. Valid
options are "sequential", "multisession" (default), "multicore", and
"cluster". See \code{\link[future:plan]{future::plan()}} and \code{\link[ecokit:set_parallel]{ecokit::set_parallel()}} for details.}

\item{width, height}{Numeric. Width and height of the output plot in
centimetres. Default: \code{30} and \code{15}, respectively.}

\item{Axis_text}{Numeric. Size of the axis text. Default: \code{4}.}
}
\description{
The \strong{\code{variance_partitioning_compute()}} function computes variance
components with respect to given grouping of fixed effects and levels of
random effects. This function inherits the main functionality from the
\code{Hmsc::computeVariancePartitioning} function, but with the added
functionality of parallel computation and using \code{TensorFlow}.\if{html}{\out{<br/>}}The
\strong{\code{variance_partitioning_plot()}} function generates plots for variance
partitioning as JPEG files. It allows for sorting the predictors and species;
e.g., by the mean value per predictor; and by original species order. It also
plots the raw variance partitioning (relative variance partitioning
multiplied by the Tjur-R\if{html}{\out{<sup>}}2\if{html}{\out{</sup>}} value).
}
\author{
Ahmed El-Gabbas
}
