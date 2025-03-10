% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mod_VarPar_Compute.R, R/Mod_VarPar_Plot.R
\name{Variance_partitioning}
\alias{Variance_partitioning}
\alias{VarPar_Compute}
\alias{VarPar_Plot}
\title{Computes and visualize variance partitioning of Hmsc models}
\usage{
VarPar_Compute(
  Path_Model,
  group = NULL,
  groupnames = NULL,
  start = 1L,
  na.ignore = FALSE,
  NCores = 8L,
  UseTF = TRUE,
  TF_Environ = NULL,
  TF_use_single = FALSE,
  Temp_Cleanup = TRUE,
  Chunk_size = 50L,
  Verbose = TRUE,
  VarParFile = "VarPar",
  VP_Commands_Only = FALSE
)

VarPar_Plot(
  Path_Model,
  EnvFile = ".env",
  VarParFile = "VarPar",
  UseTF = TRUE,
  TF_Environ = NULL,
  NCores = 1,
  Fig_width = 30,
  Fig_height = 15,
  Axis_text = 4
)
}
\arguments{
\item{Path_Model}{Character. Path to fitted \code{Hmsc} model object.}

\item{group}{vector of numeric values corresponding to group identifiers in
groupnames. If the model was defined with \code{XData} and \code{XFormula}, the
default is to use model terms.}

\item{groupnames}{vector of names for each group of fixed effect. Should
match \code{group}. If the model was defined with \code{XData} and \code{XFormula}, the
default is to use the labels of model terms.}

\item{start}{index of first MCMC sample included}

\item{na.ignore}{Logical. If \code{TRUE}, covariates are ignored for sites where
the focal species is NA when computing variance-covariance matrices for
each species.}

\item{NCores}{Integer. Number of CPU cores to use for computing variance
partitioning using TensorFlow. This is only effective when \code{UseTF} is
\code{TRUE}. Default: \code{1}.}

\item{UseTF}{Logical. Whether to use TensorFlow for calculations. Defaults to
\code{TRUE}.}

\item{TF_Environ}{Character. Path to the Python environment. This argument is
required if \code{UseTF} is \code{TRUE} under Windows. Defaults to \code{NULL}.}

\item{TF_use_single}{Logical. Whether to use single precision for the
TensorFlow calculations. Defaults to \code{FALSE}.}

\item{Temp_Cleanup}{Logical. Whether to delete temporary files after
processing. Default: \code{TRUE}.}

\item{Chunk_size}{Integer. Size of each chunk of samples to process in
parallel. Only relevant for TensorFlow. Default: \code{50}.}

\item{Verbose}{Logical. Whether to print progress messages. Default: \code{TRUE}.}

\item{VarParFile}{Character. Name of the output file to save the results.
Default: \code{VarPar}.}

\item{VP_Commands_Only}{Logical. If \code{TRUE}, returns the commands to run the
Python script. Default is \code{FALSE}. Only relevant when \code{UseTF} is \code{TRUE}.}

\item{EnvFile}{Character. Path to the environment file containing paths to
data sources. Defaults to \code{.env}.}

\item{Fig_width, Fig_height}{Numeric. Width and height of the output plot in
centimeters. Default: \code{30} and \code{15}, respectively.}

\item{Axis_text}{Numeric. Size of the axis text. Default: \code{4}.}
}
\description{
The \strong{\code{VarPar_Compute()}} function computes variance components with respect
to given grouping of fixed effects and levels of random effects. This
function inherits the main functionality from the
\code{Hmsc::computeVariancePartitioning} function, but with the added
functionality of parallel computation and using TensorFlow.\if{html}{\out{<br/>}}The
\strong{\code{VarPar_Plot()}} function generates plots for variance partitioning as
JPEG files. It allows for sorting the predictors and species; e.g., by the
mean value per predictor; and by original species order. It also plots the
raw variance partitioning (relative variance partitioning multiplied by the
Tjur-R\if{html}{\out{<sup>}}2\if{html}{\out{</sup>}} value).
}
\author{
Ahmed El-Gabbas
}
