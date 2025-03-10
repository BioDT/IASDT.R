% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mod_Predict_Plot.R
\name{Mod_Predict_Plot}
\alias{Mod_Predict_Plot}
\title{Plot species and level of invasion predictions as JPEG files using \code{ggplot2}}
\usage{
Mod_Predict_Plot(ModelDir = NULL, EnvFile = ".env", NCores = 8L)
}
\arguments{
\item{ModelDir}{Character. Path to the model directory containing
predictions.}

\item{EnvFile}{Character. Path to the environment file containing paths to
data sources. Defaults to \code{.env}.}

\item{NCores}{Integer. Number of CPU cores to use for parallel processing.
Default: 8.}
}
\value{
Saves prediction plots as JPEG files in the specified output
directory.
}
\description{
Generate predictions for species and habitat models and saves the output as
JPEG files.
}
\author{
Ahmed El-Gabbas
}
