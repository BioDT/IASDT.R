% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Mod_Parameter_Heatmap.R
\name{Parameter_Heatmap}
\alias{Parameter_Heatmap}
\alias{Mod_Heatmap_Beta}
\alias{Mod_Heatmap_Omega}
\title{Heatmaps for the \code{beta} and \code{omega} parameters of the Hmsc model}
\usage{
Mod_Heatmap_Beta(
  Path_Model = NULL,
  SupportLevel = 0.95,
  PlotWidth = 25,
  PlotHeight = 35
)

Mod_Heatmap_Omega(
  Path_Model,
  SupportLevel = 0.95,
  PlotWidth = 26,
  PlotHeight = 22.5
)
}
\arguments{
\item{Path_Model}{Character. Path to the fitted Hmsc model object.}

\item{SupportLevel}{Numeric. The posterior support threshold for determining
which values are considered significant in the heatmap. Defaults to 0.95,
indicating 95\% posterior support. Values above this threshold (or below 1 -
threshold for negative associations) are considered significant and will be
plotted (see \link[Hmsc:plotBeta]{Hmsc::plotBeta}).}

\item{PlotWidth, PlotHeight}{Integer. The width and height of the generated
heatmaps in centimeters. Defaults to 26×22.5 for \code{omega}; 25×35
for \code{beta}.}
}
\value{
Both functions do not return a value but saves heatmap plots as JPEG
files in the \code{Model_Postprocessing/Parameters_Summary} subdirectory.
}
\description{
The \code{Mod_Heatmap_Beta()} and \code{Mod_Heatmap_Omega()} functions generate
heatmaps using \code{ggplot2} to visualize parameter estimates or posterior
support values for species' environmental responses (\code{beta} parameters, which
describes how species (\emph{Y}) respond to various covariates (\emph{X}); see
\link[Hmsc:plotBeta]{Hmsc::plotBeta}) and residual associations (\code{omega} parameter),
respectively.
}
\details{
The functions exports three types of visualizations (see
\link[Hmsc:plotBeta]{Hmsc::plotBeta}):
\itemize{
\item \code{Mean}: posterior mean estimate,
\item \code{Support}: statistical support level, measured by the posterior
probability for a positive or negative response,
\item \code{Sign}: indicates whether the response is positive, negative, or neither
of these based on the chosen \code{SupportLevel}.
}

For the \code{omega} parameter, the \code{Mod_Heatmap_Omega()} function generates two
JPEG files: signs and mean values. While for the \code{beta} parameter, the
\code{Mod_Heatmap_Beta()} function generates four JPEG files : support, signs,
mean values (including and excluding the intercept).
}
\author{
Ahmed El-Gabbas. The \code{Mod_Heatmap_Beta()} function is adapted from
\link[Hmsc:plotBeta]{Hmsc::plotBeta}
}
