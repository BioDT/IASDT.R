% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod_trim_hmsc.R
\name{trim_hmsc}
\alias{trim_hmsc}
\title{Trim an Hmsc Model Object by Removing Specified Components}
\usage{
trim_hmsc(model, names_to_remove = NULL)
}
\arguments{
\item{model}{An object of class \code{Hmsc}, containing fitted Hmsc model. Must
not be \code{NULL}.}

\item{names_to_remove}{A character vector specifying the names of components
to remove from the model (e.g., \code{"postList"}, \code{"Y"}). If \code{NULL}, no
trimming is implemented. Must be non-empty and match names in the model.}
}
\value{
An \code{Hmsc} object with the specified components removed.
}
\description{
Removes specified components from an \code{Hmsc} model object one at a time. This
is used to keep only a smaller \code{Hmsc} object, containing only needed list
items. This is useful in situations in which the fitted Hmsc model is large
(e.g. in GBs) and downstream computations are implemented on parallel.
}
\details{
This function is used to reduce the memory footprint of \code{Hmsc}
models by removing unnecessary components before further processing. It
converts the model to a plain list to avoid S3 method overhead, removes
components iteratively, and restores the \code{Hmsc} class. The iterative
approach is slower than vectorized subsetting but may be preferred for
specific use cases requiring step-by-step removal. The simple trimming of
list items; e.g. \verb{model[c("postList", X")] <- NULL} took comparably too
much time for trimming done inside functions for jobs submitted via SLURM
}
\examples{
require(Hmsc)

(model <- Hmsc::TD$m)

(trimmed_model <- trim_hmsc(
  model, names_to_remove = c("postList", "rL", "ranLevels")))

setdiff(names(model), names(trimmed_model))

lobstr::obj_size(model)
lobstr::obj_size(trimmed_model)

}
\author{
Ahmed El-Gabbas
}
