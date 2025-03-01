% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/General_SplitVector.R
\name{SplitVector}
\alias{SplitVector}
\title{Split a vector into smaller chunks}
\usage{
SplitVector(Vector = NULL, NSplit = NULL, Prefix = "Chunk")
}
\arguments{
\item{Vector}{A numeric or character vector that you want to split.}

\item{NSplit}{Integer. Number of chunks to split the vector into. It must not
exceed the length of the vector.}

\item{Prefix}{Character. Prefix for the names of the chunks in the returned
list. Defaults to \code{"Chunk"}.}
}
\value{
A list of vectors, where each vector represents a chunk of the
original vector. The names of the list elements are generated using the
specified prefix followed by an underscore and the chunk number.
}
\description{
This function divides a given vector into a specified number of smaller
chunks. It is useful for partitioning data into more manageable pieces or for
parallel processing tasks.
}
\examples{
SplitVector(Vector = seq_len(100), NSplit = 3)

# -------------------------------------------

SplitVector(Vector = seq_len(100), NSplit = 2, Prefix = "T")

}
\author{
Ahmed El-Gabbas
}
