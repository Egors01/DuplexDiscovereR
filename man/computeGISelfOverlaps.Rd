% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_clustering.R
\name{computeGISelfOverlaps}
\alias{computeGISelfOverlaps}
\title{Find overlaps between entries in \code{GInteractions}}
\usage{
computeGISelfOverlaps(
  gi,
  id_column = "duplex_id",
  maxgap = 40,
  minoverlap = 10
)
}
\arguments{
\item{gi}{input gi object}

\item{id_column}{column which use for using as ids for entries}

\item{maxgap}{parameter for call of \code{\link[InteractionSet:overlaps]{InteractionSet::findOverlaps()}}}

\item{minoverlap}{parameter for call \code{\link[InteractionSet:overlaps]{InteractionSet::findOverlaps()}}}
}
\value{
dataframe with indexes of pairwise overlapsin input and
columns for span, overlap, ratios of either read arm
}
\description{
Utility function to find overlapping reads in the input and calculate
overlap scores. Removes self-hits. Computes overlap/span ratios for
each interaction arm. Sum of the scores is recorded in 'weight' field
}
\examples{
data("RNADuplexesSmallGI")
computeGISelfOverlaps(SampleSmallGI)
}
