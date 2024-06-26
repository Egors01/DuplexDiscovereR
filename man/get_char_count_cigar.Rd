% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_general.R
\name{get_char_count_cigar}
\alias{get_char_count_cigar}
\title{Count the length of the key type in CIGAR string}
\usage{
get_char_count_cigar(strings, s)
}
\arguments{
\item{strings}{CIGAR string vector}

\item{s}{CIGAR operands}
}
\value{
vector with length values
}
\description{
Takes CIGAR operands i.e M,N,S and sums the associated blocks length
It is vectorized. i.e supports vector with CIGAR strings
}
\examples{
# From a vector
get_char_count_cigar(c("4S18M22S", "25S26M"), "S")
get_char_count_cigar(c("18M22S", "20M20S"), "M")
}
