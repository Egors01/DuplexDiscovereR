% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_preproc.R
\name{preproc_chim_junction_out_se}
\alias{preproc_chim_junction_out_se}
\title{Processing of of the STAR SE Chimeric.junction.out}
\usage{
preproc_chim_junction_out_se(dt, keep_all_columns = FALSE)
}
\arguments{
\item{dt}{Chimeric.out.junction with the correct column names}

\item{keep_all_columns}{\itemize{
\item TRUE or FALSE. Keep CIGAR strings and junction coordinate columns
}}
}
\value{
tibble with annotated reads
}
\description{
Calculates alignment coordinates and returns reads with categories
}
\details{
#' \describe{
\item{multimap}{ multi-mapped read}
\item{multigap}{ more than one junction (more than two 'N' in CIGAR string)}
\item{bad junction}{ Artifacts. I.e alignments for both arms are continious, but with 'backward' chimeric junction was wrongly put}
}
}
\seealso{
\code{\link[=col_check_rename]{col_check_rename()}}
}
\keyword{interal}
