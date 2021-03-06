% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runDMRcate.R
\name{runDMRcate}
\alias{runDMRcate}
\title{Run DMRcate}
\usage{
runDMRcate(set, model, coefficient = 2, resultSet = FALSE, ...)
}
\arguments{
\item{set}{\code{GenomicRatioSet}, \code{eSet} derived object or 
\code{SummarizedExperiment}}

\item{model}{Model matrix or formula to get model matrix from \code{set}.}

\item{coefficient}{Numeric with the column of model matrix used in the analysis.
(Default: 2)}

\item{resultSet}{Should results be encapsulated in a \code{resultSet}? (Default: TRUE)}

\item{...}{Further arguments passed to \code{cpg.annotate} or \code{dmrcate}.}
}
\value{
data.frame or \code{resultSet} with the result of \code{bumphunter}
}
\description{
Run DMRcate
}
\details{
This function has been deprecated and will be defunct in the new version.
}
\seealso{
\code{\link[DMRcate]{dmrcate}}, \code{\link[DMRcate]{cpg.annotate}}
}
