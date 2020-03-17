#' Run DMRcate
#' 
#' @details This function has been deprecated and will be defunct in the new version.
#' @export
#' 
#' @param set \code{GenomicRatioSet}, \code{eSet} derived object or 
#' \code{SummarizedExperiment}
#' @param model Model matrix or formula to get model matrix from \code{set}. 
#' @param coefficient Numeric with the column of model matrix used in the analysis.
#' (Default: 2)
#' @param resultSet Should results be encapsulated in a \code{resultSet}? (Default: TRUE)
#' @param ... Further arguments passed to \code{cpg.annotate} or \code{dmrcate}.
#' @return data.frame or \code{resultSet} with the result of \code{bumphunter}
#' @seealso \code{\link[DMRcate]{dmrcate}}, \code{\link[DMRcate]{cpg.annotate}}
runDMRcate <- function(set, model, coefficient = 2, resultSet = FALSE, ...){
  
  .Deprecated()
}