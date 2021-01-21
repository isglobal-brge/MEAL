#' Run different DMR detection methods
#' 
#' @export runRegionAnalysis
#' @details This function has been deprecated and will be defunct in the new version.
#'  
#' 
#' @param set \code{GenomicRatioSet}, \code{eSet} derived object or 
#' \code{SummarizedExperiment}
#' @param model Model matrix representing a linear model.
#' @param methods Character vector with the names of the methods used to estimate 
#' the regions. Valid names are: "blockFinder", "bumphunter" and "DMRcate".
#' @param coefficient Numeric with the index of the model matrix used to perform
#' the analysis.
#' @param bumphunter_params List with other parameter passed to \code{runBumphunter} 
#' function.
#' @param blockFinder_params List with other parameter passed to \code{runBlockFinder} 
#' function.
#' @param dmrcate_params List with other parameter passed to \code{runDMRcate} 
#' function.
#' @param verbose Logical value. Should the function be verbose? (Default: FALSE)
#' @param resultSet Should results be encapsulated in a \code{resultSet}? (Default: TRUE)
#' @return List or \code{resultSet} with the result of the DMR detection methods.
#' @seealso \code{\link[minfi]{bumphunter}}, \code{\link[minfi]{blockFinder}}, 
#' \code{\link[DMRcate]{dmrcate}} 
#' @examples
#' if (require(minfiData)){
#' set <- ratioConvert(mapToGenome(MsetEx[1:10,]))
#'  model <- model.matrix(~Sample_Group, data = pData(MsetEx))
#'  res <- runRegionAnalysis(set, model) 
#'  res
#' }
runRegionAnalysis <- function (set, model, methods = c("blockFinder", "bumphunter", "DMRcate"), 
                            coefficient = 2, bumphunter_params = NULL, 
                            blockFinder_params = NULL, dmrcate_params = NULL,
                            verbose = FALSE, resultSet = TRUE)
{
  .Deprecated()
  
}