#' Run blockFinder
#' 
#' @details This function has been deprecated and will be defunct in the new version.
#' 
#' @export
#' 
#' @param set \code{GenomicRatioSet}, \code{eSet} derived object or 
#' \code{SummarizedExperiment}
#' @param model Model matrix or formula to get model matrix from \code{set}. 
#' @param coefficient Numeric with the column of model matrix used in the analysis.
#' (Default: 2)
#' @param blockfinder_cutoff Numeric with the minimum cutoff to include a probe
#' in a block. (Default: 0.1)
#' @param num_permutations Numeric with the number of permutations run to compute 
#' the blocks p-value. (Default: 0)
#' @param verbose Logical value. Should the function be verbose? (Default: FALSE)
#' @param resultSet Should results be encapsulated in a \code{resultSet}? (Default: TRUE)
#' @param ... Further arguments passed to \code{blockFinder}.
#' @return data.frame or \code{resultSet} with the result of \code{blockFinder}
#' @seealso \code{\link[minfi]{blockFinder}}

runBlockFinder <- function(set, model, coefficient = 2, blockfinder_cutoff = 0.1, 
                           num_permutations = 0, resultSet = FALSE, 
                           verbose = FALSE, ...) {
  
  .Deprecated()
}