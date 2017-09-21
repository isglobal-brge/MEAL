#' Run blockFinder
#' 
#' Run blockFinder to a methylation dataset. This function contains all steps of
#' blockFinder analysis, from model.matrix creation to running the analysis. 
#' 
#' @details runBlockFinder is a wrapper for minfi \code{blockFinder}. This function 
#' runs all the steps required prior running blockFinder from the methylation set
#' and the formula of the model. This implementation allows running blockFinder to other objects than 
#' \code{GenomicRatioSet}. The result can be encapsulated in a \code{ResultSet} to
#' take adavantege of its plotting capabilities.
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
  
  ## Create model matrix from formula
  if (is(model, "formula")){
    model <- createModel(set, model)
    set <- set[, rownames(model)]
  }
  
  
  ## Convert eSets and SummarizedExperiments to GenomicRatioSet
  if (is(set, "eSet")){
    granges <- GenomicRanges::makeGRangesFromDataFrame(fData(set), 
                                                       start.field = "position", end.field = "position")
    beta <- Biobase::exprs(set)
    set <- minfi::GenomicRatioSet(gr = granges, Beta = beta, colData = pData(set),
                                  CN = matrix(1, ncol = ncol(beta), nrow = nrow(beta)),
                                  annotation = c(array = "IlluminaHumanMethylation450k",
                                                 annotation = "ilmn12.hg19"))
  } else if (is(set, "SummarizedExperiment") & !is(set, "GenomicRatioSet")){
    set <- as(set, "RangedSummarizedExperiment")
    set <- minfi::GenomicRatioSet(gr = rowRanges(set), Beta = assays(set), 
                                  colData = colData(set),
                                  CN = matrix(1, ncol = assays(set), nrow = assays(set)),
                                  annotation = c(array = "IlluminaHumanMethylation450k",
                                                 annotation = "ilmn12.hg19"))
  } else if (!is(set, "GenomicRatioSet")){
    stop("set must be an eSet, GenomicRatioSet or SummarizedExperiment.")
  }
  set <- sort(set)
  
  ## Prepare object to run blockFinder
  blockset <- minfi::cpgCollapse(set, returnBlockInfo = FALSE, verbose = verbose )
  
  res <- tryCatch(minfi::blockFinder(blockset, design = model, coef = coefficient,
                                     cutoff = blockfinder_cutoff,
                                     B = num_permutations,
                                     nullMethod = "bootstrap",
                                     verbose = verbose, ...)$table, 
                  error = function(err) {
                    if(startsWith(as.character(err), "Error in res$table$indexStart")) {
                      warning("No blocks were found.")
                    } else {
                      stop(err)
                    }
                    return(data.frame())
                  })
  
  if (resultSet)
  {
    fFun <- getFeatureDataFun(set)
    
    res <- create_resultset("runBlockFinder", lResults = 
                              list(blockFinder = list(result = res, error = NA)),  
                            fData = list(main = fFun(set)), lOptions = list())
  }
  res
}