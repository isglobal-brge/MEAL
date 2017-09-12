#' Run DMRcate
#' 
#' Run DMRcate to a methylation dataset. This function contains all steps of
#' DMRcate analysis, from model.matrix creation to running the analysis. 
#' 
#' @details runDMRcate is a wrapper for \code{dmrcate} function. runDMRcate  
#' runs all the steps required prior running blockFinder from the methylation set
#' and the formula of the model. This implementation allows running blockFinder to 
#' other objects than \code{GenomicRatioSet}. The result can be encapsulated in a 
#' \code{ResultSet} to take adavantege of its plotting capabilities.
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
  
  ## Create model matrix from formula
  if (is(model, "formula")){
    model <- createModel(set, model)
    set <- set[, rownames(model)]
  }
  
  
  ## Get matrix
  if (is(set, "GenomicRatioSet")){
    mat <- set
  } else if (is(set, "MethylationSet")){
    mat <- MultiDataSet::getMs(set)
  } else if (is(set, "SummarizedExperiment")){
    mat <- Biobase::assays(set)
  } else {
    stop("set must be a MethylationSet, GenomicRatioSet or SummarizedExperiment.")
  }
  
  annParams <- formals(DMRcate::cpg.annotate)
  annParams <- annParams[!names(annParams) %in% c("object", "design", "what", "coef")]
  
  dmrParams <- formals(DMRcate::dmrcate)
  dmrParams <- dmrParams[!names(dmrParams) %in% c("object")]
  
  dots <- list(...)

  annParams[intersect(names(dots), names(annParams))] <- 
    dots[intersect(names(dots), names(annParams))]
  
  dmrParams[intersect(names(dots), names(dmrParams))] <- 
    dots[intersect(names(dots), names(dmrParams))]
  
  myannotation <- do.call(DMRcate::cpg.annotate, 
                          c(list(object = mat, design = model, 
                                 what = "M", coef = coefficient), annParams))

  dmrcoutput <- do.call(DMRcate::dmrcate,
                        c(list(object = myannotation), dmrParams))
  
  res <- dmrcoutput$results
  if (resultSet)
  {
    fFun <- getFeatureDataFun(set)
    
    res <- create_resultset("runDMRcate", lResults = 
                              list(dmrcate = list(result = res, error = NA)),  
                            fData = list(main = fFun(set)), lOptions = list())
  }
  res
}