#' Obtain probe results from a \code{ResultSet}
#' 
#' It computes the statistics from the \code{MArrayLM} computed with 
#' \code{DiffMeanAnalysis} or \code{DiffVarAnalysis}. This function allows to 
#' specify the contrasts and to get F-statistics for a group of variables. 
#' 
#' @name getProbeResults
#' @aliases getProbeResults 
#' 
#' @param object ResultSet
#' @param rid Name of the results: "DiffMean" for mean differences, "DiffVar" for
#' variance differences. (Default: DiffMean)
#' @param coef Number of the coefficient used to compute the statistics. If a vector
#' is supplied, F-statistics evaluating the global effect of the coefficients are
#' computed. (Default: 2).
#' @param contrast Matrix of contrasts 
#' @param fNames Names of the columns of \code{object} fData that will be added to 
#' the results data.frame.
#' @param ... Further arguments passed to \code{getAssociation}.
#' @return data.frame with the probe results. 
#' @export
getProbeResults <- function(object, rid = "DiffMean", coef = 2, contrast = NULL, 
                        fNames = c("chromosome", "start"), robust = FALSE, ...) {
    
  stopifnot(is(object, "ResultSet"))
  
    if (is.numeric(rid)){
      rid <- names(object)[rid]
    }
    
    if (rid == "DiffMean") {
      res <- getAssociation(object, rid = rid, coef = coef, contrast = contrast, 
                            fNames = fNames, robust = robust, ...)
    } else if (rid == "DiffVar"){
      res <- object@results[[rid]]$result
      
      if (!is.null(contrast)){
        res <- missMethyl::contrasts.varFit(res, contrast)
      }
      res <- missMethyl::topVar(res, coef = coef, number = nrow(res))
      
      ## Add fData to results
      if (!is.null(fNames)){
        fData <- object@fData[[1]]
        
        if (!all(fNames %in% colnames(fData))){
          stop("All fNames must be present in ResultSet fData.")
        }
        
        res <- cbind(res, fData[rownames(res), fNames])  
      }
    } else {
      stop("This ResultSet does not contain probe results")
    }
    
    return(res)
  }
