#' Get a summary of RDA results
#' 
#' Get statistics from RDA result.
#'  
#'  
#' @aliases getRDAresults
#' @export
#' @name getRDAresults
#' @param object \code{ResultSet}
#' @return Numeric vector with the RDA statistics
getRDAresults <- function(object) {
    
  stopifnot(is(object, "ResultSet"), "RDA" %in% names(object))
    rda <- getAssociation(object, rid = "RDA")
    res <- c(R2 = rda$rdaR2, pval = rda$pval, global.R2 = unname(rda$globalR2), 
                global.pval = unname(rda$globalpval))
    
    return(res)
  }
