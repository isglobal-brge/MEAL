#' Get a summary of RDA results
#' 
#' Get statistics from RDA result.
#'  
#' @param object \code{ResultSet}
#' @name getRDAresults
#' @aliases getRDAresults
#' @export

getRDAresults <- function(object) {
    
  stopifnot(is(object, "ResultSet"), "RDA" %in% names(object))
    rda <- getAssociation(object, rid = "RDA")
    res <- c(R2 = rda$rdaR2, pval = rda$pval, global.R2 = unname(rda$globalR2), 
                global.pval = unname(rda$globalpval))
    
    return(res)
  }
