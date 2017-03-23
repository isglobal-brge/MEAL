#' @export
#' 
#' @name analysisRegionResults
#' @rdname AnalysisRegionResults-class
#' @aliases AnalysisRegionResults-methods
#' 
#' @param analysisResults \code{AnalysisResults}
#' @param range \code{GenomicRanges}
#' @param rdaRes List with RDA results
#' @return An \code{AnalysisRegionResults}
#' @examples
#' showClass("AnalysisRegionResults")
analysisRegionResults <- function(analysisResults, range, rdaRes){
  
  rdaobj <- rdaRes$rda
  class(rdaobj) <- c("list", class(rdaobj))
  results <- new(Class = "AnalysisRegionResults", analysisResults, range = range, 
                 rda = rdaobj, regionR2 = rdaRes$rdaR2, RDAPval = rdaRes$rdapval, 
                 globalR2 = rdaRes$globalR2, globalPval = rdaRes$globalpval)
  return(results)
}