#' @describeIn AnalysisRegionResults Get range where the analyses was performed
#' @aliases AnalysisRegionResults-methods getRange
#' @param object \code{MethylationResults} 
setMethod(
  f = "getRange",
  signature = "AnalysisRegionResults",
  definition = function(object) {
    return(object@range)
  }
)

#' @describeIn AnalysisRegionResults Get rda object.
#' @aliases AnalysisRegionResults-methods getRDA
setMethod(
  f = "getRDA",
  signature = "AnalysisRegionResults",
  definition = function(object) {
    return(object@rda)
  }
)

#' @describeIn AnalysisRegionResults Get global p-value.
#' @aliases AnalysisRegionResults-methods globalPval
setMethod(
  f = "globalPval",
  signature = "AnalysisRegionResults",
  definition = function(object) {
    return(object@globalPval)
  }
)

#' @describeIn AnalysisRegionResults Get global R2.
#' @aliases AnalysisRegionResults-methods globalR2
setMethod(
  f = "globalR2",
  signature = "AnalysisRegionResults",
  definition = function(object) {
    return(object@globalR2)
  }
)

#' @describeIn AnalysisRegionResults Get p-value of RDA. 
#' @aliases AnalysisRegionResults-methods RDAPval
setMethod(
  f = "RDAPval",
  signature = "AnalysisRegionResults",
  definition = function(object) {
    return(object@RDAPval)
  }
)

#' @describeIn AnalysisRegionResults Get R2 of the RDA model
#' @aliases AnalysisRegionResults-methods regionR2
setMethod(
  f = "regionR2",
  signature = "AnalysisRegionResults",
  definition = function(object) {
    return(object@regionR2)
  }
)