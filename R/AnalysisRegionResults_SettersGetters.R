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


#' @describeIn AnalysisRegionResults Get R2 values of cpgs vs variables.
#' @aliases AnalysisRegionResults-methods regionLM
setMethod(
  f = "regionLM",
  signature = "AnalysisRegionResults",
  definition = function(object) {
    return(object@regionLM)
  }
)

#' @describeIn AnalysisRegionResults Get p-value of lineal model R2. 
#' @aliases AnalysisRegionResults-methods regionPval
setMethod(
  f = "regionPval",
  signature = "AnalysisRegionResults",
  definition = function(object) {
    return(object@regionPval)
  }
)

#' @describeIn AnalysisRegionResults Get R2 of the region vs variables lineal model 
#' @aliases AnalysisRegionResults-methods regionR2
setMethod(
  f = "regionR2",
  signature = "AnalysisRegionResults",
  definition = function(object) {
    return(object@regionR2)
  }
)

#' @describeIn AnalysisRegionResults Get SNPs data
#' @aliases AnalysisRegionResults-methods snps
setMethod(
  f = "snps",
  signature = "AnalysisRegionResults",
  definition = function(object) {
    return(object@snps)
  }
)

#' @describeIn AnalysisRegionResults Get p-values of correlations of snps-cpgs
#' pairs
#' @aliases AnalysisRegionResults-methods snpsPvals
setMethod(
  f = "snpsPvals",
  signature = "AnalysisRegionResults",
  definition = function(object) {
    return(object@snpsPvals)
  }
)

#' @describeIn AnalysisRegionResults Get variance of SNP matrix present in the 
#' component used to adjusting.
#' @aliases AnalysisRegionResults-methods snpsVar
setMethod(
  f = "snpsVar",
  signature = "AnalysisRegionResults",
  definition = function(object) {
    return(object@snpsVar)
  }
)