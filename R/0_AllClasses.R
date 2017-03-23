#' AnalysisResults instances
#' 
#' Container with the results of per probe and per region analyses.
#' 
#' @export
#' @rdname AnalysisResults-class
#' @name AnalysisResults
#' @aliases AnalysisResults-class AnalysisResults-methods
#' 
#' @slot originalclass Character with the class of the object used to perform the 
#' analysis
#' @slot features Matrix with the values of the most significant features.
#' @slot phenotypes AnnotatedDataFrame with the phenotypes.
#' @slot model Matrix with the model used in the analysis
#' @slot sampleNames Character vector with the names of the samples
#' @slot variableNames Character vector with the names of the variables used in the
#' analysis. Names are equal to those find in phenotypes.
#' @slot covariableNames Character vector with the names of the covariables used in the
#' analysis. Names are equal to those find in phenotypes.
#' @slot results List of data.frames with the results of per probe analysis. Names are
#' those of the model.
#' @slot DMRcate List of data.frames with the results of \code{DMRcate}. Names are
#' those of the model.
#' @slot Bumphunter List of data.frames with the results of \code{Bumphunter}. Names are
#' those of the model.
#' @slot BlockFinder List of data.frames with the results of \code{BlockFinder}. Names are
#' those of the model.
setClass (
  Class = "AnalysisResults",
  representation(
    originalclass = "character", 
    features = "matrix",
    phenotypes = "AnnotatedDataFrame",
    model = "matrix",
    sampleNames = "character",
    variableNames = "character",
    covariableNames = "character",
    results = "list",
    DMRcate = "list",
    Bumphunter = "list",
    BlockFinder = "list"
  )
)

#' AnalysisRegionResults instances
#' 
#' \code{AnalysisResults} heir with the analyses performed in a range of the 
#' whole genome.
#'  
#' @export
#' @rdname AnalysisRegionResults-class
#' @name AnalysisRegionResults
#' @aliases AnalysisRegionResults-class AnalysisRegionResults-methods
#'  
#' @slot range \code{GenomicRanges} used to perform the analysis.
#' @slot rda \code{rda} object from \code{vegan} package with the results of RDA
#' analysis in the range.
#' @slot regionR2 Numeric with the R2 of the region calculated using a redundancy analysis.
#' @slot RDAPval Numeric with the p-value of the RDA. 
#' @slot globalR2 Numeric with the global R2.
#' @slot globalPval Numeric with the probability of finding a region with the same number of probes with
#' a bigger R2.
setClass (
  Class = "AnalysisRegionResults",
  contains = "AnalysisResults",
  representation(
    range = className("GRanges","GenomicRanges"),
    rda = "list",
    regionR2 = "numeric",
    RDAPval = "numeric",
    globalR2 = "numeric",
    globalPval = "numeric"
  )
)