#' @describeIn AnalysisRegionResults Plot R2 region values
#' @aliases AnalysisRegionResults-methods
#' @param feat Numeric with the index of the cpg or character with its name.
#' @param ... Further arguments passed to plotLM
setMethod(
  f = "plotRegionR2",
  signature = "AnalysisRegionResults",
  definition = function(object, feat, ...){
  if (!length(snps(object))){
    warning("There are no relevant snps to any of the cpgs.")
  }
  regionVariance <- regionLM(object)
  if (!length(regionVariance)){
    stop("There is no information of R squared in the set.")
  }
  if (length(feat) != 1){
    stop("feat must contain only one value.")
  }
  if (is.numeric(feat)) {
    if (feat < 1 | feat > length(feats(object))){
      stop("feat index must be greater than 0 and smaller than the number of cpgs.")
    }
    feat <- feats(object)[feat]
  }
  if (!feat %in% feats(object)){
    stop("feat name is not present in the set.")
  }
  plotLM(regionVariance[[feat]], feat_name = feat, ...)
})