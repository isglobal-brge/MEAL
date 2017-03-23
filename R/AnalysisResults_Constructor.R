#' @export 
#' @name analysisResults
#' @rdname AnalysisResults-class
#' @aliases AnalysisResults-methods
#' 
#' @param set \code{MethylationSet} or \code{ExpressionSet} used to perform the analysis
#' @param model Model matrix used to produce the calculations
#' @param regionResults List with the region results
#' @param probeResults List with the probe results
#' @param num_feat Numeric with the minimum number of feature values to be included.
#' @param num_vars Numeric with the number of columns of the pData table 
#' that should be considered as variables.
#' @return \code{AnalysisResults}
#' @examples
#' showClass("AnalysisResults")
analysisResults <- function(set, model, regionResults, probeResults, num_feat = 50, 
                               num_vars = ncol(pData(set))){
  if (!all(names(regionResults) %in% c("bumphunter", "blockFinder", "DMRcate"))){
    varnames <- names(regionResults)
    regionResults$bumphunter <- lapply(varnames,
                                       function(x) regionResults[[x]][["bumphunter"]])
    names(regionResults$bumphunter) <- varnames
    regionResults$blockFinder <- lapply(varnames,
                                        function(x) regionResults[[x]][["blockFinder"]])
    names(regionResults$blockFinder) <- varnames
    regionResults$DMRcate <- lapply(varnames,
                                    function(x) regionResults[[x]][["DMRcate"]])
    names(regionResults$DMRcate) <- varnames
  }else{
    varname <- colnames(pData(set))[1]
    regionResults$bumphunter <- list(regionResults$bumphunter)
    names(regionResults$bumphunter) <- varname
    regionResults$blockFinder <- list(regionResults$blockFinder)
    names(regionResults$blockFinder) <- varname
    regionResults$DMRcate <- list(regionResults$DMRcate)
    names(regionResults$DMRcate) <- varname
  }
  maxfeat <- nrow(set)  
  if (is(probeResults, "list")){
    sigfeat <- lapply(names(probeResults),
                      function(x) rownames(probeResults[[x]])[probeResults[[x]]$adj.P.Val < 0.05])
    feats <- unique(unlist(sigfeat)) 
    num_feat <- min(maxfeat, max(num_feat, length(feats)))
    if (num_feat > length(feats)){
      feats <- c(feats, rownames(probeResults[[1]])[1:num_feat])
      feats <- feats[1:num_feat]
    }
  }else{
    sigfeat <- sum(probeResults$adj.P.Val < 0.05)
    num_feat <- min(maxfeat, max(num_feat, sigfeat))
    feats <- as.character(rownames(probeResults)[1:num_feat])
    probeResults <- list(probeResults)
    names(probeResults) <- colnames(pData(set))[1]
  }
  sampleNames <- colnames(set)
  variablesNames <- colnames(pData(set))[1:num_vars]
  if (ncol(pData(set)) > num_vars){
    covariableNames <- colnames(pData(set))[(num_vars+1):ncol(pData(set))]
  }else{
    covariableNames <- character()
  }
  phenotypes <- phenoData(set)
  if (is(set, "ExpressionSet")){
    features <- exprs(set[feats, ])
  }else{
    features <- MultiDataSet::betas(set[feats, ])
  }
  results <- new(Class = "AnalysisResults", features = features, phenotypes = phenotypes,
                 model = model, sampleNames = sampleNames, variableNames = variablesNames, 
                 covariableNames = covariableNames, results = probeResults, 
                 DMRcate = regionResults$DMRcate, Bumphunter = regionResults$bumphunter, 
                 BlockFinder = regionResults$blockFinder, originalclass = class(set)) 
  return(results)
}