#' @describeIn AnalysisResults Get BlockFinder analysis results
#' @aliases AnalysisResults-methods blocks
#' @param object \code{AnalysisResults}
setMethod(
  f = "blocks",
  signature = "AnalysisResults",
  definition = function(object) {
    return(object@BlockFinder)
  }
)

#' @describeIn AnalysisResults Get Bumphunter analysis results
#' @aliases AnalysisResults-methods bumps
setMethod(
  f = "bumps",
  signature = "AnalysisResults",
  definition = function(object) {
    return(object@Bumphunter)
  }
)

#' @describeIn AnalysisResults Get covariable names
#' @aliases AnalysisResults-methods covariableNames
setMethod(
  f = "covariableNames",
  signature = "AnalysisResults",
  definition = function(object) {
    return(object@covariableNames)
  }
)

#' @describeIn AnalysisResults Get dmrCate analysis results
#' @aliases AnalysisResults-methods dmrCate
setMethod(
  f = "dmrCate",
  signature = "AnalysisResults",
  definition = function(object) {
    return(object@DMRcate)
  }
)

#' @describeIn AnalysisResults Get features names
#' @aliases AnalysisResults-methods feats
setMethod(
  f = "feats",
  signature = "AnalysisResults",
  definition = function(object) {
    return(rownames(featvals(object)))
  }
)

#' @describeIn AnalysisResults Get features values matrix
#' @aliases AnalysisResults-methods featvals
setMethod(
  f = "featvals",
  signature = "AnalysisResults",
  definition = function(object) {
    return(object@features)
  }
)

#' @describeIn AnalysisResults Get probe results of a gene
#' @aliases AnalysisResults-methods 
#' @param gene Character with the name of the gene
setMethod(
  f = "getGeneVals",
  signature = "AnalysisResults",
  definition = function(object, gene){
    if (!is(gene, "character")){
      stop("gene must be a character vector")
    }
    if (!"genes" %in% colnames(probeResults(object)[[1]])){
      stop("Results must have a column called genes to extract the probes.")
    }
    res <- lapply(probeResults(object), function(x){
      genelist <- strsplit(x$genes, ";")
      mask <- vapply(genelist, function(y) gene %in% y, logical(length(gene)))
      if (length(gene) > 1){
        mask <- as.logical(colSums(mask))
      }
      return(x[mask, ])
    }) 
    if (nrow(res[[1]]) == 0){
      warning("Gene name is not present in the results. A list with empty data.frames will be returned.")
    }
    res
  }
)


#' @describeIn AnalysisResults Get Ms values
#' @aliases AnalysisResults-methods 
#' @param threshold Numeric with the threshold to avoid 0s and 1s. 
setMethod(
  f = "getMs",
  signature = "AnalysisResults",
  definition = function(object, threshold = 0.0001){
    M <- featvals(object)
    M[M==0] <- threshold
    M[M==1] <- 1 - threshold
    M <- minfi::logit2(M)
  }
)

#' @describeIn AnalysisResults Get model used to perform the analysis
#' @aliases AnalysisResults-methods model
setMethod(
  f = "model",
  signature = "AnalysisResults",
  definition = function(object) {
    return(object@model)
  }
)


#' @describeIn AnalysisResults Get names of the variables in the model matrix
#' @aliases AnalysisResults-methods modelVariables
setMethod(
  f = "modelVariables",
  signature = "AnalysisResults",
  definition = function(object) {
    return(names(object@results))
  }
)

#' @describeIn AnalysisResults Get phenotypes data (AnnotatedDataFrame)
#' @aliases AnalysisResults-methods 
setMethod(
  f = "phenoData",
  signature = "AnalysisResults",
  definition = function(object) {
    return(object@phenotypes)
  }
)

#' @describeIn AnalysisResults Set phenotypes data (AnnotatedDataFrame)
#' @aliases AnalysisResults-methods  
#' @param value AnnotatedDataFrame or data.frame with the phenotype
setMethod(
  f = "phenoData<-",
  signature = "AnalysisResults",
  definition = function(object, value) {
    object@phenotypes <- value
    return(object)
  }
)

#' @describeIn AnalysisResults Get phenotypes data (data.frame)
#' @aliases AnalysisResults-methods 
setMethod(
  f = "pData",
  signature = "AnalysisResults",
  definition = function(object) {
    return(pData(object@phenotypes))
  }
)

#' @describeIn AnalysisResults Set phenotypes data (data.frame)
#' @aliases AnalysisResults-methods
setMethod(
  f = "pData<-",
  signature = "AnalysisResults",
  definition = function(object, value) {
    pData(object@phenotypes) <- value
    return(object)
  }
)

#' @describeIn AnalysisResults Get per probe analysis results
#' @aliases AnalysisResults-methods probeResults 
setMethod(
  f = "probeResults",
  signature = "AnalysisResults",
  definition = function(object) {
    return(object@results)
  }
)

#' @describeIn AnalysisResults Get all per region analysis results
#' @aliases AnalysisResults-methods regionResults
setMethod(
  f = "regionResults",
  signature = "AnalysisResults",
  definition = function(object) {
    results <- list()
    results[["DMRCate"]] <- dmrCate(object)
    results[["Bumphunter"]] <- bumps(object)
    results[["BlockFinder"]] <- blocks(object)
    return(results)
  }
)

#' @describeIn AnalysisResults Get sample names
#' @aliases AnalysisResults-methods 
setMethod(
  f = "sampleNames",
  signature = "AnalysisResults",
  definition = function(object) {
    return(object@sampleNames)
  }
)

#' @describeIn AnalysisResults Get variable names
#' @aliases AnalysisResults-methods variableNames
setMethod(
  f = "variableNames",
  signature = "AnalysisResults",
  definition = function(object) {
    return(object@variableNames)
  }
)