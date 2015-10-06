#' @describeIn MethylationSet  Modify a \code{MethylationSet} to 
#' only contain common samples
#' @aliases MethylationSet-methods
setMethod(
  f = "checkSamples",
  signature = "MethylationSet",
  definition = function(object) {
    bothsamples <- intersect(colnames(betas(object)), rownames(pData(object)))
    if(!length(bothsamples)){
      stop("There are no common samples between the beta matrix and the phenotypes table. Please, check sample names.")
    }      
    
    object <- object[, bothsamples]
    phenoData(object) <- phenoData(object)[bothsamples, , drop=FALSE]
    sampleNames(object) <- bothsamples
    return(object)
  }
)
