#' @describeIn MultiDataSet Get sample names
#' @aliases MultiDataSet-methods
setMethod(
  f = "sampleNames",
  signature = "MultiDataSet",
  definition = function(object) {
    tables <- names(object)
    if (is.null(tables)){
      return(NULL)
    }
    return(rownames(object@phenoData[[1]]))
  }
)