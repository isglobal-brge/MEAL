#' @describeIn MultiDataSet Get an eSet from a slot
#' @aliases MultiDataSet-methods
#' @param x \code{MultiDataSet}
#' @param i slot
setMethod(
  f = "[[",
  signature = "MultiDataSet",
  definition = function(x, i) {
    if (i %in% names(x)) {
      set <- x@return_method[[i]](x@assayData[[i]], x@phenoData[[i]], x@featureData[[i]])
      validObject(set)
      return(set)
    }
  }
)

#' @describeIn MultiDataSet Subset a MultiDataSet
#' @aliases MultiDataSet-methods [
#' @param j samples
#' @param drop Logical indicating if dropped will be applied. 
setMethod(
  f = "[",
  signature = c(x = "MultiDataSet", i = "ANY"), 
  definition = function(x, i, j, drop = TRUE) {
    if (!missing(i)) {
      if (length(i) == 1 & drop){
        x <- x@return_method[[i]](x@assayData[[i]], x@phenoData[[i]], x@featureData[[i]])
        validObject(x)
        } else{
        x@return_method <- x@return_method[i]
        x@phenoData <- x@phenoData[i]
        x@featureData <- x@featureData[i]
        x@assayData <- x@assayData[i]
      }
    }
    if(!missing(j)) {
      if (is(x, "MultiDataSet")){
        sets <- lapply(names(x), function(y) x[y][ , j])
        names(sets) <- names(x)
        x <- new("MultiDataSet")
        for (i in 1:length(sets)){
          x <- add.set(x, sets[[i]], names(sets)[i])
        }
      }else{
        x <- x[ , j]
      }
    }
    return(x)
  }
)