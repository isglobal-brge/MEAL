setMethod(
  f = "initialize",
  signature = "MethylationSet",
  definition = function(.Object, assayData, ...) {
    if (missing(assayData)){
      assayData <- assayDataNew(storage.mode = "lockedEnvironment",  
                                meth = matrix(ncol = 0, nrow = 0))
    } 
    .Object <- callNextMethod(.Object, assayData = assayData, ...)
    return(.Object)
  })