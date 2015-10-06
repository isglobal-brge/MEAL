#' @describeIn MultiDataSet Get names of slots
#' @aliases MultiDataSet-methods
setMethod(
  f = "names",
  signature = "MultiDataSet",
  definition = function(x) {
    return(names(x@assayData))
  }
)
