#' @describeIn MultiDataSet Method to add a slot to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param object \code{MultiDataSet}
#' @param set Object derived from \code{eSet} to be used to fill the slot.
#' @param dataset.name Character with the name of the slot to be filled.
#' @param warnings Logical to indicate if warnings will be displayed.
#' @return \code{MultiDataSet}
setMethod(
  f = "add.set",
  signature = c("MultiDataSet", "eSet"),
  definition = function(object, set, dataset.name, warnings = TRUE) {
    validObject(set)
    if (!is.null(sampleNames(object))){
      commonSamples <- intersect(sampleNames(set), sampleNames(object))
      if (length(commonSamples) == 0){
        stop("There are no common sample between the sets of the MultiDataSet and the new set.")
      }
      if (length(commonSamples) < length(sampleNames(object)) | length(commonSamples) < length(sampleNames(set))){
        warning("The samples of multiDataSet are not the same that samples in the new set. Only common samples will be retained.")
        object <- object[ , commonSamples]
        set <- set[ , commonSamples]
      }
    }
    if(dataset.name %in% names(object) & warnings) {
      warning("Slot '", dataset.name, "' is already set in 'MultiDataSet'. Previous content will be overwritten.")
    }
    
    object@assayData[[dataset.name]] <- assayData(set)
    object@phenoData[[dataset.name]] <- phenoData(set)
    object@featureData[[dataset.name]] <- featureData(set)
    
    returnfunc <- function(env, phe, fet) {
      new(class(set), assayData = env, phenoData = phe, featureData = fet)
    }
    
    object@return_method[[dataset.name]] <- returnfunc
    return(object)
  }
)
