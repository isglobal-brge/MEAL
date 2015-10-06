#' @describeIn MultiDataSet Method to add a slot of methylation to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param methySet \code{MethylationSet} to be used to fill the slot.
setMethod(
  f = "add.methy",
  signature = c("MultiDataSet", "MethylationSet"),
  definition = function(object, methySet, warnings = TRUE) {
    
    if (!all(c("chromosome", "position") %in% fvarLabels(methySet))){
      stop("methyset must contain a fData with columns chromosome and position.")
    }
    
    object <- add.set(object, methySet, "methylation")
    
    return(object)
  }
)
