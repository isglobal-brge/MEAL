#' @describeIn MultiDataSet Method to add a slot of expression to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param gexpSet \code{ExpressionSet} to be used to fill the slot.
setMethod(
  f = "add.genexp",
  signature = c("MultiDataSet", "ExpressionSet"),
  definition = function(object, gexpSet, warnings = TRUE) {
    if (!all(c("chromosome", "start", "end") %in% fvarLabels(gexpSet))){
      stop("gexpSet must contain a fData with columns chromosome, start and end.")
    }
    object <- add.set(object, gexpSet, "expression")
    
    return(object)
  }
)
