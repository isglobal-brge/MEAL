#' @describeIn MultiDataSet Method to add a slot of SNPs to \code{MultiDataSet}.
#' @aliases MultiDataSet-methods
#' @param snpSet \code{SnpSet} to be used to fill the slot.
setMethod(
  f = "add.snps",
  signature = c("MultiDataSet", "SnpSet"),
  definition = function(object, snpSet, warnings = TRUE) {
    if (!all(c("chromosome", "position") %in% fvarLabels(snpSet))){
      stop("snpSet must contain a fData with columns chromosome and position.")
    }
    object <- add.set(object, snpSet, "snps")
    
    return(object)
  }
)
