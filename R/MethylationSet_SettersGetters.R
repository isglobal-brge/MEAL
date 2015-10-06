#' @describeIn MethylationSet Get beta matrix
#' @aliases MethylationSet-methods betas
#' @param object \code{MethylationSet} 
setMethod(
  f = "betas",
  signature = "MethylationSet",
  definition = function(object) {
    return(assayDataElement(object, "meth"))
  }
)

#' @describeIn MethylationSet Get Ms values
#' @param threshold Numeric with the threshold to avoid 0s and 1s. 
#' @aliases MethylationSet-methods 
setMethod(
  f = "getMs",
  signature = "MethylationSet",
  definition = function(object, threshold = 0.0001){
    M <- betas(object)
    M[M==0] <- threshold
    M[M==1] <- 1 - threshold
    M <- minfi::logit2(M)
  }
)