#' @name methylationSet
#' @rdname MethylationSet-class
#' @aliases MethylationSet MethylationSet-methods
#' 
#' @param betas Matrix of beta values
#' @param phenotypes Data.frame or AnnotatedDataFrame with the phenotypes
#' @param annotationDataFrame Data.frame or AnnotatedDataFrame with the phenotypes with 
#' the annotation of the methylation sites. A column with the chromosomes named chr 
#' and a column with the positions names pos are required.
#' @param annoString Character with the name of the annotation used. 
#' @return \code{MethylationSet}
#' @examples
#' showClass("MethylationSet")
methylationSet <- function(betas, phenotypes, annotationDataFrame, annoString = "custom"){
  
  meth <- assayDataNew(storage.mode = "lockedEnvironment",  meth = betas)
  set <- new(Class = "MethylationSet", assayData = meth)
  
  if (is(phenotypes, "data.frame")){
    phenotypes <- AnnotatedDataFrame(phenotypes)
  }
  if (!is(phenotypes, "AnnotatedDataFrame")){
    stop("phenotypes must be a data.frame or an AnnotatedDataFrame.")
  }
  Biobase::phenoData(set) <- phenotypes
  
  if (is(annotationDataFrame, "data.frame")){
    annotationDataFrame <- AnnotatedDataFrame(annotationDataFrame)
  } 
  if (!is(annotationDataFrame, "AnnotatedDataFrame")){
    stop("annotationDataFrame must be a data.frame or an AnnotatedDataFrame.")
  }
  featureData(set) <- annotationDataFrame
  
  annotation(set) <- annoString
  return(set)
}