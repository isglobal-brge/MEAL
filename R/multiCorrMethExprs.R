#' Computes the correlation between methylation and expression in a genomic range 
#' 
#' Estimates the correlation between methylation and expression in a range. First, 
#' the sets are filtered to only contain the features of the range. Then, a multivariate
#' approach (redundancy analysis) is applied. 
#' 
#' When there are known variables that affect methylation and/or expression, their effect 
#' can be substracted using a linear model and then the residuals are used. 
#' 
#' @export multiCorrMethExprs
#' 
#' @param multiset \code{MultiDataSet} containing a \code{methylation} and an 
#' \code{expression} slots.
#' @param vars_meth Character vector with the names of the variables that will be
#' used to obtain the methylation residuals. By default, none is used and residuals 
#' are not computed.
#' @param vars_meth_types Character vector with the types of the methylation variables. 
#' By default, variables type won't be changed. 
#' @param vars_exprs Character vector with the names of the variables that will
#' be used to obtain the expression residuals. By default, none is used and
#' residuals are not computed.
#' @param vars_exprs_types Character vector with the types of the expression variables. 
#' By default, variables type won't be changed.
#' @param range \code{GenomicRanges} with the desired range. 
#' @return An \code{rda} object
multiCorrMethExprs <- function(multiset, vars_meth = NULL, vars_exprs = NULL, 
                                 vars_meth_types = rep(NA, length(vars_meth)), 
                                 vars_exprs_types = rep(NA, length(vars_exprs)), 
                                 range = NULL){
  
  if (!is(multiset, "MultiDataSet")){
    stop("multiset must be a MultiDataSet")
  }
  
  if (!all(c("methylation", "expression") %in% names(multiset))){
    stop("multiset must contain methylation and expression data.")
  }
  
  mset <- multiset[["methylation"]]
  eset <- multiset[["expression"]]
  
  if (ncol(mset) == 0 | nrow(mset) == 0){
    stop("The mset has no beta values")
  }
  
  if (ncol(eset) == 0 | nrow(eset) == 0){
    stop("The eset has no expression values")
  }
  
  if (!all(c("chromosome", "start", "end") %in% fvarLabels(eset))){
    stop("eset must contain a featureData with columns chromosome, start and end")
  }
  
  mset <- filterSet(set = mset, range = range)
  if (!nrow(mset)){
    stop("There are no cpgs in the specified range")
  }  

  eset <- filterSet(set = eset, range = range)
  if (!nrow(eset)){
    stop("There are no expression probes in the specified range")
  }  
  exprsres <- setResidues(eset, vars_names = vars_exprs, vars_types = vars_exprs_types)
  exprsres <- t(data.matrix(exprsres))
  
  methres <- setResidues(mset, vars_names = vars_meth, vars_types = vars_meth_types)
  methres <- t(data.matrix(minfi::logit2(methres)))
  
  ## Given that methylation data will have more probes than samples in expression set, we reduce its dimensionality
  pca <- prcomp(methres)
  methpc <- pca$x[, which(pca$sdev > mean(pca$sdev)), drop = FALSE]
  
  rd <- vegan::rda(exprsres, methpc)
  
  ## The biplot is computed using the correlation with the original variables
  rd$CCA$biplot <- cor(methres, rd$CCA$u)
  return(rd)
}