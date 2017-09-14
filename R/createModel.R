# Create model.matrix from formula and omic set
# 
# Given an \code{eSet} or \code{SummarizedExperiment}, create a model matrix using
# formula. If there are missing values in the phenotype, the original set is subsetted
# to complete cases.
# 
#  
# @param set \code{eSet} or \code{SummarizedExperiment}
# @param model Formula with the model
# @param warnings Should warnings be displayed? (Default:TRUE)
# @return Matrix with the model matrix. 
createModel <- function(set, model, warnings = TRUE){
  ## Change getters depending on set class
  if (is(set, "eSet")){
    pFun <- Biobase::pData
  }
  else if (is(set, "SummarizedExperiment")){
    pFun <- SummarizedExperiment::colData
  } else{
    stop("set must be an eSet or a SummarizedExperiment derived object")
  }
  
  ## Create model matrix
  model <- model.matrix(model, pFun(set))
  
  ## Check missing data
  if (nrow(model) != ncol(set)){
    if (warnings){
      warning("There are some missing values in the samples data. Only complete cases will be used.")
    }
    set <- set[, rownames(model)]      
  }
  
  ## Check other paramaters
  if (ncol(set) == 0 | nrow(set) == 0){
    stop("The set is empty.")
  }
  if (ncol(model) == 0 | nrow(model) == 0){
    stop("The model matrix is empty.")
  }
  if (ncol(set) != nrow(model)){
    stop("The number of samples is different in the set and in the model")
  }
  
  return(model)
}