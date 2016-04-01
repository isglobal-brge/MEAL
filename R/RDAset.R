#' Calculate RDA for a set
#' 
#' Perform RDA calculation for a \code{AnalysisRegionResults}. Feature values will
#' be considered the matrix X and phenotypes the matrix Y. Adjusting for covariates 
#' is done using covariable_names stored in the object.
#'  
#' @export RDAset
#' @param set \code{AnalysisResults}
#' @param equation Character with the equation used in the analysis
#' @return Object of class \code{rda}
#' @seealso \code{\link[vegan]{rda}}
#' @examples
#' if (require(minfiData)){
#' set <- prepareMethylationSet(getBeta(MsetEx)[1:50,], pheno = pData(MsetEx))
#' methyOneVar <- DAPipeline(set, variable_names = "sex", probe_method = "ls")
#' rda <- RDAset(methyOneVar)
#' rda
#' }
RDAset <- function(set, equation = NULL){
  if (!is(set, "AnalysisResults")){
    stop("set must be a AnalysisResults")
  }
  if (ncol(featvals(set)) == 0 || nrow(featvals(set)) == 0){
    stop("set has no beta values.")
  }
  if (ncol(pData(set)) == 0 || nrow(pData(set)) == 0){
    stop("set has no phenotypic information.")
  }
  if (length(variableNames(set)) == 0){
    stop("set must contain variable names.")
  }
  
  if (set@originalclass == "ExpressionSet"){
    vals <- t(featvals(set))
  } else{
    vals <- t(getMs(set))
  }
  
  
  varmatrix <- model(set)
  covs <- covariableNames(set)
  if (length(covs)){
    covsind <- lapply(covs, function(x) grep(x, colnames(varmatrix)))
    covsind <- unique(unlist(covsind))
    covmatrix <- varmatrix[ , covsind, drop = FALSE]
    varmatrix <- varmatrix[ , -covsind, drop = FALSE]
    rdaRes <- vegan::rda(vals, varmatrix, covmatrix)
  }else{
    rdaRes <- vegan::rda(vals, varmatrix)
  }
  
  return(rdaRes)
}