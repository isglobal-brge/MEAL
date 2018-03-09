#' Run differential variance analysis
#' 
#' Run differential variance analysis. This analysis can only be run with categorical 
#' variables. This function relies on \code{varFit} from missMethyl package.
#' 
#' @export runDiffVarAnalysis
#' 
#' @param set Matrix, \code{GenomicRatioSet}, \code{SummarizedExperiment} or 
#' \code{ExpressionSet}. 
#' @param model Model matrix or formula to get model matrix from \code{set}. 
#' @param coefficient Numeric with the coefficients used to make the groups. If NULL, all
#' possible groups will be computed. 
#' @param resultSet Should results be encapsulated in a \code{resultSet}? (Default: TRUE)
#' @param betas If \code{set} is a \code{GenomicRatioSet}, should beta values be
#' used? (Default: TRUE)
#' @param warnings Should warnings be displayed? (Default:TRUE) 
#' @param ... Further arguments passed to \code{varFit}.
#' @return \code{MArrayLM} or \code{resultSet} with the result of the differential
#' variance analysis.
#' @examples
#' if (require(minfiData)){
#'  mvalues <- getM(MsetEx)[1:100, ]
#'  model <- model.matrix(~ Sample_Group, data = pData(MsetEx)) 
#'  res <- runDiffVarAnalysis(mvalues, model)
#'  res
#' }
runDiffVarAnalysis <- function(set, model, coefficient = NULL,  
                            resultSet = TRUE, betas = TRUE, warnings = TRUE, ...) {
  
  ## Create model matrix from formula
  if (is(model, "formula")){
    model <- createModel(set, model)
    set <- set[, rownames(model)]
  }
  
  
  ## Get matrix
  if (is(set, "ExpressionSet")){
    mat <- Biobase::exprs(set)
  } else if (is(set, "GenomicRatioSet")){
    mat <- minfi::getBeta(set)
    if (!betas) {
      mat[mat == 0] <- 1e-3
      mat[mat == 1] <- 1 - 1e-3
      mat <- minfi::logit2(mat)
    }
  }else if (is(set, "SummarizedExperiment")){
    mat <- SummarizedExperiment::assay(set)
  } else if (is.matrix(set)){
    mat <- set
  } else {
    stop("set must be a matrix, ExpressionSet, MethylationSet, GenomicRatioSet or SummarizedExperiment.")
  }

  fFun <- getFeatureDataFun(set)
  
  fit <- missMethyl::varFit.default(mat, model, coef = coefficient, ...)
  
  if (resultSet){
    fit <- MultiDataSet::create_resultset(fOrigin = "DiffVarAnalysis", fData = list(main = fFun(set)), 
                                          lResults = list(DiffVar = list(result = fit)))
  } 
  fit
}