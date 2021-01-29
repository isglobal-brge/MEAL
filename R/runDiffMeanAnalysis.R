#' Run differential mean analysis
#' 
#' Run differential mean analysis using t-moderated statistics. This function relies
#' on \code{lmFit} from limma package.
#' 
#' @export runDiffMeanAnalysis
#' 
#' @param set Matrix, \code{GenomicRatioSet}, \code{SummarizedExperiment} or 
#' \code{ExpressionSet}. 
#' @param model Model matrix or formula to get model matrix from \code{set}. 
#' @param weights weights used in the lmFit model.
#' @param method String indicating the method used in the regression: "ls" or 
#' "robust". (Default: "ls")
#' @param max_iterations Numeric indicating the maximum number of iterations
#' done in the robust method. 
#' @param betas If \code{set} is a \code{GenomicRatioSet}, should beta values be
#' used? (Default: TRUE)
#' @param resultSet Should results be encapsulated in a \code{resultSet}? (Default: TRUE)
#' @param warnings Should warnings be displayed? (Default:TRUE) 
#' @return \code{MArrayLM} or \code{resultSet} with the result of the differential
#' mean analysis.
#' @examples
#' if (require(minfiData)){
#'  mvalues <- getM(MsetEx)[1:100, ]
#'  model <- model.matrix(~ Sample_Group, data = pData(MsetEx)) 
#'  res <- runDiffMeanAnalysis(mvalues, model, method = "ls")
#'  res
#' }
runDiffMeanAnalysis <- function(set, model, weights = weights,  method = "ls", max_iterations = 100, 
                    betas = TRUE, resultSet = TRUE, warnings = TRUE, ...) {
  
  ## Create model matrix from formula
  if (is(model, "formula")){
    model <- createModel(set, model)
    set <- set[, rownames(model)]
  }
  
 method <- match.arg(method, choices = c("ls", "robust"))

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
  } else if (is(set, "SummarizedExperiment")){
    mat <- SummarizedExperiment::assay(set)
  } else if (is.matrix(set)){
    mat <- set
  } else {
    stop("set must be a matrix, ExpressionSet, MethylationSet, GenomicRatioSet or SummarizedExperiment.")
  }
  
  fFun <- getFeatureDataFun(set)

  if (length(max_iterations) == 0 || !is(max_iterations, "numeric") || max_iterations < 1){
    warning("max_iterations must be a numeric greater than 1. The default value (100) will be used.")
    max_iterations <- 100
  }

  fit <- limma::lmFit(mat, model, method = method, maxit = max_iterations, weights = weights)

  if (resultSet){
    fit <- MultiDataSet::create_resultset(fOrigin = "DAProbe", fData = list(main = fFun(set)), 
                                   lResults = list(DiffMean = list(result = fit)))
  } 
  fit
}