#' Compute signification of RDA test
#' 
#' Compare R2 obtained in our region of interest with the global R^2 and the R^2 of regions with 
#' the same number of probes.
#' @export
#' @param set \code{MethylationSet} or \code{ExpressionSet} 
#' @param res \code{AnalysisResults}
#' @param R2 Numeric with the R2 obtained in the RDA
#' @param nperm Numeric with the number of permutations.
#' @return Numeric vector with the probability of finding a region with the same number of probes
#' with a bigger R2 and the global R2.
computeRDAR2 <- function(set, res, R2, nperm = 1e6-1){

  if (is(set, "ExpressionSet")){
    matMs <- exprs(set)
  } else{
    matMs <- getMs(set)
  }
  allvar <- apply(matMs, 1, stats::var)
  ## Get the model used in the regression
  varmatrix <- model(res)
  covs <- covariableNames(res)
  
  ## If there are covariables, substract their effect prior to adjusting 
  if (length(covs)){
    covsind <- lapply(covs, function(x) grep(x, colnames(varmatrix)))
    covsind <- unique(unlist(covsind))
    covmatrix <- varmatrix[ , c(1, covsind), drop = FALSE]
    varmatrix <- varmatrix[ , -covsind, drop = FALSE]
    covfit <- limma::fitted.MArrayLM(limma::lmFit(matMs, covmatrix))
    matMs <- matMs - covfit
  }
  fit <- limma::fitted.MArrayLM(limma::lmFit(matMs, varmatrix))
  varfit <- apply(fit, 1, stats::var)
  num_feats <- length(feats(res))
  
  inds <- sample(1:(nrow(set) - num_feats), nperm, replace = TRUE)
  
  ## Compute the R2 explained by groups with the same number of features than the region of interest
  randomR2 <- vapply(inds, function(x) { 
      R2 <- sum(varfit[x:(x + num_feats - 1)])/sum(allvar[x:(x + num_feats - 1)])}, numeric(1))
  randomR2 <- c(randomR2, R2)
  pval <- sum(randomR2 >= R2)/(nperm + 1)

  globalR2 <- sum(varfit)/sum(allvar)
  c(pval, globalR2)
}
