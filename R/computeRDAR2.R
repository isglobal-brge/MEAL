#' Compute signification of RDA test
#' 
#' Compare R2 obtained in our region of interest with the global R^2 and the R^2 of regions with 
#' the same number of probes.
#' @export
#' @param fullMat Matrix with the whole genome expression or methylation values
#' @param varsmodel Matrix with the model
#' @param covarsmodel Matrix with the covariables model
#' @param featNum Numeric with the number of features of the RDA model
#' @param R2 Numeric with the R2 of the RDA model
#' @param nperm Numeric with the number of permutations.
#' @return Numeric vector with the probability of finding a region with the same number of probes
#' with a bigger R2 and the global R2.
computeRDAR2 <- function(fullMat, varsmodel, covarsmodel = NULL, featNum, R2, nperm = 1e6-1){

  allvar <- apply(fullMat, 1, stats::var)

   ## If there are covariables, substract their effect prior to adjusting 
  if (!is.null(covarsmodel)){
    covfit <- limma::fitted.MArrayLM(limma::lmFit(fullMat, covarsmodel))
    fullMat <- fullMat - covfit
  }
  fit <- limma::fitted.MArrayLM(limma::lmFit(fullMat, varsmodel))
  varfit <- apply(fit, 1, stats::var)
  
  inds <- sample(1:(nrow(fullMat) - featNum), nperm, replace = TRUE)
  
  ## Compute the R2 explained by groups with the same number of features than the region of interest
  randomR2 <- vapply(inds, function(x) { 
      R2 <- sum(varfit[x:(x + featNum - 1)])/sum(allvar[x:(x + featNum - 1)])}, numeric(1))
  randomR2 <- c(randomR2, R2)
  pval <- sum(randomR2 >= R2)/(nperm + 1)

  globalR2 <- sum(varfit)/sum(allvar)
  c(globalpval = pval, globalR2 = globalR2)
}
