#' AnalysisRegionResults instances
#' 
#' @export
#' 
#' @name analysisRegionResults
#' @rdname AnalysisRegionResults-class
#' @aliases AnalysisRegionResults-methods
#' 
#' @param analysisResults \code{AnalysisResults}
#' @param set \code{MethylationSet} or \code{ExpressionSet}
#' @param range \code{GenomicRanges}
#' @param snpspvals Data.frame obtained from \code{calculateRelevantSNPs}
#' @param regionlm Data.frame obtained from \code{explainedVariance}
#' @param snpsVar Numeric with the variability of the SNP matrix explained by the components 
#' used to adjust the linear model. 
#' @param relevantsnps Character vector with the relevant snps names
#' @param equation Character containing the formula to be used to create the model.
#' @return An \code{AnalysisRegionResults}
#' @examples
#' showClass("AnalysisRegionResults")
analysisRegionResults <- function(analysisResults, set, range, snpspvals = data.frame(),
                                     regionlm = list(), relevantsnps = character(),
                                     snpsVar = as.numeric(NA), equation = NULL){
  
  if (is(set, "ExpressionSet")){
    features <- exprs(set)
  }else{
    features <- betas(set)
  }
  phenotypes <- phenoData(set)
  phenoData(analysisResults) <- phenoData(set)
  rda <- RDAset(analysisResults, equation)
  pval <- round(vegan::anova.cca(rda)[["Pr(>F)"]][1], 3)
  r2 <- round(vegan::RsquareAdj(rda)$r.squared, 3)
  class(rda) <- c("list", "cca", "rda")
  results <- new(Class = "AnalysisRegionResults", analysisResults, features = features, 
                 phenotypes = phenotypes, range = range, snpsPvals = snpspvals, 
                 regionLM = regionlm, snps = relevantsnps, snpsVar = snpsVar, 
                 rda = rda, regionR2 = r2, regionPval = pval, originalclass = class(set))
  return(results)
}