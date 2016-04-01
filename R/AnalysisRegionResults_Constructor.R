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
#' @param nperm Numeric with the number of permutations used to compute RDAs p-values.
#' @return An \code{AnalysisRegionResults}
#' @examples
#' showClass("AnalysisRegionResults")
analysisRegionResults <- function(analysisResults, set, range, snpspvals = data.frame(),
                                     regionlm = list(), relevantsnps = character(),
                                     snpsVar = as.numeric(NA), equation = NULL, nperm = 1e5){
  phenoData(analysisResults) <- phenoData(set)
  rda <- RDAset(analysisResults, equation)
  pval <- vegan::anova.cca(rda, permutations = permute::how(nperm = nperm))[["Pr(>F)"]][1]
  r2 <- vegan::RsquareAdj(rda)$r.squared
  res <- computeRDAR2(set, analysisResults, r2, nperm = nperm)
  globpval <- res[1]
  globr <- res[2]
  class(rda) <- c("list", "cca", "rda")
  results <- new(Class = "AnalysisRegionResults", analysisResults, range = range, snpsPvals = snpspvals, 
                 regionLM = regionlm, snps = relevantsnps, snpsVar = snpsVar, 
                 rda = rda, regionR2 = r2, RDAPval = pval, originalclass = class(set), globalR2 = globr, globalPval = globpval)
  return(results)
}