#' Plot best n cpgs
#' 
#' Wrapper of \code{plotCPG} that plots the top n features.
#' 
#' @export
#' @param set \code{AnalysisResults}, \code{AnalysisRegionResults}, \code{ExpressionSet} or 
#' \code{MethylationSet}
#' @param n Numeric with the number of features to be plotted.
#' @param variables Character vector with the names of the variables to be used 
#' in the splitting. 
#' @return Plots are created on the current graphics device.
#' @seealso \code{\link{plotFeature}}
#' @examples
#' if (require(minfiData)){
#'  set <- prepareMethylationSet(getBeta(MsetEx)[1:10, ], 
#'  pheno = data.frame(pData(MsetEx)))
#'  plotBestFeatures(set, 2, variables = "Sample_Group")
#'  }
plotBestFeatures <- function(set, n = 10, variables = variableNames(set)[1]){
  if (is(set, "AnalysisResults")){
    rows <- nrow(featvals(set))
  }else if (is(set, "ExpressionSet") || is(set, "MethylationSet")){
    rows <- nrow(set)
  }else{
    stop("set must be of class AnalysisResults, ExpressionSet or MethylationSet")
  }
  
  if (n > rows){
    warning("n is greater than the number of betas present in the set.")
    n <- rows
  }
  
  a <- lapply(1:n, function(x) plotFeature(set = set, feat = x, variables = variables))
  invisible(a)
}