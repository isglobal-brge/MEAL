#' Get all probes related to a gene
#' 
#' Given a \code{ResultSet} and a gene name returns the results of the 
#' analysis of all the probes of the gene. 
#' 
#' @name getGeneVals
#' @aliases getGeneVals
#' @export
#' @param object \code{ResultSet}
#' @param gene Character with the name of the gene
#' @param rid Name of the results: "DiffMean" for mean differences, "DiffVar" for
#' variance differences. (Default: DiffMean)
#' @param genecol Character with the column of \code{object} fData with the gene
#' information
#' @param ... Further arguments passed to \code{getProbeResults}
#' @return data.frame with the results of the analysis of the probes belonging to 
#' the gene 
#' @examples
#' \dontrun{
#' if (require(minfiData)){
#'  set <- ratioConvert(mapToGenome(MsetEx[1:10,]))
#' methyOneVar <- runPipeline(set, variable_names = "sex")
#' getGeneVals(methyOneVar, "TSPY4")
#' }
#' }
getGeneVals <- function(object, gene, rid = 1, genecol = "genes", ...){
  
  stopifnot(is(object, "ResultSet"), is(gene, "character"))
  
  res <- getProbeResults(object, rid, ...)
  
  stopifnot(genecol %in% colnames(res))
  
  genelist <- strsplit(res$genes, ";")
  mask <- vapply(genelist, function(y) gene %in% y, logical(length(gene)))
  if (length(gene) > 1){
    mask <- as.logical(colSums(mask))
  }
  return(res[mask, , drop = FALSE])
  if (nrow(res) == 0){
    warning("Gene name is not present in the results. An empty data.frame will be returned.")
  }
  res
}