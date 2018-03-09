#' Calculate RDA for a set
#' 
#' Perform RDA calculation for a \code{AnalysisRegionResults}. Feature values will
#' be considered the matrix X and phenotypes the matrix Y. Adjusting for covariates 
#' is done using a model matrix passed in covarsmodel.
#'  
#' @export runRDA
#' @param set \code{MethylationSet}, \code{ExpressionSet} or \code{matrix}
#' @param model Model matrix or formula to get model matrix from \code{set}. 
#' @param num_vars Numeric with the number of variables in the matrix for which the
#' analysis will be performed. Compulsory if equation is not null. 
#' @param range \code{GenomicRanges} with the region used for RDA
#' @param betas If \code{set} is a \code{GenomicRatioSet}, should beta values be
#' used? (Default: TRUE)
#' @param resultSet Should results be encapsulated in a \code{resultSet}? (Default: TRUE)
#' @param num_permutations Numeric with the number of permutations run to compute 
#' the p-value. (Default: 1e4)
#' @return Object of class \code{rda} or \code{resultSet} 
#' @seealso \code{\link[vegan]{rda}}
#' @examples
#' if (require(minfiData)){
#' set <- ratioConvert(mapToGenome(MsetEx[1:10,]))
#' model <- model.matrix(~set$age)
#' rda <- runRDA(set, model)
#' rda
#' }
runRDA <- function(set, model, num_vars = ncol(model), range, betas = FALSE, 
                   resultSet = TRUE, num_permutations = 1e4){
  
  ## Create model matrix from formula
  if (is(model, "formula")){
    model <- createModel(set, model)
    set <- set[, rownames(model)]
  }
  
  ## Get original matrix to compute region significance
  if (is(set, "ExpressionSet")){
    orimat <- Biobase::exprs(set)
  } else if (is(set, "GenomicRatioSet")){
    orimat <- minfi::getBeta(set)
    if (!betas) {
      orimat[orimat == 0] <- 1e-3
      orimat[orimat == 1] <- 1 - 1e-3
      orimat <- minfi::logit2(orimat)
    }
  } else if (is(set, "SummarizedExperiment")){
    orimat <- SummarizedExperiment::assay(set)
  } else if (is.matrix(set)){
    orimat <- set
  } else {
    stop("set must be a matrix, ExpressionSet, MethylationSet, GenomicRatioSet or SummarizedExperiment.")
  }
  
  
  if(!missing(range)){
    stopifnot(is(range, "GenomicRanges"))
    if (is(set, "eSet")){
      grset <- GenomicRanges::makeGRangesFromDataFrame(fData(set), 
                            start.field = c("position", "start"), 
                            end.field = c("position", "end"))
      names(grset) <- rownames(set)
      grset <- IRanges::subsetByOverlaps(grset, range)
      set <- set[names(grset), ]
    } else if (is(set, "RangedSummarizedExperiment")){
      set <- IRanges::subsetByOverlaps(set, range)
    } else if (is(set, "SummarizedExperiment")){
      grset <- GenomicRanges::makeGRangesFromDataFrame(rowData(set), 
                                                       start.field = c("position", "start"), 
                                                       end.field = c("position", "end"))
      names(grset) <- rownames(set)
      grset <- IRanges::subsetByOverlaps(grset, range)
      set <- set[names(grset), ]
    } else{
      stop("range parameter only works with eSet or SummarizedExperiment derived objects.")
    }
    
  } 
  
  ## Get region matrix
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
  
  varsmodel <- model[, seq_len(num_vars), drop = FALSE]
  if(num_vars == ncol(model)) {
    covarsmodel <- NULL
  }else {
    covarsmodel <- model[, (num_vars+1):ncol(model), drop = FALSE]
  }
  res <- vegan::rda(t(mat), varsmodel, covarsmodel)
  res$pval <- vegan::anova.cca(res, permutations = permute::how(nperm = num_permutations))[["Pr(>F)"]][1]
  res$rdaR2 <- vegan::RsquareAdj(res)$r.squared
  
  regvals <- computeRDAR2(fullMat = orimat, varsmodel = varsmodel, 
                          covarsmodel = covarsmodel, featNum = nrow(mat), 
                          R2 = res$rdaR2, num_permutations = num_permutations)
  
  res$globalR2 <- regvals["globalR2"] 
  res$globalpval<- regvals["globalpval"]

  if (resultSet)
  {
    res <- create_resultset("runRDA", lResults = 
                              list(RDA = list(result = res, error = NA)),  
                            fData = list(), lOptions = list())
  }
  res
  return(res)
}