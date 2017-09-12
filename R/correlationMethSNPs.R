#### Revisar !!!!!!!!!!!!!!



#' Computes the correlation between methylation and SNPs
#' 
#' Estimates the correlation between methylation and expression. When there are known 
#' variables that affect methylation and/or expression, their effect can be substracted
#' using a linear model and then the residuals are used. 
#' 
#' For each cpg, a range is defined by the position of the cpg plus the flank parameter
#' (upstream and downstream). Only those expression probes that are entirely in 
#' this range will be selected. For these reason, it is required that the \code{ExpressionSet} 
#' contains a featureData with the chromosome and the starting and ending positions 
#' of the probes.  
#' 
#' @export correlationMethSNPs
#' 
#' @param multiset \code{MultiDataSet} containing a \code{methylation} and an 
#' \code{expression} slots.
#' @param meth_set_name Character vector with the name of the \code{MultiDataSet}'s slot containing methylation
#' data.
#' @param snps_set_name Character vector with the name of the \code{MultiDataSet}'s slot containing SNPs
#' data.
#' @param variable_names Character vector with the names of the variables that will be
#' used to obtain the methylation residuals. By default, none is used and residuals 
#' are not computed.
#' @param covariable_names Character vector with the names of the variables that
#' will be used to adjust the model. 
#' @param range \code{GenomicRanges} with the range used in the analñysis
#' @param snps_cutoff Numerical with the threshold to consider a p-value from a SNP-cpg correlation
#' significant. 
#' @param verbose Logical value. If TRUE, it writes out some messages indicating progress. 
#' If FALSE nothing should be printed.
#' @return List with the results:
#' \itemize{
#'  \item cpg: Name of the cpg
#'  \item exprs: Name of the expression probe
#'  \item beta: coefficient of the methylation change
#'  \item se: standard error of the beta 
#'  \item P.Value: p-value of the beta coefficient
#'  \item adj.P.Val: q-value computed using B&H
#' }
correlationMethSNPs <- function(multiset, meth_set_name = NULL, snps_set_name = NULL,
                                 range, variable_names, covariable_names = NULL, snps_cutoff = 0.01, 
                                verbose = TRUE){
  
  ## Data Checking
  # Check object is a MultiDataSet
  if (!is(multiset, "MultiDataSet")){
    stop("multiset must be a MultiDataSet")
  }
  
  # Check meth_set_name and snps_set_name are characters
  if (!(is.character(meth_set_name) | is.null(meth_set_name))){
    stop("meth_set_name must be a character.")
  } 
  if (!(is.character(snps_set_name) | is.null(snps_set_name))){
    stop("snps_set_name must be a character.")
  } 
  
  # Add the dataset type to the names
  meth_set_name <- paste(c("methylation", meth_set_name), collapse = "+")
  snps_set_name <- paste(c("snps", snps_set_name), collapse = "+")
  
  
  # Check our object has the right sets
  if (!all(c(meth_set_name, snps_set_name) %in% names(multiset))){
    stop("multiset must contain meth_set_name and snps_set_name.")
  }
  
  if (!all(c(meth_set_name, snps_set_name) %in% rowRangesElements(multiset))){
    stop("multiset must contain meth_set_name and snps_set_name.")
  }

  multiset <- multiset[, c(meth_set_name, snps_set_name)]
  
  if (!is.null(range)){
    if (!is(range, "GenomicRanges")){
      stop("Range should be empty or a GenomicRanges")
    }
    multiset <- multiset[, , range]
    
    if (any(MultiDataSet:::nrows(multiset) == 0)){
      stop("There are no SNPs or CpGs in the region")
    } 
  }
  
  snpSet <- multiset[[snps_set_name]]
  methSet <- multiset[[meth_set_name]]
  
  if (!length(variable_names) | !sum(variable_names %in% colnames(pData(methSet)))){
    stop("variable_names is empty or is not a valid column of the phenoData of the methylation Set.")
  }
  
  if (verbose){
    message("Obtaining SNPs P-values")
  }
  
  snpspvals <- calculateRelevantSNPs(methSet, snps = snpSet, num_cores = num_cores)
  relevantsnps <- rownames(snpspvals)[sapply(rownames(snpspvals),
                                             function(x) any(snpspvals[x, ] < snps_cutoff))]
  snpsgeno <- snpCall(snpSet)
  snpsgeno <- snpsgeno[rownames(snpspvals), , drop = FALSE]
  variable_label <- paste0(variable_names, collapse = "")
  
  if (verbose){
    message("Calculating variables R2.")
  }
  
  if (sum(snpspvals < snps_cutoff) == 0){
    regionLM <- list()
  }else{
    regionLM <- lapply(rownames(methSet), function(x)
      explainedVariance(data = data.frame(as.vector(MultiDataSet::betas(methSet[x, ])),
                                          pData(methSet)[ , variable_names],
                                          t(snpsgeno[snpspvals[ , x ] < snps_cutoff, , drop = FALSE]),
                                          pData(methSet)[ , covariable_names]),
                        num_mainvar = length(variable_names),
                        num_covariates = length(covariable_names),
                        variable_label = variable_label))
    names(regionLM) <- rownames(methSet)
  }
  return(list(LM = regionLM, pvals = snpspvals))
}


