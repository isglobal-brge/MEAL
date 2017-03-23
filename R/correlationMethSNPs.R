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
#' @param vars_meth Character vector with the names of the variables that will be
#' used to obtain the methylation residuals. By default, none is used and residuals 
#' are not computed.
#' @param vars_meth_types Character vector with the types of the methylation variables. 
#' By default, variables type won't be changed. 
#' @param vars_exprs Character vector with the names of the variables that will
#' be used to obtain the expression residuals. By default, none is used and
#' residuals are not computed.
#' @param vars_exprs_types Character vector with the types of the expression variables. 
#' By default, variables type won't be changed.
#' @param meth_set_name Character vector with the name of the \code{MultiDataSet}'s slot containing methylation
#' data.
#' @param exprs_set_name Character vector with the name of the \code{MultiDataSet}'s slot containing expression
#' data.
#' @param sel_cpgs Character vector with the name of the CpGs used in the analysis. If empty, all the CpGs of the 
#' methylation set will be used. 
#' @param snps_cutoff Numerical with the threshold to consider a p-value from a SNP-cpg correlation
#' significant. 
#' @param flank Numeric with the number of pair bases used to define the cpg-expression 
#' probe pairs.
#' @param num_cores Numeric with the number of cores to be used.
#' @param verbose Logical value. If TRUE, it writes out some messages indicating progress. 
#' If FALSE nothing should be printed.
#' @return Data.frame with the results of the linear regression:
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

  multi <- multi[c(meth_set_name, snps_set_name)]
  
  if (!is.null(range)){
    if (!is(range, "GenomicRanges")){
      stop("Range should be empty or a GenomicRanges")
    }
    multi <- multi[, , range]
    
    if (any(nrows(multi) == 0)){
      stop("There are no SNPs or CpGs in the region")
    } 
  }
  
  snpSet <- multi[[snps_set_name]]
  methSet <- multi[[meth_set_name]]
  
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


### Testing

## test_08 ####
##Create multiset with snps
# geno <- matrix(rep(c(2, 1, 2, 1, 1, 1), 2), ncol = 6)
# colnames(geno) <- c("5723646052_R02C02", "5723646052_R04C01", "5723646052_R05C02",
#                     "5723646053_R04C02", "5723646053_R05C02", "5723646053_R06C02")
# rownames(geno) <- c("rs3115860", "SNP1-1628854")
# map <- data.frame(chromosome = c("chrY", "chr2"), position = c(4241114, 1234321),
#                   stringsAsFactors = FALSE)
# rownames(map) <- rownames(geno)
# snps <- new("SnpSet", call = geno)
# fData(snps) <- map
# 
# multiset <- new("MultiDataSet")
# # multiset <- add_methy(multiset, set)
# multiset <- add_snps(multiset, snps)
# 
# results <- DARegionAnalysis(set = multiset, variable_names = "sex", 
#                             variable_types = "categorical", range = range)
# expect_match(class(results), "AnalysisRegionResults")
# expect_error(DARegionAnalysis(set = multiset, variable_names = character(), 
#                               variable_types = character(), range = range), "variable_names is empty.")
# expect_error(DARegionAnalysis(set = multiset, variable_names = "sex", 
#                               variable_types = character(), range = range), "variable_types is empty.")
# expect_error(DARegionAnalysis(set = multiset, variable_names = character(), 
#                               covariable_names = "sex", covariable_types = "categorical",
#                               range = range),
#              "variable_names is empty or is not a valid column of the phenoData of the set.")
# expect_error(DARegionAnalysis(set = multiset, variable_names = "cot", covariable_names = "sex",
#                               covariable_types = "categorical", range = range),
#              "variable_names is empty or is not a valid column of the phenoData of the set.")

# ## test_11 ####
# ##Create multiset with snps
# geno <- matrix(rep(c(3, 1, 3, 1, 1, 1), 2), ncol = 6, byrow = T)
# colnames(geno) <- c("5723646052_R02C02", "5723646052_R04C01", "5723646052_R05C02",
#                     "5723646053_R04C02", "5723646053_R05C02", "5723646053_R06C02")
# rownames(geno) <- c("rs3115860", "SNP1-1628854")
# map <- AnnotatedDataFrame(data.frame(chromosome = c("chrY", "chr2"), position = c(4241114, 1234321),
#                                      stringsAsFactors = FALSE))
# rownames(map) <- rownames(geno)
# snps <- new("SnpSet", call = geno, featureData = map)
# 
# multiset <- new("MultiDataSet")
# multiset <- add_methy(multiset, set)
# multiset <- add_snps(multiset, snps)
# rangeSNPsCov <- DARegionAnalysis(multiset, variable_names = "sex", range = range, 
#                                  snps_cutoff = 0.05)
# 
# test_that("Plot Region R2", {  
#     expect_error(plotRegionR2(rangeSNPsCov, 12351413654), "feat index must be greater than 0 and smaller than the number of cpgs.")
#     expect_error(plotRegionR2(rangeSNPsCov, -12), "feat index must be greater than 0 and smaller than the number of cpgs.")
#     expect_error(plotRegionR2(rangeSNPsCov, "12"), "feat name is not present in the set.")
#     expect_error(plotRegionR2(rangeSNPsCov, 1:4), "feat must contain only one value")
#     expect_error(plotRegionR2(rangeSNPsCov, 1:5), "feat must contain only one value.")
#     plotRDA(rangeSNPsCov)
# })
