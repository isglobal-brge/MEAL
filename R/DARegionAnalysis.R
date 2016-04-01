#' Analyse methylation or expression in a specific range
#' 
#' Methylation analysis in a genomic range, taking into account snps.
#' 
#' @details Set is filtered to the range specified. If SNPs are present in the set, 
#' those are also filtered and then, correlation between SNPs and cpgs is tested. 
#' SNPs that are correlated to at least one cpg are added to covariables. After that,
#' \code{DAPipeline} is run. RDA test of the region is performed, returning the
#' R2 between the variables and the beta matrix and a p-value of this R2.   
#' 
#' @export DARegionAnalysis
#' @param set \code{MethylationSet}, \code{ExpressionSet} or \code{MultiDataSet}.
#' @param omicset In a \code{MultiDataSet} allows to choose between methylation and
#' expression (valid values are: "methylation" or "expression").
#' @param range \code{GenomicRanges} with the desired range.
#' @param variable_names Character vector with the names of the variables that
#' will be returned as result.
#' @param variable_types Character vector with the types of the variables. By
#' default, variables type won't be changed. 
#' @param covariable_names Character vector with the names of the variables that
#' will be used to adjust the model. 
#' @param covariable_types Character vector with the types of the covariables. By
#' default, variables type won't be changed. 
#' @param equation String containing the formula to be used to create 
#' the model.
#' @param num_var Numeric with the number of variables in the matrix for which the
#' analysis will be performed. Compulsory if equation is not null. 
#' @param labels Character vector with the labels of the variables.
#' @param sva Logical indicating if Surrogate Variable Analysis should be applied.
#' @param use_snps Logical indicating if SNPs should be used in the analysis. 
#' @param snps_cutoff Numerical with the threshold to consider a SNP-cpg correlation
#' p-value significant. 
#' @param region_methods Character vector with the methods used in \code{DARegion}. If
#' "none", region analysis is not performed.
#' @param shrinkVar Logical indicating if shrinkage of variance should be applied 
#' in probe analysis.
#' @param probe_method Character with the type of linear regression applied in 
#' probe analysis ("ls" or "robust")
#' @param max_iterations Numeric with the maximum of iterations in the robust
#' regression.
#' @param num_cores Numeric with the number of cores to be used.
#' @param verbose Logical value. If TRUE, it writes out some messages indicating progress. 
#' If FALSE nothing should be printed.
#' @param nperm Numeric with the number of permutations used to compute RDA p-values.
#' @param ... Further arguments passsed to \code{DAPipeline} function.
#' @return \code{AnalysisRegionResult} object 
#' @seealso \code{\link{preparePhenotype}}, \code{\link{DAPipeline}}
#' @examples
#' if (require(minfiData)){
#'  set <- prepareMethylationSet(getBeta(MsetEx)[1:1000, ], pheno = pData(MsetEx))
#'  range <- GenomicRanges::GRanges(seqnames=Rle("chrX"), 
#'  ranges = IRanges(30000, end = 123000000))
#'  res <- DARegionAnalysis(set, range = range, variable_names = "Sample_Group",
#'  probe_method = "ls") 
#'  res
#' }
DARegionAnalysis <- function(set, range, omicset = "methylation", variable_names, 
                             variable_types = rep(NA, length(variable_names)), 
                             covariable_names = NULL, 
                             covariable_types = rep(NA, length(covariable_names)), 
                             equation = NULL,   num_var = NULL,  labels = NULL, 
                             sva = FALSE, use_snps = TRUE, snps_cutoff = 0.01, 
                             region_methods = c("blockFinder", "bumphunter", "DMRcate"),
                             shrinkVar = FALSE, probe_method = "robust", max_iterations = 100, 
                             num_cores = 1, verbose = FALSE, nperm = 1e3,  ...){
  snps <- FALSE
  
  if (!(is(set, "MethylationSet") | is(set, "ExpressionSet") | is(set, "MultiDataSet"))){
    stop("set must be a MethylationSet, an ExpressionSet or a MultiDataSet.")
  }
  
  if (is(set, "MultiDataSet")){
    if (!"snps" %in% names(set)){
      stop("To use a MultiDataSet, snps must be present in the object.")
    }
    snpSet <- set[["snps"]]
    set <- set[[omicset]]
    snps <- TRUE
  }
  
  msg <- validObject(set, test = TRUE)
  
  options("mc.cores" = num_cores)
  
  if (is(msg, "character")){
    message(paste(msg, collapse = "\n"))
    stop("checkProbes and checkSamples might solve validity issues.")
  }
  if (ncol(set) == 0 | nrow(set) == 0){
    stop("The set has no beta values")
  }
  
  if (ncol(pData(set)) == 0 | nrow(pData(set)) == 0){
    stop("Set has no phenotypes")
  }
  
  if (!length(variable_names) | !sum(variable_names %in% colnames(pData(set)))){
    stop("variable_names is empty or is not a valid column of the phenoData of the set.")
  }
  
  if (!is(range, "GenomicRanges")){
    stop("range must be a GenomicRanges object.")
  }
  
  if (length(range) == 0){
    stop("range is empty.")
  }
  
  if (length(range) != 1){
    stop("range must be a GenomicRanges with only one range.")
  }
  
  
  if (verbose){
    message("Filtering the set.")
  }
  
  additives <- FALSE
  allvar <- c(variable_types, covariable_types)
  if (length(allvar) != sum(is.na(allvar))){
    addmask <- allvar == "additive"
    if (sum(addmask, na.rm = TRUE) > 0){
      allvar[addmask] <- "categorical"
      originalpheno <- preparePhenotype(phenotypes = pData(set), 
                                        variable_names = c(variable_names, covariable_names),
                                        variable_types = allvar)
      additives <- TRUE
    }
  }
  pData(set) <- preparePhenotype(phenotypes = pData(set), 
                                 variable_names = c(variable_names, covariable_names),
                                 variable_types = c(variable_types, covariable_types))
  
  ### Create model first to adjust for SVA using the whole dataset
  model <- createModel(data = pData(set), equation = equation, names = labels)
  if (sva){
    if (verbose){
      message("Computing SVA")
    }
    if (is(set, "ExpressionSet")){
      vals <- exprs(set)
    }else{
      vals <- betas(set)
    }
    if (nrow(vals) > 1){
      if (length(covariable_names) == 0){
        model0 <- model.matrix(~ 1, pData(set))
      }else{
        model0 <- model.matrix(~ ., pData(set)[ , covariable_names, drop = FALSE])
      }
      n.sv <- sva::num.sv(vals, model)
      if (n.sv > 0){
        svobj <- sva::sva(vals, model, model0, n.sv = n.sv)
        pData(set) <- cbind(pData(set), svobj$sv)
      }
    }
  }
  
  filt.set <- filterSet(set, range)
  if (!nrow(filt.set)){
    stop("There are no features in the range specified.")
  }  
  if (snps){
    if (verbose){
      message("Filtering the SNPs.")
    }
    snpSet <- filterSet(snpSet, range)
    
    if (nrow(snpSet) != 0){
      if (verbose){
        message("Obtaining SNPs P-values")
      }
      
      snpspvals <- calculateRelevantSNPs(filt.set, snps = snpSet, num_cores = num_cores)
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
        regionLM <- lapply(rownames(filt.set), function(x) 
          explainedVariance(data = data.frame(as.vector(getMs(filt.set[x, ])), 
                                              pData(filt.set)[ , variable_names],
                                              t(snpsgeno[snpspvals[ , x ] < snps_cutoff, , drop = FALSE]),
                                              pData(filt.set)[ , covariable_names]), 
                            num_mainvar = length(variable_names),
                            num_covariates = length(covariable_names),
                            variable_label = variable_label))
        names(regionLM) <- rownames(filt.set)
      }  
      #Prepare snps to be included in the phenotypes matrix
      snpsgeno <- snpsgeno[relevantsnps, , drop = FALSE]
      if (nrow(snpsgeno) != 0){
        snpsgeno <- snpsgeno[!duplicated(snpsgeno, MARGIN = 1), , drop = FALSE]
        snpsgeno <- normalSNP(snpsgeno)
        pc <- prcomp(snpsgeno)
        snpsgeno <- pc$rotation
        if (ncol(snpsgeno) > 6){
          snpsgeno <- snpsgeno[ , 1:6]
          snpsVar <- (cumsum(pc$sdev^2)/sum(pc$sdev^2))[6]
        } else {
          snpsVar <- 1
        }
        colnames(snpsgeno) <- paste0("snp", colnames(snpsgeno))
        pheno <- merge(pData(filt.set), snpsgeno, by=0)
        rownames(pheno) <- pheno$Row.names
        pheno <- pheno[colnames(betas(filt.set)), ]
        pheno <- pheno[ , -1]
        pData(filt.set) <- pheno
        covariable_names <- c(covariable_names, colnames(snpsgeno))
      } else {
        snpsVar <- as.numeric(NA)
      }
    }
  }
  results <- DAPipeline(filt.set, variable_names = variable_names, 
                        covariable_names = covariable_names, 
                        verbose = verbose, num_cores = num_cores, 
                        num_feat = nrow(filt.set), region_methods = region_methods, 
                        equation = equation, num_var = num_var, ...)
  if (additives){
    if (snps){
      originalpheno <- merge(originalpheno, snpsgeno, by = 0)
      originalpheno <- originalpheno[ , -1]
      rownames(originalpheno) <- originalpheno$Row.names
    }
    pData(set) <- originalpheno
  }
  if (snps){  
    output <- analysisRegionResults(analysisResults = results, set = set, range = range,
                                    snpspvals = snpspvals, snpsVar = snpsVar, 
                                    regionlm = regionLM, relevantsnps = relevantsnps, 
                                    equation = equation, nperm = nperm)
  }else{
    output <- analysisRegionResults(analysisResults = results, set = set, range = range,
                                    equation = equation, nperm = nperm)                           
  }
  return(output)
}