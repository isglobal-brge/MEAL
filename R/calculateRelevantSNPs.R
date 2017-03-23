#' Calculate the SNPs correlated to cpgs
#' 
#' This function estimates the correlation between the snps and the cpgs. For each
#' pair cpg-SNP the p-value is returned.
#' 
#' @export calculateRelevantSNPs
#' @param set \code{MethylationSet}
#' @param snps \code{SnpSet}
#' @param num_cores Numeric with the number of cores to be used. 
#' @return Data.frame with the pvalues for pairs SNPs-cpgs. SNPs are in the rows 
#' and cpgs in the columns. 
#' @examples
#' \dontrun{
#' ## betamatrix: matrix of beta values
#' ## phenodf: data.frame with the phenotypes
#' ## snpsobject: SnpSet
#' set <- prepareMethylationSet(matrix = betamatrix, phenotypes = phenodf)
#' relevantSNPs <- calculateRelevantSNPs(set, snpsobject)
#' }

calculateRelevantSNPs <- function(set, snps, num_cores = 1){
  if (!is(set, "MethylationSet")){
    stop("set must be a MethylationSet")
  }
  if (nrow(set) == 0 || ncol(set) == 0){
    stop("set must contain beta values.")
  }
  
  if (nrow(snps) == 0 || ncol(snps) == 0){
    stop("SnpSet is empty.")
  }
  
  Ms <- data.frame(t(betas(set)))
  snps <- new("SnpMatrix", t(snpCall(snps)))
  
  if (nrow(snps) != nrow(Ms)){
    stop("Number of samples must be the same in beta matrix and in SnpMatrix")
  }
  
  snpsvscpgs <- parallel::mclapply(colnames(Ms), 
                                   function(x) SNPsforCPG(cpgname = x,
                                                          methy.data = Ms, snp.data = snps))
  snpsvscpgs <- data.frame(snpsvscpgs)
  colnames(snpsvscpgs) <- colnames(Ms)
  snpsvscpgs <- snpsvscpgs[apply(!is.na(snpsvscpgs), 1, all),]
  snpsvscpgs
}