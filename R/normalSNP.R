#' Normalize SNPs values
#' 
#' SNPs values, introduced as numerical, are normalized to be used in lineal models.
#' 
#' @export normalSNP
#' @param snps Numerical vector or matrix representing the SNPs in the form: 0 homozygote
#' recessive, 1 heterozygote, 2 homozygote dominant.
#' @return Numerical vector or matrix with the snps normalized.
#' @examples
#' snps <- c(1, 0, 0, 1, 0, 0, 2, 1, 2)
#' normSNPs <- normalSNP(snps)
#' normSNPs

normalSNP <- function(snps)
{
  if (is(snps, "data.frame")){
    snps <- data.matrix(snps)
  }
  if (is(snps, "matrix")){
    means <- matrix(rep(rowMeans(snps), rep(ncol(snps), nrow(snps))), ncol = ncol(snps), byrow = TRUE)
    f <- sqrt(means/2*(1-means/2))
    ans <- (snps - means) / f
  }else{
  mu <- mean(snps)
  f <- mu/2    
  ans <- (snps - mu)/sqrt(f*(1-f)) 
  }
  return(ans)
}
