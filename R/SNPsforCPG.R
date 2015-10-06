# Calculate correlation cpg-SNP
# 
# Checks if SNPs are correlated with beta values of a cpg.
# 
# @param cpgname Character with the name of the cpg
# @param methy.data Matrix with the beta values.
# @param snp.data \code{SnpMatrix} with the genotypes
# @return Vector with the p-values of the correlation of the SNPs with the cpg. 
SNPsforCPG <- function (cpgname, methy.data, snp.data)
{
  cpg <- formula(paste(cpgname, "~ 1"))
  ans <- snpStats::snp.rhs.tests(cpg, data = methy.data, snp.data = snp.data, family = "Gaussian")
  pvals <- snpStats::p.value(ans)
  names(pvals) <- names(ans)
  pvals
}