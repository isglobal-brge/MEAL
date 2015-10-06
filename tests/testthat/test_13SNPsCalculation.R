context("Snps p-vals")
 
library(minfiData)

set <- prepareMethylationSet(getBeta(MsetEx)[1:10,], pheno = pData(MsetEx))

##Create multiset with snps
geno <- matrix(rep(c(2, 1, 2, 1, 1, 1), 2), ncol = 6)
colnames(geno) <- c("5723646052_R02C02", "5723646052_R04C01", "5723646052_R05C02",
                    "5723646053_R04C02", "5723646053_R05C02", "5723646053_R06C02")
rownames(geno) <- c("rs3115860", "SNP1-1628854")
map <- data.frame(chromosome = c("chrY", "chr2"), position = c(4241114, 1234321),
                  stringsAsFactors = FALSE)
rownames(map) <- rownames(geno)
snps <- new("SnpSet", call = geno)
fData(snps) <- map

range <- GenomicRanges::GRanges(seqnames=Rle("chrY"), 
                                ranges = IRanges(3000000, end=12300000))
snps <- filterSet(snps, range)

test_that("p-vals calculation works", {
  pvals <- calculateRelevantSNPs(set, snps)
  expect_equal(nrow(pvals), 1)
  expect_equal(ncol(pvals), as.integer(nrow(set)))
})

test_that("Empty variables", {
  emptyset <- new(Class = "MethylationSet")
  emptySnps <- new(Class = "SnpSet")
  expect_error(calculateRelevantSNPs(emptyset, snps), "set must contain beta values.")
  expect_error(calculateRelevantSNPs(set, emptySnps), "SnpSet is empty.")
  library(minfiData)
  expect_error(calculateRelevantSNPs(MsetEx, snps), "set must be a MethylationSet")
  
})