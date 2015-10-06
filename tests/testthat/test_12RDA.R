context("RDA calculation")

set.seed(0)
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

multiset <- new("MultiDataSet")
multiset <- add.methy(multiset, set)
multiset <- add.snps(multiset, snps)

range <- GenomicRanges::GRanges(seqnames=Rle("chrY"), 
                                ranges = IRanges(3000000, end=12300000))
rangeSNPs <- DARegionAnalysis(multiset, variable_names = "sex", range = range)



test_that("RDA calculation works", {
  rda <- RDAset(rangeSNPs)
  expect_equal(class(rda), c("rda", "cca"))
})

test_that("wrong variables", {
  expect_error(RDAset(multiset), "set must be a AnalysisResults")
  
  emptyset <- new(Class = "AnalysisRegionResults")
  expect_error(RDAset(emptyset), "set has no beta values.")
  
  pData(rangeSNPs) <- data.frame()
  expect_error(RDAset(rangeSNPs), "set has no phenotypic information.")
})