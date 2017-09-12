context("Snps p-vals")
 
library(minfiData)

set <- prepareMethylationSet(getBeta(MsetEx)[1:10,], pheno = data.frame(pData(MsetEx)))

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
snps <- new("SnpSet", call = geno)
fData(snps) <- map

test_that("p-vals calculation works", {
  pvals <- calculateRelevantSNPs(set, snps)
  expect_equal(nrow(pvals), 1)
  expect_equal(ncol(pvals), as.integer(nrow(set)))
})

multi <- createMultiDataSet()
multi <- add_methy(multi, set)
multi <- add_snps(multi, snps)


test_that("Correlation Function", {
  res <- correlationMethSNPs(multi, range = range, variable_names = "sex", 
                             snps_cutoff = 0.5)
})


test_that("Empty variables", {
  emptyset <- new(Class = "MethylationSet")
  emptySnps <- new(Class = "SnpSet")
  expect_error(calculateRelevantSNPs(emptyset, snps), "set must contain beta values.")
  expect_error(calculateRelevantSNPs(set, emptySnps), "SnpSet is empty.")
  library(minfiData)
  expect_error(calculateRelevantSNPs(MsetEx, snps), "set must be a MethylationSet")
  
})

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

