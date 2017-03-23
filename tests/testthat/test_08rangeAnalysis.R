context("Range Analysis")

library(minfiData)
set <- prepareMethylationSet(getBeta(MsetEx)[1:10,], pheno = pData(MsetEx))
range <- GenomicRanges::GRanges(seqnames=Rle("chrY"), 
                                ranges = IRanges(3000000, end=12300000))

eset <- ExpressionSet(matrix(runif(12, max = 15), 3))
annot <- data.frame(chromosome = rep("chrY", 3), start = c(10000, 4500000, 5000000), 
                    end = c(100090, 6000000, 90000000))
fData(eset) <- annot
pData(eset) <- data.frame(sex = c("H", "M", "H", "M"))


test_that("One variable categorical, two levels", {
  results <- DARegionAnalysis(set = set, variable_names = "sex", 
                                  variable_types = "categorical", range = range)
  expect_match(class(results), "AnalysisRegionResults")
 results <- DARegionAnalysis(set = eset, variable_names = "sex", 
                                  variable_types = "categorical", range = range)
  expect_match(class(results), "AnalysisRegionResults")
  results <- DARegionAnalysis(set = eset, variable_names = "sex", sva = TRUE,
                                  variable_types = "categorical", range = range)
  expect_match(class(results), "AnalysisRegionResults")
  results <- DARegionAnalysis(set = set, variable_names = "sex", sva = TRUE,
                                  variable_types = "categorical", range = range)
  expect_match(class(results), "AnalysisRegionResults")
})

test_that("One genetic variable", {
  miniset <- MsetEx[1:10, ]
  pheno <- pData(MsetEx)
  pheno$geno <- c("aa", "AA", "aA", "AA", "AA", "aA")
  set2 <- prepareMethylationSet(miniset, pheno)
  results <- DARegionAnalysis(set = set2, variable_names = "geno", 
                                  variable_types = "additive", range = range)
  expect_match(class(results), "AnalysisRegionResults")
  respheno <- pData(results)
  expect_match(class(respheno[ , 1]), "factor")
})

test_that("Empty variables", {
  emptyset <- new(Class = "MethylationSet")
  expect_error(DARegionAnalysis(set = emptyset, variable_names = character(), 
                                    variable_types = character(), range = range), "The set has no beta values")
   dataset <- matrix(runif(20), 5)
  expect_error(DARegionAnalysis(set = dataset, variable_names = character(), 
                                    variable_types = character(), range = range),
               "set must be a MethylationSet.")
  emptyrange <- GRanges()
  expect_error(DARegionAnalysis(set = set, variable_names = "sex", 
                                    variable_types = "categorical", range = emptyrange), "range is empty")
  outrange <- GenomicRanges::GRanges(seqnames = Rle("chrM"), 
                                     ranges = IRanges(30000, end=123000000))
  expect_error(DARegionAnalysis(set = set, variable_names = "sex", 
                                    variable_types = "categorical", range = outrange), 
               "There are no features in the range specified.")
})