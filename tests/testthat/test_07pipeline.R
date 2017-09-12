context("Pipeline")

library(minfiData)
library(BRGEdata)
# miniset <- MsetEx[1:100, ]
# pheno <- data.frame(pData(MsetEx))
# pheno$geno <- c("aa", "AA", "aA", "AA", "AA", "aA")
# set <- prepareMethylationSet(miniset, pheno)

set <- mapToGenome(ratioConvert(MsetEx), mergeManifest = TRUE)
set <- set[apply(getBeta(set), 1, function(x) all(!is.na(x))), ]

eset <- ExpressionSet(matrix(runif(400, max = 15), 100))
annot <- data.frame(chr = rep("chr1", 100), start = 1:100, 
                    end = 1001:1100)
fData(eset) <- annot
pData(eset) <- data.frame(sex = rep(c("H", "M", "H", "M")))

range2 <- GRanges("chr1:1-2100")

test_that("run Pipeline", {
  results <- runPipeline(set = set, variable_names = "sex")
  expect_match(class(results), "ResultSet")
  results <- runPipeline(set = eset, variable_names = "sex", verbose = TRUE)
  expect_match(class(results), "ResultSet")

  results <- runPipeline(set = eset, variable_names = "sex", sva = TRUE)
  expect_match(class(results), "ResultSet")
  
  results <- runPipeline(set = eset, variable_names = "sex", range = range)
  expect_match(class(results), "ResultSet")
  
})

# test_that("One variable categorical, three levels", {
#   results <- DAPipeline(set = set, variable_names = "person", variable_types = "categorical")
#   expect_match(class(results), "AnalysisResults")
#   expect_equal(length(probeResults(results)), 2)
#   expect_equal(names(probeResults(results)), c("personid2", "personid3"))
# })
# 
# test_that("One variable categorical, two levels, one covariate", {
#   results <- DAPipeline(set = set, variable_names = "sex", variable_types = "categorical", 
#                             covariable_names = "status", covariable_types = "categorical")
#   expect_match(class(results), "AnalysisResults")
# })
# 
# test_that("Two categorical variables", {
#   results <- DAPipeline(set = set, variable_names = c("sex", "status"), 
#                             variable_types = c("categorical", "categorical"))
#   expect_match(class(results), "AnalysisResults")
#   expect_equal(names(probeResults(results)), c("sexM", "statusnormal"))
# })
# 
# test_that("One genetic variable", {
#   results <- DAPipeline(set = set, variable_names = "geno", variable_types = "additive")
#   expect_match(class(results), "AnalysisResults")
#   respheno <- pData(results)
#   expect_match(class(respheno[ , 1]), "factor")
# })
# 
# test_that("Empty variables", {
#   emptyset <- new(Class = "MethylationSet")
#   expect_error(DAPipeline(set = emptyset, variable_names = character(), 
#                               variable_types = character()), "The set has no beta values")
#   expect_error(DAPipeline(set = set, variable_names = character(), 
#                               variable_types = character()), "variable_names is empty.")
#   expect_error(DAPipeline(set = set, variable_names = "sex", 
#                               variable_types = character()), "variable_types is empty.")
#   expect_error(DAPipeline(set = set, variable_names = character(), covariable_names = "sex", covariable_types = "categorical"),
#                "variable_names is empty or is not a valid column of the phenoData of the set.")
#   expect_error(DAPipeline(set = set, variable_names = "cot", covariable_names = "sex", covariable_types = "categorical"),
#                "variable_names is empty or is not a valid column of the phenoData of the set.")
#   dataset <- matrix(runif(20), 5)
#   expect_error(DAPipeline(set = dataset, variable_names = character(), 
#                               variable_types = character()), "set must be a MethylationSet.")
#   
#   })