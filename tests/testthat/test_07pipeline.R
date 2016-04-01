context("Pipeline")

library(minfiData)
miniset <- MsetEx[1:10, ]
pheno <- pData(MsetEx)
pheno$geno <- c("aa", "AA", "aA", "AA", "AA", "aA")
set <- prepareMethylationSet(miniset, pheno)

eset <- ExpressionSet(matrix(runif(12, max = 15), 3))
annot <- data.frame(chromosome = rep("chrY", 3), start = c(10000, 4500000, 5000000), 
                    end = c(100090, 6000000, 90000000))
fData(eset) <- annot
pData(eset) <- data.frame(sex = c("H", "M", "H", "M"))

test_that("One variable categorical, two levels", {
  results <- DAPipeline(set = set, variable_names = "sex", variable_types = "categorical")
  expect_match(class(results), "AnalysisResults")
  results <- DAPipeline(set = eset, variable_names = "sex", variable_types = "categorical",
                            verbose = TRUE)
  expect_match(class(results), "AnalysisResults")
  results <- DAPipeline(set = set, variable_names = "sex", variable_types = "categorical", sva = TRUE)
  expect_match(class(results), "AnalysisResults")
})

test_that("One variable categorical, three levels", {
  results <- DAPipeline(set = set, variable_names = "person", variable_types = "categorical")
  expect_match(class(results), "AnalysisResults")
  expect_equal(length(probeResults(results)), 2)
  expect_equal(names(probeResults(results)), c("personid2", "personid3"))
})

test_that("One variable categorical, two levels, one covariate", {
  results <- DAPipeline(set = set, variable_names = "sex", variable_types = "categorical", 
                            covariable_names = "status", covariable_types = "categorical")
  expect_match(class(results), "AnalysisResults")
})

test_that("Two categorical variables", {
  results <- DAPipeline(set = set, variable_names = c("sex", "status"), 
                            variable_types = c("categorical", "categorical"))
  expect_match(class(results), "AnalysisResults")
  expect_equal(names(probeResults(results)), c("sexM", "statusnormal"))
})

test_that("One genetic variable", {
  results <- DAPipeline(set = set, variable_names = "geno", variable_types = "additive")
  expect_match(class(results), "AnalysisResults")
  respheno <- pData(results)
  expect_match(class(respheno[ , 1]), "factor")
})

test_that("Empty variables", {
  emptyset <- new(Class = "MethylationSet")
  expect_error(DAPipeline(set = emptyset, variable_names = character(), 
                              variable_types = character()), "The set has no beta values")
  expect_error(DAPipeline(set = set, variable_names = character(), 
                              variable_types = character()), "variable_names is empty.")
  expect_error(DAPipeline(set = set, variable_names = "sex", 
                              variable_types = character()), "variable_types is empty.")
  expect_error(DAPipeline(set = set, variable_names = character(), covariable_names = "sex", covariable_types = "categorical"),
               "variable_names is empty or is not a valid column of the phenoData of the set.")
  expect_error(DAPipeline(set = set, variable_names = "cot", covariable_names = "sex", covariable_types = "categorical"),
               "variable_names is empty or is not a valid column of the phenoData of the set.")
  dataset <- matrix(runif(20), 5)
  expect_error(DAPipeline(set = dataset, variable_names = character(), 
                              variable_types = character()), "set must be a MethylationSet.")
  
  })