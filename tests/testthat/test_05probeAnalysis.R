context("Probe Analysis")

library(minfiData)
library(minfi)
miniset <- MsetEx[1:10, ]
set <- ratioConvert(mapToGenome(miniset))

test_that("DiffMean", {
  probes <- runDiffMeanAnalysis(set = set, model = ~ status)
  expect_match(class(probes), "ResultSet")
  
  probes <- runDiffMeanAnalysis(set = set, model = ~ status, resultSet = FALSE)
  expect_match(class(probes), "MArrayLM")
  expect_equal(as.integer(nrow(set)), nrow(probes))
  
  beta <- getBeta(set)
  model <- model.matrix(~ status, pData(set))
  probes <- runDiffMeanAnalysis(set = beta, model = model)
  expect_match(class(probes), "ResultSet")
})

test_that("DiffVar", {
  probes <- runDiffVarAnalysis(set = set, model = ~ status)
  expect_match(class(probes), "ResultSet")
  
  probes <- runDiffVarAnalysis(set = set, model = ~ status, resultSet = FALSE)
  expect_match(class(probes), "MArrayLM")
  expect_equal(as.integer(nrow(set)), nrow(probes))
  
  beta <- getBeta(set)
  model <- model.matrix(~ status, pData(set))
  probes <- runDiffVarAnalysis(set = beta, model = model)
  expect_match(class(probes), "ResultSet")
})


eset <- ExpressionSet(matrix(rnorm(1000), ncol = 10))
pData(eset) <- data.frame(sex = rep(c("H", "M"), each = 5))


test_that("Probe Analysis works with ExpressionSet", {
  probes <- runDiffMeanAnalysis(set = eset, model = ~sex)
  expect_match(class(probes), "ResultSet")
  
  probes <- runDiffMeanAnalysis(set = eset, model = ~ sex, resultSet = FALSE)
  expect_match(class(probes), "MArrayLM")
  expect_equal(as.integer(nrow(eset)), nrow(probes))
  
  probes <- runDiffVarAnalysis(set = eset, model = ~sex)
  expect_match(class(probes), "ResultSet")
  
  probes <- runDiffVarAnalysis(set = eset, model = ~ sex, resultSet = FALSE)
  expect_match(class(probes), "MArrayLM")
  expect_equal(as.integer(nrow(eset)), nrow(probes))
})


# 
# test_that("Empty variables", {
#   emptyset <- new(Class = "MethylationSet")
#   emptyeset <- new(Class = "ExpressionSet")
#   emptymodel <-  matrix(ncol = 0, nrow = 0)
#   expect_error(DiffMeanAnalysis(set = emptyset, model = emptymodel), "The set is empty.")
#   expect_error(DiffMeanAnalysis(set = set, model = emptymodel), "The model matrix is empty.")
#   expect_error(DiffMeanAnalysis(set = emptyset, model = model), "The set is empty.")
#   expect_error(DiffMeanAnalysis(set = emptyeset, model = model), "The set is empty.")
# })
# 
# test_that("Wrong variables", {
#   expect_error(DiffMeanAnalysis(set = set, model= model[1:2, ]), 
#                "The number of samples is different in the set and in the model")
#   expect_error(DiffMeanAnalysis(set = set, model= model, method = character()), 
#                "method should be one of \"ls\", \"robust\"")
#   expect_error(DiffMeanAnalysis(set = set, model= model, method = "character"), 
#                "method should be one of \"ls\", \"robust\"")
#   expect_error(DiffMeanAnalysis(set = set, model= model, method = c("ls", "robust")), 
#                "method should be one of \"ls\", \"robust\"")
#   expect_warning(DiffMeanAnalysis(set = set, model= model, max_iterations = -1), 
#                "max_iterations must be a numeric greater than 1. The default value will be used.")
#   expect_warning(DiffMeanAnalysis(set = set, model= model, max_iterations = "1"), 
#                  "max_iterations must be a numeric greater than 1. The default value will be used.")
#   expect_warning(DiffMeanAnalysis(set = set, model= model, max_iterations = numeric()), 
#                  "max_iterations must be a numeric greater than 1. The default value will be used.")
# })