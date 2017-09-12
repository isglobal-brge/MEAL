context("Probe Analysis")

library(minfiData)
miniset <- MsetEx[1:10, ]
set <- prepareMethylationSet(miniset, data.frame(pData(MsetEx)))
# pData(set) <- preparePhenotype(pData(set), variable_names = c("status", "age"),
#                                    variable_types = c("categorical","continuous"))
# model <- createModel(pData(set))

test_that("DiffMean", {
  probes <- DiffMeanAnalysis(set = set, model = ~ status)
  expect_match(class(probes), "ResultSet")
  
  probes <- DiffMeanAnalysis(set = set, model = ~ status, resultSet = FALSE)
  expect_match(class(probes), "MArrayLM")
  expect_equal(as.integer(nrow(set)), nrow(probes))
  
  beta <- betas(set)
  model <- model.matrix(~ status, pData(set))
  probes <- DiffMeanAnalysis(set = beta, model = model)
  expect_match(class(probes), "ResultSet")
})

test_that("DiffVar", {
  probes <- DiffVarAnalysis(set = set, model = ~ status)
  expect_match(class(probes), "ResultSet")
  
  probes <- DiffVarAnalysis(set = set, model = ~ status, resultSet = FALSE)
  expect_match(class(probes), "MArrayLM")
  expect_equal(as.integer(nrow(set)), nrow(probes))
  
  beta <- betas(set)
  model <- model.matrix(~ status, pData(set))
  probes <- DiffVarAnalysis(set = beta, model = model)
  expect_match(class(probes), "ResultSet")
})


eset <- ExpressionSet(matrix(runif(400, max = 15), 4))
annot <- data.frame(chr = rep("chr1", 4), start = c(10, 1000, 3000, 3100), 
                    end = c(900, 2000, 3300, 4000))
fData(eset) <- annot
pData(eset) <- data.frame(sex = rep(c("H", "M", "H", "M"), 25))
rownames(eset) <- letters[1:4]

test_that("Probe Analysis works with ExpressionSet", {
  probes <- DiffMeanAnalysis(set = eset, model = ~sex)
  expect_match(class(probes), "ResultSet")
  
  probes <- DiffMeanAnalysis(set = eset, model = ~ sex, resultSet = FALSE)
  expect_match(class(probes), "MArrayLM")
  expect_equal(as.integer(nrow(eset)), nrow(probes))
  
  probes <- DiffVarAnalysis(set = eset, model = ~sex)
  expect_match(class(probes), "ResultSet")
  
  probes <- DiffVarAnalysis(set = eset, model = ~ sex, resultSet = FALSE)
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

