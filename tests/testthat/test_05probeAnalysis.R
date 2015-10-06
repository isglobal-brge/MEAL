context("Probe Analysis")

library(minfiData)
miniset <- MsetEx[1:10, ]
set <- prepareMethylationSet(miniset, pData(MsetEx))
pData(set) <- preparePhenotype(pData(set), variable_names = c("status", "age"),
                                   variable_types = c("categorical","continuous"))
model <- createModel(pData(set))

test_that("Probe Analysis works", {
  probes <- DAProbe(set = set, model = model)
  expect_match(class(probes), "data.frame")
  expect_equal(as.integer(nrow(set)), nrow(probes))
  beta <- betas(set)
  probes <- DAProbe(set = beta, model = model)
  expect_match(class(probes), "data.frame")
  expect_equal(as.integer(nrow(beta)), nrow(probes))
})

test_that("Multiple Variables", {
  probes <- DAProbe(set = set, model = model, coefficient = 2:3)
  expect_match(class(probes), "list")
  expect_equal(names(probes), c("statusnormal", "age"))
  expect_match(class(probes[[1]]), "data.frame")
  expect_equal(as.integer(nrow(set)), nrow(probes[[1]]))
})

test_that("Empty variables", {
  emptyset <- new(Class = "MethylationSet")
  emptyeset <- new(Class = "ExpressionSet")
  emptymodel <-  matrix(ncol = 0, nrow = 0)
  expect_error(DAProbe(set = emptyset, model = emptymodel), "The set is empty.")
  expect_error(DAProbe(set = set, model = emptymodel), "The model matrix is empty.")
  expect_error(DAProbe(set = emptyset, model = model), "The set is empty.")
  expect_error(DAProbe(set = emptyeset, model = model), "The set is empty.")
})

test_that("Wrong variables", {
  expect_error(DAProbe(set = set, model= model[1:2, ]), 
               "The number of samples is different in the set and in the model")
  expect_error(DAProbe(set = set, model= model, method = character()), 
               "method should be one of \"ls\", \"robust\"")
  expect_error(DAProbe(set = set, model= model, method = "character"), 
               "method should be one of \"ls\", \"robust\"")
  expect_error(DAProbe(set = set, model= model, method = c("ls", "robust")), 
               "method should be one of \"ls\", \"robust\"")
  expect_warning(DAProbe(set = set, model= model, max_iterations = -1), 
               "max_iterations must be a numeric greater than 1. The default value will be used.")
  expect_warning(DAProbe(set = set, model= model, max_iterations = "1"), 
                 "max_iterations must be a numeric greater than 1. The default value will be used.")
  expect_warning(DAProbe(set = set, model= model, max_iterations = numeric()), 
                 "max_iterations must be a numeric greater than 1. The default value will be used.")
})

eset <- ExpressionSet(matrix(runif(12, max = 15), 3))
annot <- data.frame(chr = rep("chrY", 3), start = c(10000, 4500000, 5000000), 
                    end = c(100090, 6000000, 90000000))
fData(eset) <- annot
pData(eset) <- data.frame(sex = c("H", "M", "H", "M"))

model <- createModel(pData(eset))

test_that("Probe Analysis works with ExpressionSet", {
  probes <- DAProbe(set = eset, model = model)
  expect_match(class(probes), "data.frame")
  expect_equal(as.integer(nrow(eset)), nrow(probes))
})