context("Region Analysis")

library(minfiData)
miniset <- MsetEx[1:100, ]
set <- prepareMethylationSet(miniset, pData(MsetEx))
pData(set) <- preparePhenotype(pData(set), variable_names = c("status", "age"),
                                   variable_types = c("categorical","continuous"))
model <- createModel(pData(set))

test_that("Region Analysis works", {
  regions <- DARegion(set = set, model = model, methods = c("DMRcate", "bumphunter", "blockFinder"))
  expect_match(class(regions), "list")
  expect_equal(length(regions), 3)
  expect_equal(names(regions), c("bumphunter", "blockFinder", "DMRcate"))
})

test_that("Multiple variables", {
  regions <- DARegion(set = set, model = model, coefficient = 2:3)
  expect_match(class(regions), "list")
  expect_equal(length(regions), 2)
  expect_equal(names(regions), c("statusnormal", "age"))
  expect_equal(length(regions[[1]]), 3)
  expect_equal(names(regions[[1]]), c("bumphunter", "blockFinder", "DMRcate"))
})

test_that("Check Permutations", {
  regions <- DARegion(set = set, model = model, num_permutations = 10, verbose = TRUE, methods = c("DMRcate", "bumphunter", "blockFinder"))
  expect_match(class(regions), "list")
  expect_equal(length(regions), 3)
  expect_equal(names(regions), c("bumphunter", "blockFinder", "DMRcate"))
})

test_that("Select methods", {
  regions <- DARegion(set = set, model = model, methods = "bumphunter")
  expect_match(class(regions), "list")
  expect_equal(length(regions), 3)
  expect_equal(names(regions), c("bumphunter", "blockFinder", "DMRcate"))
  expect_equal(regions$blockFinder, NA)
  expect_equal(regions$DMRcate, NA)
  regions <- DARegion(set = set, model = model, methods = "blockFinder")
  expect_equal(names(regions), c("bumphunter", "blockFinder", "DMRcate"))
  expect_equal(regions$DMRcate, NA)
  regions <- DARegion(set = set, model = model, methods = "DMRcate")
  expect_equal(names(regions), c("bumphunter", "blockFinder", "DMRcate"))
  expect_equal(regions$bumphunter, NA)
  regions <- DARegion(set = set, model = model, methods = c("DMRcate", "bumphunter"))
  expect_equal(names(regions), c("bumphunter", "blockFinder", "DMRcate"))
  expect_equal(regions$blockFinder, NA)
  expect_error(DARegion(set = set, model = model, methods = character()), 
               "Method variable is empty or none of the methods introduced is valid.")
  expect_error(DARegion(set = set, model = model, methods = "character"), 
               "Method variable is empty or none of the methods introduced is valid.")
  regions <- DARegion(set = set, model = model, methods = c("DMRcate", "character"))
  expect_equal(names(regions), c("bumphunter", "blockFinder", "DMRcate"))
  expect_equal(regions$bumphunter, NA)
})

test_that("Empty variables", {
  emptyset <- new(Class = "MethylationSet")
  emptymodel <-  matrix(ncol = 0, nrow = 0)
  expect_error(DARegion(set = emptyset, model = emptymodel), "The set is empty.")
  expect_error(DARegion(set = set, model = emptymodel), "The model matrix is empty.")
  expect_error(DARegion(set = emptyset, model = model), "The set is empty.")
})

test_that("Check Errors", {
  model2 <- model[1:3, ]
  expect_error(DARegion(set = set, model = model2), "The number of samples is different in the set and in the model.")
  eset <- ExpressionSet(matrix(runif(12, max = 15), 3))
  expect_error(DARegion(set = eset, model = model), "set must be a MethylationSet.")
  emptyset <- new(Class = "MethylationSet")
  expect_error(DARegion(set = emptyset, model = model), "The set is empty.")
})
