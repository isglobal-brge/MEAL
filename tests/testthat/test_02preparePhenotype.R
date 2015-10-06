context("Prepare Phenotype")

test_that("Perform right filtering", {
  pheno <- data.frame(a = c("a", "b", "a", "b", "a"), b = 1:5, c = runif(5))
  expect_equal(ncol(preparePhenotype(pheno, variable_names = "a", variable_types = NA)), 1)
  expect_equal(ncol(preparePhenotype(pheno, variable_names = c("a", "b"))), 2)  
  expect_equal(ncol(preparePhenotype(pheno, variable_names = 1)), 1)  
  expect_equal(ncol(preparePhenotype(pheno, variable_names = 1:2)), 2)    
})

test_that("Perform right transformation", {
  pheno <- data.frame(a = c("a", "b", "a", "b", "a"), b = 1:5,
                      c = c("a/a", "a/b", "a/a", "b/b", "a/a"))
  pheno1 <- preparePhenotype(pheno, variable_names = "a", variable_types = "categorical")
  expect_match(class(pheno1[,1]), "factor")
  pheno1 <- preparePhenotype(pheno, variable_names = "b", variable_types = "continuous")
  expect_match(class(pheno1[,1]), "numeric")
  pheno1 <- preparePhenotype(pheno, variable_names = "c", variable_types = "dominant")
  expect_match(class(pheno1[,1]), "factor")
  expect_equal(nlevels(pheno1[,1]), 2)
  expect_equal(levels(pheno1[,1]), c("a/a", "a/b-b/b"))
  pheno1 <- preparePhenotype(pheno, variable_names = "c", variable_types = "additive")
  expect_match(class(pheno1[,1]), "numeric")
  pheno1 <- preparePhenotype(pheno, variable_names = "c", variable_types = "recessive")
  expect_match(class(pheno1[,1]), "factor")
  expect_equal(nlevels(pheno1[,1]), 2)
  expect_equal(levels(pheno1[,1]), c("a/a-a/b", "b/b"))
})

test_that("Wrong variables", {
  pheno <- data.frame()
  expect_error(preparePhenotype(pheno, variable_names = "a", variable_types = NA), 
               "Phenotypes data.frame is empty.")
  pheno <- runif(5)
  expect_error(preparePhenotype(pheno, variable_names = "a", variable_types = NA), 
               "Phenotypes is not a valid data.frame")
  pheno <- data.frame(a = c("a", "b", "a", "b", "a"), b = 1:5,
                      c = c("a/a", "a/b", "a/a", "b/b", "a/a"))  
  expect_error(preparePhenotype(pheno, variable_names = data.frame(a = c("a", "b")), 
                                variable_types = NA), "variable_names must be a vector.")
  expect_error(preparePhenotype(pheno, variable_names = "a", variable_types = data.frame(a = c("a", "b"))), 
               "variable_types must be a vector.")
  expect_error(preparePhenotype(pheno, variable_names = character(), variable_types = character()), 
               "variable_names is empty.")
  expect_error(preparePhenotype(pheno, variable_names = "a", variable_types = character()), 
               "variable_types is empty.")  
})

test_that("Wrong variable names", {
  pheno <- data.frame(a = c("a", "b", "a", "b", "a"), b = 1:5,
                      c = c("a/a", "a/b", "a/a", "b/b", "a/a"))  
  expect_error(preparePhenotype(pheno, variable_names = "A", variable_types = NA), 
               "A are not columns of the phenotypes table.")
  expect_error(preparePhenotype(pheno, variable_names = c("A", "a"), variable_types = NA), 
               "A are not columns of the phenotypes table.")
  expect_error(preparePhenotype(pheno, variable_names = c("A", "B"), variable_types = NA), 
               "A, B are not columns of the phenotypes table.")
  expect_error(preparePhenotype(pheno, variable_names = 0, variable_types = NA), 
               "variable_names indexes must be between 1 and 3")
  expect_error(preparePhenotype(pheno, variable_names = c(0, 1), variable_types = NA), 
               "variable_names indexes must be between 1 and 3")
  expect_error(preparePhenotype(pheno, variable_names = c(0, 6), variable_types = NA), 
               "variable_names indexes must be between 1 and 3")
  expect_error(preparePhenotype(pheno, variable_names = c(1, 6), variable_types = NA), 
               "variable_names indexes must be between 1 and 3")
  expect_error(preparePhenotype(pheno, variable_names = 6, variable_types = NA), 
               "variable_names indexes must be between 1 and 3")
 })


test_that("Wrong variable types", {
  pheno <- data.frame(a = c("a", "b", "a", "b", "a"), b = 1:5,
                      c = c("a/a", "a/b", "a/a", "b/b", "a/a"))  
  expect_warning(preparePhenotype(pheno, variable_names = c("a", "b"), variable_types = "categorical"), 
               "The number of types of variables is smaller than the number of variables. Variables without type won't be changed.")
  })