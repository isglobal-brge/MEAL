context("Explained Variance")

data(mtcars)

test_that("Basic Explained Variance",{
  expect_equal(length(explainedVariance(mtcars)), 10)
  expect_equal(length(explainedVariance(mtcars, num_covariates = 2)), 9)
})

test_that("Wrong variables",{
  expect_error(explainedVariance(1:4), "data must be a data.frame or a matrix.")
  expect_error(explainedVariance(data.frame()), "data must contain at lest two columns.")
  expect_error(explainedVariance(data.frame(1:4)), "data must contain at lest two columns.")
  expect_error(explainedVariance(mtcars, num_covariates = 20), 
               "Maximum number of covariates is number of columns minus number of main variables minus 1.")
  expect_error(explainedVariance(mtcars, num_covariates = -1), 
               "Number of covariates must be a positive numerical.")
  expect_error(explainedVariance(mtcars, num_mainvar = 20), 
               "Maximum number of variables is number of columns minus 1.")
  expect_error(explainedVariance(mtcars, num_mainvar = -1), 
               "Number of main variables must be a positive numerical.")
})