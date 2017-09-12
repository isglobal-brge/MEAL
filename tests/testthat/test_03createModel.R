context("Create Model")

test_that("Check colnames", {
  pheno <- data.frame(a = c("a", "b", "c", "b", "a"), d = 1:5, e = runif(5))
  model <- createModel(pheno)
  expect_equal(ncol(model), 5)
  names <- letters[1:5]
  model <- createModel(pheno, names = names)
  expect_equal(colnames(model), names)
})


test_that("Check equation", {
  pheno <- data.frame(a = c("a", "b", "a", "b", "a"), d = 1:5, e = runif(5))
  model <- createModel(pheno, equation = "~ a")
  expect_equal(ncol(model), 4)
  model <- createModel(pheno, equation = "~ a + d")
  expect_equal(ncol(model), 4)
  model <- createModel(pheno, equation = "~ a + d + a:d")
  expect_equal(ncol(model), 5)
  expect_equal(colnames(model), c("(Intercept)", "ab", "d", "e", "ab:d"))
  names <- letters[1:5]
  model <- createModel(pheno, equation = "~ a + d + a:d", names = names)
  expect_equal(colnames(model), names)
})

test_that("Check wrong variables", {
  pheno <- data.frame()
  expect_error(createModel(pheno), "data is empty.")
  pheno <- 1:5
  expect_error(createModel(pheno), "data is not a valid data.frame")
  pheno <- data.frame(a = c("a", "b", "a", "b", "a"), d = 1:5, e = runif(5))
  equation <- character()
  expect_error(createModel(pheno, equation), "Equation is an invalid formula.")
  equation <- "a"
  expect_error(createModel(pheno, equation), "Equation is an invalid formula.")
  equation <- "~ A"
  expect_error(createModel(pheno, equation), "A model can't be created with this equation.")
  expect_warning(createModel(pheno, names = character()), 
                 "names variable is empty.")
  expect_warning(createModel(pheno, equation = "~ a + d", names = character()), 
                 "names variable is empty.")
  expect_warning(createModel(pheno, names = "a"), 
                 "Length of names is different than column names of model. Names won't be changed.")
  expect_warning(createModel(pheno, equation = "~ a + d", names = "a"), 
                 "Length of names is different than column names of model. Names won't be changed.")
})


