context("Correlation")

library(minfiData)
set <- mapToGenome(ratioConvert(MsetEx[1:10,  1:4]), mergeManifest = TRUE)
colData(set)$geno <-  c("aa", "AA", "aA", "AA")

eset <- ExpressionSet(matrix(runif(12, max = 15), 3))
fData(eset) <- data.frame(chromosome = rep("chrY", 3), start = c(10000, 6700000, 5000000), 
                          end = c(100090, 6710000, 90000000))
featureNames(eset) <- letters[1:3]
rownames(fData(eset)) <- featureNames(eset)
pData(eset) <- data.frame(sex = c("H", "M", "H", "M"))
sampleNames(eset) <- sampleNames(set)

multi <- createMultiDataSet()
multi <- add_genexp(multi, eset)
multi <- add_methy(multi, set)

test_that("Check lineal correlation is working", {
  results <- correlationMethExprs(multiset = multi)
  expect_match(class(results), "data.frame")
  expect_equal(nrow(results), 1)  
  results <- correlationMethExprs(multiset = multi, vars_meth = "person")
  expect_equal(nrow(results), 1)  
  results <- correlationMethExprs(multiset = multi, vars_exprs = "sex")
  expect_equal(nrow(results), 1)  
  results <- correlationMethExprs(multiset = multi, vars_meth = "person",
                                  vars_exprs = "sex")
  expect_equal(nrow(results), 1)  
  results <- correlationMethExprs(multiset = multi, flank = 5000000)
  fData(eset)[] <- data.frame(chromosome = rep("chrY", 3), start = c(10000, 6700000, 5000000), 
                            end = c(100090, 7710000, 90000000))
  multi2 <- add_genexp(multi, eset, overwrite = TRUE)
  expect_warning(results <- correlationMethExprs(multiset = multi2),  
                 "There are no expression probes in the range of the cpgs. An empty data.frame will be returned.")
  expect_equal(nrow(results), 0)
})

test_that("Empty and wrong variables linear", {
  expect_error(correlationMethExprs(multiset = eset), "multiset must be a MultiDataSet")
  emptymulti <- createMultiDataSet()
  expect_error(correlationMethExprs(multiset = emptymulti), "multiset must contain meth_set_name and exprs_set_name.")

 
  multi <- add_genexp(multi, eset, overwrite = TRUE)
  expect_error(correlationMethExprs(multiset = multi, flank = -5000000), "flank must be a positive integer")
  expect_error(correlationMethExprs(multiset = multi, flank = "233"), "flank must be a positive integer")
  expect_error(correlationMethExprs(multiset = multi, flank = c(50000, 100)), "flank must be a positive integer")
  # expect_warning(correlationMethExprs(multiset = multi, vars_meth = "cot"),  
  #                "vars_meth is/are not a valid column of the mset phenoData. No residues will be computed")  
  # expect_warning(correlationMethExprs(multiset = multi, vars_exprs = "cot"),  
  #                "vars_exprs is/are not a valid column of the eset phenoData. No residues will be computed")  
})