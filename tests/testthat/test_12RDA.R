context("RDA calculation")

set.seed(0)
library(minfiData)
set <- prepareMethylationSet(getBeta(MsetEx)[1:10,], pheno = pData(MsetEx))

range <- GenomicRanges::GRanges(seqnames=Rle("chrY"), 
                                ranges = IRanges(3000000, end=12300000))
set <- filterSet(set, range)
vars <- DAPipeline(set, variable_names = "sex")
eqset <- DAPipeline(set, variable_names = "sex", covariable_names = "status", equation = "~sex + status", num_var = 1)

test_that("RDA calculation works", {
  rda <- RDAset(vars)
  expect_equal(class(rda), c("rda", "cca"))
  rda <- RDAset(eqset)
  expect_equal(class(rda), c("rda", "cca"))
})

test_that("wrong variables", {
  expect_error(RDAset(set), "set must be a AnalysisResults")
  
  emptyset <- new(Class = "AnalysisRegionResults")
  expect_error(RDAset(emptyset), "set has no beta values.")
  
  pData(vars) <- data.frame()
  expect_error(RDAset(vars), "set has no phenotypic information.")
})