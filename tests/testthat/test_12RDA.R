context("RDA calculation")

library(minfiData)
set <- prepareMethylationSet(getBeta(MsetEx)[1:10,], pheno = data.frame(pData(MsetEx)))

range <- GenomicRanges::GRanges(seqnames=Rle("chrY"), 
                                ranges = IRanges(3000000, end=12300000))
set <- filterSet(set, range)
model <- model.matrix(~ sex, pData(set))
covarsmodel <- model.matrix(~ sex + status, pData(set))[, 3, drop = FALSE]

test_that("RDA calculation works", {
  rda <- RDAset(set = set, varsmodel = model)
  expect_equal(class(rda), c("rda", "cca"))
  rda <- RDAset(set = set, varsmodel = model, covarsmodel = covarsmodel)
  expect_equal(class(rda), c("rda", "cca"))
})

test_that("wrong variables", {
  # expect_error(RDAset(vars), "set must be a Methylation")
  
  # emptyset <- new(Class = "AnalysisRegionResults")
  # expect_error(RDAset(emptyset), "set has no beta values.")
  
})