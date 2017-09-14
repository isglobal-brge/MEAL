# context("Get genes")
# 
# library(minfiData)
# set <- prepareMethylationSet(getBeta(MsetEx)[1:10,], pheno = data.frame(pData(MsetEx)))
# methyOneVar <- DAPipeline(set, variable_names = "sex", probe_method = "ls")
# methyTwoVar <- DAPipeline(set, variable_names =  "person", probe_method = "ls")
# 
# test_that("One category", {
#   res <- getGeneVals(methyOneVar, "TSPY4") 
#   expect_equal(length(res), 1)
#   expect_equal(nrow(res[[1]]), 2)  
#   expect_equal(nrow(getGeneVals(methyOneVar, "TTTY14")[[1]]), 1)  
# })
# 
# test_that("Two categories", {
#   res <- getGeneVals(methyTwoVar, "TSPY4") 
#   expect_equal(length(res), 2)
#   expect_equal(nrow(res[[1]]), 2)  
#   expect_equal(nrow(getGeneVals(methyTwoVar, "TTTY14")[[1]]), 1)  
# })
# 
# test_that("Wrong variables", {
#  expect_warning(getGeneVals(methyTwoVar, "TS12"), "Gene name is not present in the results. A list with empty data.frames will be returned.")
#   expect_error(getGeneVals(methyTwoVar, 12), "gene must be a character vector")
#   })