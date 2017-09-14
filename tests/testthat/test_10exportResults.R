# context("Export files")
# 
# library(minfiData)
# set <- prepareMethylationSet(getBeta(MsetEx)[1:10,], pheno = data.frame(pData(MsetEx)))
# methyOneVar <- DAPipeline(set, variable_names = "sex", probe_method = "ls")
# methyTwoVar <- DAPipeline(set, variable_names =  "person", probe_method = "ls")
# dir.create("./testfolder")
# setwd("./testfolder")
# 
# test_that("One variable, export same folder", {
#   exportResults(methyOneVar)
#   csvs <- dir( pattern = ".csv")
#   expect_equal(length(csvs), 4)
#   expect_equal(csvs, c("blockFinderResults.csv","bumphunterResults.csv", "dmrCateResults.csv",
#                        "probeResults.csv"))  
# })
# 
# test_that("One variable, export different folder", {
#   exportResults(methyOneVar, dir = "./test2")
#   csvs <- dir(path = "./test2",  pattern = ".csv")
#   expect_equal(length(csvs), 4)
#   expect_equal(csvs, c("blockFinderResults.csv","bumphunterResults.csv", "dmrCateResults.csv",
#                        "probeResults.csv"))  
# })
# 
# test_that("Two variables, export different folder", {
#   exportResults(methyTwoVar, dir = "./test3")
#   dirs <- dir(path = "./test3")
#   expect_equal(length(dirs), 2)
#   expect_equal(dirs, c("personid2", "personid3"))  
# })
# 
# 
# test_that("Two variables, choose one, export different folder", {
#   exportResults(methyTwoVar, dir = "./test4", vars = "personid2")
#   csvs <- dir(path = "./test4",  pattern = ".csv")
#   expect_equal(length(csvs), 4)
#   expect_equal(csvs, c("blockFinderResults.csv","bumphunterResults.csv", "dmrCateResults.csv",
#                        "probeResults.csv"))  
# })
# 
# test_that("Wrong variables", {
#   expect_error(exportResults(methyTwoVar, vars = character()), "Vars introduced are not present in the results object.")
#   expect_error(exportResults(methyTwoVar, vars = "cot"), "Vars introduced are not present in the results object.")
# })
# 
# test_that("One variable, add prefix export different folder", {
#   exportResults(methyOneVar, dir = "./test5", prefix = "p")
#   csvs <- dir(path = "./test5",  pattern = ".csv")
#   expect_equal(length(csvs), 4)
#   expect_equal(csvs, c("pblockFinderResults.csv","pbumphunterResults.csv", "pdmrCateResults.csv",
#                        "pprobeResults.csv"))  
# })
# 
# setwd("../")
# unlink("./testfolder", recursive = TRUE)