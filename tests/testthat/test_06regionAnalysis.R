context("Region Analysis")

library(minfiData)
library(minfi)
miniset <- MsetEx[1:1000, ]
set <- ratioConvert(mapToGenome(miniset))

test_that("Bumphunter", {
  regions <- runBumphunter(set = set, model = ~ status)
  expect_match(class(regions), "data.frame")
})

test_that("BlockFinder", {
  regions <- runBlockFinder(set = set, model = ~ status)
  expect_match(class(regions), "data.frame")
})

test_that("DMRcate", {
  regions <- runDMRcate(set = set, model = ~ status, pcutoff = 0.9)
  expect_match(class(regions), "data.frame")
  expect_equal(nrow(regions), 9)
})

test_that("Region Analysis Wrapper", {
  regions <- runRegionAnalysis(set = set, model = ~ status)
  expect_match(class(regions), "ResultSet")
  
  regions <- runRegionAnalysis(set = set, model = ~ status, 
                            bumphunter_params = list(bumphunter_cutoff = 0.01), 
                            blockFinder_params = list(blockfinder_cutoff = 0.01), 
                            dmrcate_params = list(pcutoff = 0.75))
  expect_match(class(regions), "ResultSet")
})