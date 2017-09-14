context("RDA calculation")

library(minfiData)
set <- mapToGenome(ratioConvert(MsetEx[1:10, ]), mergeManifest = TRUE)

range <- GenomicRanges::GRanges(seqnames=Rle("chrY"), 
                                ranges = IRanges(3000000, end=12300000))
test_that("runRDA", {
  rda <- runRDA(set = set, ~ sex)
  expect_true(is(rda, "ResultSet"))
  # rda2 <- runRDA(set = set, ~ sex + status, num_vars = 2)
  # expect_equal(class(rda), c("rda", "cca"))
  # rda <- runRDA(set = set, ~ sex, range = range)
  # expect_equal(class(rda), c("rda", "cca"))
})