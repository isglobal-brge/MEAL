# context("Filtering")
# 
# 
# library(minfiData)
# 
# range <- GenomicRanges::GRanges(seqnames=Rle("chrY"), 
#                                 ranges = IRanges(3000000, end = 12300000))
# set <- prepareMethylationSet(MsetEx[1:100, ], data.frame(pData(MsetEx)))
# 
# eset <- ExpressionSet(matrix(runif(6, max = 15), 4))
# annot <- data.frame(chromosome = c(rep("chrY", 3), "chr1" ), start = c(10000, 4500000, 5000000, 10000), 
#                     end = c(100090, 6000000, 90000000, 100090))
# fData(eset) <- annot
# 
# 
# test_that("Results filtering", {
#     methyOneVar <- DAPipeline(set, variable_names = "sex", probe_method = "ls")
#     
#   filteredResults <- filterResults(probeResults(methyOneVar), range)
#   expect_match(class(filteredResults), "data.frame")
#   expect_equal(nrow(filteredResults), 32)
#   filteredResults <- filterResults(probeResults(methyOneVar), range)
#   expect_match(class(filteredResults), "data.frame")
#   expect_equal(nrow(filteredResults), 32)
# })
# 
# test_that("Results filtering, empty variables", {
#     methyOneVar <- DAPipeline(set, variable_names = "sex", probe_method = "ls")
#     
#   emptyRes <- data.frame()
#   emptyrange <- GRanges()
#   expect_error(filterResults(emptyRes, emptyrange), "results is empty")
#   expect_error(filterResults(probeResults(methyOneVar), emptyrange), "range is empty")
# })
# 
# test_that("Results filtering, wrong variables", {
#     methyOneVar <- DAPipeline(set, variable_names = "sex", probe_method = "ls")
#     
#   outrange <- GenomicRanges::GRanges(seqnames=Rle("chrM"), 
#                                      ranges = IRanges(30000, end=123000000))
#   expect_warning(filterResults(probeResults(methyOneVar), outrange), "There are no cpgs in the range. An empty data.frame will be returned.")
#   outrange <- data.frame(chr = "chr1", start = 19051, end = 1251905)
#   expect_error(filterResults(probeResults(methyOneVar), outrange), "range must be a GenomicRanges object.")
#   outrange <- GenomicRanges::GRanges(seqnames=Rle(c("chr1", "chr2")), 
#                                      ranges = IRanges(c(30000, 2149421), end=c(123000000, 234841418)))
#   expect_error(filterResults(probeResults(methyOneVar), outrange), "range must be a GenomicRanges with only one range.")
# }