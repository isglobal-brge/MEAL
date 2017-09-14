context("Pipeline")

library(minfiData)
# miniset <- MsetEx[1:100, ]
# pheno <- data.frame(pData(MsetEx))
# pheno$geno <- c("aa", "AA", "aA", "AA", "AA", "aA")
# set <- prepareMethylationSet(miniset, pheno)

set <- mapToGenome(ratioConvert(MsetEx[1:10, ]), mergeManifest = TRUE)
set <- set[apply(getBeta(set), 1, function(x) all(!is.na(x))), ]

eset <- ExpressionSet(matrix(runif(400, max = 15), 100))
annot <- data.frame(chr = rep("chr1", 100), start = 1:100, 
                    end = 1001:1100)
fData(eset) <- annot
pData(eset) <- data.frame(sex = rep(c("H", "M", "H", "M")))

range <- GRanges("chr1:1-2100")

test_that("run Pipeline", {
  results <- runPipeline(set = set, variable_names = "sex")
  expect_match(class(results), "ResultSet")
  results <- runPipeline(set = eset, variable_names = "sex", verbose = TRUE)
  expect_match(class(results), "ResultSet")

  # results <- runPipeline(set = eset, variable_names = "sex", sva = TRUE)
  expect_match(class(results), "ResultSet")
  
  results <- runPipeline(set = eset, variable_names = "sex", range = range)
  expect_match(class(results), "ResultSet")
  
})
