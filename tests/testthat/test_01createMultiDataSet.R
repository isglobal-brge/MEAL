context("Create MultiSet")

test_that("MultiSet", {
  multi <- new("MultiDataSet")
  eset <- new("ExpressionSet", exprs = matrix(runif(4), ncol = 2))
  fData(eset) <- data.frame(chromosome = c("chr1", "chr2"), start = c(12414, 1234321),
                            end = c(121241, 12412412414), stringsAsFactors = FALSE)
  multi <- add.genexp(multi, eset)
  expect_equal(names(multi), "expression")
  
  multiset <- new("MultiDataSet")
  geno <- matrix(c(3,1,2,1), ncol = 2)
  colnames(geno) <- c("VAL0156", "VAL0372")
  rownames(geno) <- c("rs3115860", "SNP1-1628854")
  map <- AnnotatedDataFrame(data.frame(chromosome = c("chr1", "chr2"), position = c(12414, 1234321),
                                       stringsAsFactors = FALSE))
  rownames(map) <- rownames(geno)
  snpSet <- new("SnpSet", call = geno, featureData = map)
  
  beta_matrix <- matrix(runif(10), nrow = 5)
  rownames(beta_matrix) <- c("cg00050873", "cg00212031", "cg00213748", "cg00214611", "cg00455876")
  colnames(beta_matrix) <- c("VAL0156", "VAL0372")
  phenotypes <- data.frame(age = c(12, 23))
  rownames(phenotypes) <- c("VAL0156", "VAL0372")
  set <- prepareMethylationSet(matrix = beta_matrix, phenotypes = phenotypes)
  
  multiset <- add.methy(multiset, set)
  expect_equal(names(multiset), "methylation")
  multiset <- add.snps(multiset, snpSet)
  expect_equal(names(multiset), c("methylation", "snps"))
  
  expect_is(multiset[["methylation"]], "MethylationSet")
  expect_is(multiset[["snps"]], "SnpSet")
  
  library(minfiData)
  sampleNames(MsetEx)[1:2] <- sampleNames(set)
  multiset <- add.set(multiset, MsetEx, "cot")
  expect_equal(names(multiset), c("methylation", "snps", "cot"))
  expect_equal(names(multiset[c("cot", "snps")]), c("cot", "snps"))
  expect_is(multiset["snps"], "SnpSet")
  multiset[ , sampleNames(set)[1]]
  
})

test_that("wrong variables",{
  multiset <- new("MultiDataSet")
  
  geno <- matrix(c(3,1,2,1), ncol = 2)
  snpSet <- new("SnpSet", call = geno)
  expect_error(add.snps(multiset, snpSet), "snpSet must contain a fData with columns chromosome and position.")
  fData(snpSet) <- data.frame(chr = c("chr1", "chr2"), pos = c(12414, 1234321),
                                       stringsAsFactors = FALSE)
  expect_error(add.snps(multiset, snpSet), "snpSet must contain a fData with columns chromosome and position.")
  fData(snpSet) <- data.frame(chromosome = c("chr1", "chr2"), pos = c(12414, 1234321),
                              stringsAsFactors = FALSE)
  expect_error(add.snps(multiset, snpSet), "snpSet must contain a fData with columns chromosome and position.")
  fData(snpSet) <- data.frame(chr = c("chr1", "chr2"), position = c(12414, 1234321),
                              stringsAsFactors = FALSE)
  expect_error(add.snps(multiset, snpSet), "snpSet must contain a fData with columns chromosome and position.")
  
  
  eset <- new("ExpressionSet", exprs = matrix(runif(4), ncol = 2))
  expect_error(add.genexp(multiset, eset), "gexpSet must contain a fData with columns chromosome, start and end.")
  fData(eset) <- data.frame(chr = c("chr1", "chr2"), pos = c(12414, 1234321),
                              stringsAsFactors = FALSE)
  expect_error(add.genexp(multiset, eset), "gexpSet must contain a fData with columns chromosome, start and end.")
  fData(eset) <- data.frame(chromosome = c("chr1", "chr2"), pos = c(12414, 1234321),
                              stringsAsFactors = FALSE)
  expect_error(add.genexp(multiset, eset), "gexpSet must contain a fData with columns chromosome, start and end.")
  fData(eset) <- data.frame(chr = c("chr1", "chr2"), position = c(12414, 1234321),
                              stringsAsFactors = FALSE)
  expect_error(add.genexp(multiset, eset), "gexpSet must contain a fData with columns chromosome, start and end.")
})