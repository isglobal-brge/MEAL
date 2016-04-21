context("Plotting")

library(minfiData)
set <- prepareMethylationSet(getBeta(MsetEx)[1:10,], pheno = pData(MsetEx))
methyOneVar <- DAPipeline(set, variable_names = "sex", probe_method = "ls")
methyTwoVar <- DAPipeline(set, variable_names =  "person", probe_method = "ls")

range <- GenomicRanges::GRanges(seqnames=Rle("chrY"), 
                                ranges = IRanges(3000000, end=12300000))
rangeNoSNPs <- DARegionAnalysis(set, variable_names = "sex", range = range)

##Create multiset with snps
geno <- matrix(rep(c(3, 1, 3, 1, 1, 1), 2), ncol = 6, byrow = T)
colnames(geno) <- c("5723646052_R02C02", "5723646052_R04C01", "5723646052_R05C02",
                    "5723646053_R04C02", "5723646053_R05C02", "5723646053_R06C02")
rownames(geno) <- c("rs3115860", "SNP1-1628854")
map <- AnnotatedDataFrame(data.frame(chromosome = c("chrY", "chr2"), position = c(4241114, 1234321),
                  stringsAsFactors = FALSE))
rownames(map) <- rownames(geno)
snps <- new("SnpSet", call = geno, featureData = map)

multiset <- new("MultiDataSet")
multiset <- add_methy(multiset, set)
multiset <- add_snps(multiset, snps)
rangeSNPsCov <- DARegionAnalysis(multiset, variable_names = "sex", range = range, 
                                    snps_cutoff = 0.05)

eset <- ExpressionSet(matrix(runif(12, max = 15), 3))
annot <- data.frame(chromosome = rep("chrY", 3), start = c(10000, 4500000, 5000000), 
                    end = c(100090, 6000000, 90000000))
fData(eset) <- annot
pData(eset) <- data.frame(sex = c("H", "M", "H", "M"))
esetRes <- DAPipeline(eset, variable_names = "sex", probe_method = "ls")

test_that("Plot CPGs",{
  expect_error(plotFeature(methyOneVar, 12351413654), "feat index must be greater than 0 and smaller than the number of features.")
  expect_error(plotFeature(methyOneVar, -12), "feat index must be greater than 0 and smaller than the number of features.")
  expect_error(plotFeature(methyOneVar, 1, variables = character()), "variables must have one or two values.")
  expect_error(plotFeature(methyOneVar, 1, variables = "character"),"Not all variables are present in set phenodata.")
  expect_error(plotFeature(methyOneVar, 1, variables = c("age", "sex", "status")), "variables must have one or two values.")
  expect_error(plotFeature(methyOneVar, "12"), "feat name is not present in the set.")
  expect_error(plotFeature(methyOneVar, numeric()), "feat must contain only one value.")
  expect_error(plotFeature(methyOneVar, 1:5), "feat must contain only one value.")
  expect_error(plotFeature(MsetEx, 12, "age"), "set must be of class AnalysisResults, MethylationSet or ExpressionSet")
  expect_error(plotBestFeatures(MsetEx, 12, "age"), "set must be of class AnalysisResults, ExpressionSet or MethylationSet")
  expect_warning(plotBestFeatures(rangeNoSNPs, 12421531461, "sex"), "n is greater than the number of betas present in the set.")
})

test_that("Plot Ewas", {
  expect_error(plotEWAS(methyOneVar, variable = character()), "variable must have one value.")
  expect_error(plotEWAS(methyOneVar, variable = c("a", "b")), "variable must have one value.")
  expect_error(plotEWAS(methyOneVar, variable = "a"), "Variable is not present in the model.")
  expect_warning(plotEWAS(rangeNoSNPs), "For a better visualization of the region, use plotRegion.")
  plotEWAS(methyOneVar)
})

test_that("Plot Volcano", {
  plotVolcano(methyOneVar)
  expect_error(plotVolcano(methyOneVar, variable = character()), "variable must have one value.")
  expect_error(plotVolcano(methyOneVar, variable = c("a", "b")), "variable must have one value.")
  expect_error(plotVolcano(methyOneVar, variable = "a"), "Variable is not present in the model.")
})

test_that("QQPlot", {
  expect_error(plotQQ(methyOneVar, variable = character()), "variable must have one value.")
  expect_error(plotQQ(methyOneVar, variable = c("a", "b")), "variable must have one value.")
  expect_error(plotQQ(methyOneVar, variable = "a"), "Variable is not present in phenodata.")
  expect_warning(plotQQ(rangeNoSNPs), "QQplot is not recommended in range analyis.")
})

test_that("plotRDA", {
  expect_error(plotRDA(rangeNoSNPs, n_feat = -1), "n_feat must be greater than 1.")
  expect_warning(plotRDA(rangeNoSNPs, n_feat = 9999), "n_feat is greater than the total number of features in the range.")
})

test_that("Plot Region", {
  expect_error(plotRegion(rangeNoSNPs, variable = character()), "variable must have one value.")
  expect_error(plotRegion(rangeNoSNPs, variable = c("a", "b")), "variable must have one value.")
  expect_error(plotRegion(rangeNoSNPs, variable = "a"), "Variable is not present in modelVariables.\nValid variables are sex")
  expect_error(plotRegion(methyOneVar), "range must be present to use plotRegion with a AnalysisResults")
  expect_error(plotRegion(esetRes), "Results must have a column called position to perform this plot.")
  plotRegion(rangeNoSNPs)
})

test_that("Plot Region R2", {  
  expect_error(plotRegionR2(rangeSNPsCov, 12351413654), "feat index must be greater than 0 and smaller than the number of cpgs.")
  expect_error(plotRegionR2(rangeSNPsCov, -12), "feat index must be greater than 0 and smaller than the number of cpgs.")
  expect_error(plotRegionR2(rangeSNPsCov, "12"), "feat name is not present in the set.")
  expect_error(plotRegionR2(rangeSNPsCov, 1:4), "feat must contain only one value")
  expect_error(plotRegionR2(rangeSNPsCov, 1:5), "feat must contain only one value.")
  plotRDA(rangeSNPsCov)
})

dev.off()
file.remove("Rplots.pdf")