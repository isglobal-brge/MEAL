#' MEAL (Methylation and Expression AnaLizer): Package for analysing methylation and expression data
#' 
#' MEAL has three different categories of important functions: processing, 
#' analysing and plotting.
#' 
#' @section processing: 
#' Functions used to create \code{MEAL} objects and to modify
#' them. Main functions are \link{prepareMethylationSet} and \link{preparePhenotype}
#' 
#' @section analysing:
#' Functions used to perform the analysis of methylation data. \link{DAProbe} 
#' performs per probe analysis and \link{DARegion} performs per region analysis.
#' There are two wrappers: \link{DAPipeline} and \link{DARegionAnalysis} 
#' that performs per probe and per region analysis. The first one analyses the whole
#' methylation sites and the second one only a given region. Finally, 
#' \link{correlationMethExprs} and \link{multiCorrMethExprs} compute the correlation 
#' between methylation and expression probes
#' 
#' @section plotting:
#' Functions used to plot the results of the analysis. Some are interesting for
#' whole methylome analysis (e.g. \link{plotEWAS}) and others for analysis of one 
#' genomic region (e.g. \link{plotRDA})
#' 
#' @docType package
#' @name MEAL
#' 
#' @import BiocGenerics
#' @import Biobase
#' @import ggplot2
#' @import methods
#' @importFrom doParallel registerDoParallel
#' @importFrom DMRcate dmrcate
#' @importFrom GenomicRanges seqnames start end makeGRangesFromDataFrame findOverlaps subjectHits queryHits GRanges
#' @importFrom IRanges IRanges
#' @importFrom limma lmFit contrasts.fit eBayes topTable 
#' @importFrom minfi logit2 ilogit2 bumphunter GenomicRatioSet cpgCollapse blockFinder getBeta
#' @importFrom parallel mclapply
#' @importFrom S4Vectors Rle runValue
#' @importFrom scales alpha
#' @importFrom SNPassoc dominant recessive additive
#' @importFrom snpStats p.value snp.rhs.tests
#' @importFrom sva num.sv sva
#' @importFrom vegan anova.cca RsquareAdj rda
#' @importClassesFrom snpStats SnpMatrix
NULL