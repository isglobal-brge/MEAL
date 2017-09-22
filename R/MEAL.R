#' MEAL (Methylation and Expression AnaLizer): Package for analysing methylation 
#' and expression data
#' 
#' MEAL is a package designed to facilitate the analysis methylation and expression 
#' data. The package can analyze one dataset and can find correlations between
#' methylation and expression data. MEAL has a vignette that explains the main
#' functionalities of the package. 
#' 
#' @docType package
#' @name MEAL
#' 
#' @import BiocGenerics
#' @import Biobase
#' @import methods
#' @import MultiDataSet
#' @importFrom doParallel registerDoParallel
#' @importFrom DMRcate cpg.annotate dmrcate
#' @importFrom GenomicRanges seqnames start end makeGRangesFromDataFrame findOverlaps GRanges
#' @importFrom ggplot2 aes aes_string alpha facet_grid geom_bar geom_errorbar geom_hline geom_line geom_point geom_polygon geom_smooth geom_text geom_vline ggplot ggtitle position_jitter scale_colour_manual scale_fill_manual scale_x_continuous scale_y_continuous theme
#' @importFrom graphics legend points text
#' @importFrom IRanges IRanges
#' @importFrom limma lmFit contrasts.fit eBayes topTable 
#' @importFrom matrixStats rowVars
#' @importFrom minfi logit2 ilogit2 bumphunter GenomicRatioSet cpgCollapse blockFinder getBeta
#' @importFrom parallel mclapply
#' @importFrom permute how
#' @importFrom S4Vectors Rle runValue subjectHits queryHits
#' @importFrom SNPassoc dominant recessive additive
#' @importFrom snpStats p.value snp.rhs.tests
#' @importFrom stats contrasts cor.test formula lm model.matrix p.adjust prcomp qbeta qt
#' @importFrom sva num.sv sva
#' @importFrom utils write.csv2
#' @importFrom vegan anova.cca RsquareAdj rda
#' @importClassesFrom snpStats SnpMatrix
NULL