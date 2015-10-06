#' MethylationSet instances
#' 
#' Container with the data needed to perform methylation analysis. \code{MethylationSet}
#' inherits from \code{eSet} and contains \code{meth} matrix as assay data member.
#' 
#' FeatureData, which contains annotation data, is required to perform any of the
#' analysis.
#' 
#' @export
#' @rdname MethylationSet-class
#' @name MethylationSet
#' @aliases MethylationSet-class MethylationSet-methods
#' @slot assayData Contains matrices with equal dimensions, and with column number 
#' equal to nrow(phenoData). assayData must contain a matrix meth with rows representing 
#' features (e.g., methylation probes sets) and columns representing samples. 
#' @slot phenoData See \linkS4class{eSet}
#' @slot annotation See \linkS4class{eSet}
#' @slot featureData See \linkS4class{eSet}. fData should contain at least chromosome
#' and positions columns. 
setClass (
  Class = "MethylationSet",
  contains = "eSet", 
  prototype = prototype(new("VersionedBiobase",
                            versions = c(classVersion("eSet"), MethylationSet = "1.0.0")))
)

#' Class \code{MultiDataSet}
#'
#' The class \code{MultiDataSet} is a superior class to store multiple
#' datasets in form of triplets (assayData-phenoData-featureData). The
#' restriction is that the samples of the multiple datasets must be the same.
#'
#' The names of the three lists (\code{assayData}, \code{phenoData} and
#' \code{featureData})must be the same.
#'
#' @name MultiDataSet-class
#' @rdname MultiDataSet-class
#' @exportClass MultiDataSet
#' @slot assayData List of \code{assayData} elements.
#' @slot phenoData List of \code{AnnotatedDataFrame} containing the phenoData
#' of each \code{assayData}.
#' @slot featureData List of \code{AnnotatedDataFrame} containing the featureData
#' of each \code{assayData}.
#' @slot return_method List of functions used to create the original \code{eSet} objects.
setClass(
  Class = "MultiDataSet",
  representation = representation(
    assayData = "list",
    phenoData = "list",
    featureData = "list",
    return_method = "list"
  )
)

#' AnalysisResults instances
#' 
#' Container with the results of per probe and per region analyses.
#' 
#' @export
#' @rdname AnalysisResults-class
#' @name AnalysisResults
#' @aliases AnalysisResults-class AnalysisResults-methods
#' 
#' @slot originalclass Character with the class of the object used to perform the 
#' analysis
#' @slot features Matrix with the values of the most significant features.
#' @slot phenotypes AnnotatedDataFrame with the phenotypes.
#' @slot model Matrix with the model used in the analysis
#' @slot sampleNames Character vector with the names of the samples
#' @slot variableNames Character vector with the names of the variables used in the
#' analysis. Names are equal to those find in phenotypes.
#' @slot covariableNames Character vector with the names of the covariables used in the
#' analysis. Names are equal to those find in phenotypes.
#' @slot results List of data.frames with the results of per probe analysis. Names are
#' those of the model.
#' @slot DMRcate List of data.frames with the results of \code{DMRcate}. Names are
#' those of the model.
#' @slot Bumphunter List of data.frames with the results of \code{Bumphunter}. Names are
#' those of the model.
#' @slot BlockFinder List of data.frames with the results of \code{BlockFinder}. Names are
#' those of the model.
setClass (
  Class = "AnalysisResults",
  representation(
    originalclass = "character", 
    features = "matrix",
    phenotypes = "AnnotatedDataFrame",
    model = "matrix",
    sampleNames = "character",
    variableNames = "character",
    covariableNames = "character",
    results = "list",
    DMRcate = "list",
    Bumphunter = "list",
    BlockFinder = "list"
  )
)

#' AnalysisRegionResults instances
#' 
#' \code{AnalysisResults} heir with the analyses performed in a range of the 
#' whole genome.
#'  
#' @export
#' @rdname AnalysisRegionResults-class
#' @name AnalysisRegionResults
#' @aliases AnalysisRegionResults-class AnalysisRegionResults-methods
#'  
#' @slot range \code{GenomicRanges} used to perform the analysis.
#' @slot snps Character vector with the snps that are correlated to at least one cpg.
#' @slot snpsPvals Data.frame with the results of the correlation test SNP-cpg.
#' @slot snpsVar Numeric with the variability of the SNP matrix explained by the components 
#' used to adjust the linear model. 
#' @slot rda \code{rda} object from \code{vegan} package with the results of RDA
#' analysis in the range.
#' @slot regionLM List with the R2 of the linear model of beta values against our 
#' variable of interest and against significant SNPs for each cpg.
#' @slot regionR2 Numeric with the R2 of the region calculated using a redundancy analysis.
#' @slot regionPval Numeric with the pval of the region's R2. 
setClass (
  Class = "AnalysisRegionResults",
  contains = "AnalysisResults",
  representation(
    range = className("GRanges","GenomicRanges"),
    snps = "character",
    snpsPvals = "data.frame",
    snpsVar = "numeric",
    rda = "list",
    regionLM = "list",
    regionR2 = "numeric",
    regionPval = "numeric"   
  )
)