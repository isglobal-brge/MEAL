#' Method to add a slot of expression to \code{MultiDataSet}.
#'
#' This method adds or overwrites the slot \code{"expression"} of an
#' \code{MultiDataSet} with the content of the given \code{ExpressionSet}.
#' 
#' @aliases add.genexp
#' @rdname add.genexp-methods
#' @exportMethod add.genexp
#' @param object \code{MultiDataSet} that will be filled.
#' @param gexpSet \code{ExpressionSet} to be used to fill the slot.
#' @param warnings Logical to indicate if warnings will be displayed.
#' @return A new \code{MultiDataSet} with the slot \code{"expression"}
#' filled.
#' @examples 
#' multi <- new("MultiDataSet")
#' eset <- new("ExpressionSet", exprs = matrix(runif(4), 2))
#' fData(eset) <- data.frame(chromosome = c("chr1", "chr2"), start = c(12414, 1234321),
#'  end = c(121241, 12412412414), stringsAsFactors = FALSE)
#' multi <- add.genexp(multi, eset)
setGeneric("add.genexp", function(object, gexpSet, warnings = TRUE) {
  standardGeneric("add.genexp")
})

#' Method to add a slot of methylation to \code{MultiDataSet}.
#'
#' This method adds or overwrites the slot \code{"methylation"} of an
#' \code{MultiDataSet} with the content of the given \code{MethylationSet}.
#'
#' @rdname add.methy-methods
#' @aliases add.methy-methods
#' @param object \code{MultiDataSet} that will be filled.
#' @param methySet \code{MethylationSet} to be used to fill the slot.
#' @param warnings Logical to indicate if warnings will be displayed.
#' @return A new \code{MultiDataSet} with the slot \code{"methylation"}
#' filled.
#' @exportMethod add.methy
#' @examples 
#' if (require(MEALData)){
#'  multi <- new("MultiDataSet")
#'  betavals <- betavals[1:100, ]  ## To speed up the example, the beta values are reduced
#'  methy <- prepareMethylationSet(betavals, pheno)
#'  multi <- add.methy(multi, methy)
#' }
setGeneric("add.methy", function(object, methySet, warnings = TRUE) {
  standardGeneric("add.methy")
})

#' Method to add a slot to \code{MultiDataSet}.
#'
#' This method adds or overwrites a slot of a \code{MultiDataSet} with the content 
#' of the given \code{eSet}.
#'
#' @rdname add.set-methods
#' @aliases add.set-methods
#' @param object \code{MultiDataSet} that will be filled.
#' @param set Object derived from \code{eSet} to be used to fill the slot.
#' @param dataset.name Character with the name of the slot to be filled.
#' @param warnings Logical to indicate if warnings will be displayed.
#' @return A new \code{MultiDataSet} with a slot filled.
#' @exportMethod add.set
#' @examples 
#' multi <- new("MultiDataSet")
#' eset <- new("ExpressionSet", exprs = matrix(runif(10), 5))
#' multi <- add.set(multi, eset, "exampledata")
#' @examples 
#' multi <- new("MultiDataSet")
#' snps <- new("SnpSet", call = matrix(runif(10), 5))
#' multi <- add.snps(multi, snps)
setGeneric("add.set", function(object, set, dataset.name, warnings = TRUE) {
  standardGeneric("add.set")
})

#' Method to add a slot of SNPs to \code{MultiDataSet}.
#'
#' This method adds or overwrites the slot \code{"snps"} of an
#' \code{MultiDataSet} with the content of the given \code{SnpSet}.
#'
#' @rdname add.snps-methods
#' @aliases add.snps
#' @param object \code{MultiDataSet} that will be filled.
#' @param snpSet \code{SnpSet} to be used to fill the slot.
#' @param warnings Logical to indicate if warnings will be displayed.
#' @return A new \code{MultiDataSet} with the slot \code{"snps"}
#' filled.
#' @exportMethod add.snps
setGeneric("add.snps", function(object, snpSet, warnings = TRUE) {
  standardGeneric("add.snps")
})

#' @export 
setGeneric("betas", function(object){
  standardGeneric("betas")
})

#' @export
setGeneric("blocks", function(object){
  standardGeneric("blocks")
})

#' @export
setGeneric("bumps", function(object){
  standardGeneric("bumps")
})

#' Filter \code{MethylationSet} probes
#' 
#' This function selects probes present in the annotation matrix. Probes without
#' annotation and annotation values without beta values are discarded. 
#' 
#' @rdname checkProbes-methods
#' @export 
#' @aliases checkProbes
#' @param object \code{MethylationSet}
#' @return \code{MethylationSet} containing the common samples.
#' @examples 
#' if (require(MEALData)){
#'  betavals <- betavals[1:100, ]  ## To speed up the example, the beta values are reduced
#'  methy <- prepareMethylationSet(betavals, pheno)
#'  checkProbes(methy)
#' }
setGeneric("checkProbes", function(object){
  standardGeneric("checkProbes")
})

#' Modify a \code{MethylationSet} to only contain common samples
#' 
#' This function removes samples that have beta values but no phenotypes and vice versa.
#' If snps object is present, only samples present in the three set are retained.
#' 
#' @name checkSamples
#' @rdname checkSamples-methods
#' @aliases checkSamples 
#' @export 
#' 
#' @param object \code{MethylationSet}
#' @return \code{MethylationSet} containing the common samples.
#' @examples 
#' if (require(MEALData)){
#'  betavals <- betavals[1:100, ]  ## To speed up the example, the beta values are reduced
#'  methy <- prepareMethylationSet(betavals, pheno)
#'  checkSamples(methy)
#' }
setGeneric("checkSamples", function(object){
  standardGeneric("checkSamples")
})

#' @export
setGeneric("covariableNames", function(object){
  standardGeneric("covariableNames")
})

#' @export
setGeneric("dmrCate", function(object){
  standardGeneric("dmrCate")
})

#' Exports results data.frames to csv files.
#' 
#' Exports results to csv files. If more than one variable is present, subfolders
#' with the name of the variable are created. For each variable, four files will 
#' be generated: probeResults.csv, dmrCateResults.csv, bumphunterResults.csv
#' and blockFinderResults.csv
#' 
#' @name exportResults
#' @rdname exportResults-methods
#' @aliases exportResults 
#' @export 
#' 
#' @param object \code{MethylationResults} or \code{MethylationRegionResults} 
#' @param dir Character with the path to export.
#' @param prefix Character with a prefix to be added to all file names.
#' @param vars Character vector with the names of the variables to be exported. Note: 
#' names should be that of the model. 
#' @return Files are saved into the given folder.
#' @examples
#' if (require(minfiData)){
#' set <- prepareMethylationSet(getBeta(MsetEx)[1:10,], pheno = pData(MsetEx))
#' methyOneVar <- DAPipeline(set, variable_names = "sex", probe_method = "ls")
#' exportResults(methyOneVar)
#' }
setGeneric("exportResults", function(object, dir = "./", prefix = NULL, 
                                     vars = modelVariables(object)){
  standardGeneric("exportResults")
})

#' @export 
setGeneric("feats", function(object){
  standardGeneric("feats")
})

#' @export 
setGeneric("featvals", function(object){
  standardGeneric("featvals")
})

#' Get all probes related to gene
#' 
#' Given a \code{MethylationResults} and a gene name returns the results of the 
#' analysis of all the probes of the gene. 
#' 
#' @name getGeneVals
#' @rdname getGeneVals-methods
#' @export
#' @param object \code{MethylationResults}
#' @param gene Character with the name of the gene
#' @return List of data.frames with the results of the analysis of the probes belonging to 
#' the gene 
#' @examples
#' if (require(minfiData)){
#' set <- prepareMethylationSet(getBeta(MsetEx)[1:10,], pheno = pData(MsetEx))
#' methyOneVar <- DAPipeline(set, variable_names = "sex", probe_method = "ls")
#' getGeneVals(methyOneVar, "TSPY4")
#' }
setGeneric("getGeneVals", function(object, gene){
  standardGeneric("getGeneVals")
})


#' Transforms beta values to M-values
#' 
#' Given a \code{MethylationSet} or a \code{AnalysisResults} returns 
#' the matrix of M values using a logit2 transformation. Betas equal to 0 will 
#' be transformed to threshold and betas equal to 1, to 1 - threshold. 
#' 
#' @name getMs
#' @rdname getMs-methods
#' @export
#' @param object \code{MethylationSet} or \code{AnalysisResults}
#' @param threshold Numeric with the threshold to avoid 0s and 1s. 
#' @return Matrix with the M values.
#' @examples
#' if (require(minfiData)){
#' set <- prepareMethylationSet(MsetEx[1:100, ], pData(MsetEx))
#' mvalues <- getMs(set)
#' head(mvalues)
#' }
setGeneric("getMs", function(object, threshold = 0.0001){
  standardGeneric("getMs")
})

#' @export 
setGeneric("getRange", function(object){
  standardGeneric("getRange")
})

#' @export 
setGeneric("getRDA", function(object){
  standardGeneric("getRDA")
})

#' @export
setGeneric("model", function(object){
  standardGeneric("model")
})

#' @export
setGeneric("modelVariables", function(object){
  standardGeneric("modelVariables")
})

#' Plot a Manhattan plot with the probe results
#' 
#' Plot log p-value for each chromosome positions. Highlighting cpgs inside a range
#' is allowed.
#' 
#' @name plotEWAS
#' @rdname plotEWAS-methods
#' @aliases plotEWAS
#' @export
#' 
#' @param object \code{AnalysisResults} or \code{AnalysisRegionResults} 
#' @param variable Character with the variable name used to obtain the probe results.
#' Note: model name should be used. Original variable name might not be valid.  
#' @param range GenomicRange whose cpgs will be highlighted
#' @return A plot is generated on the current graphics device.
#' @examples
#' if (require(minfiData)){
#' betas <- getBeta(MsetEx)[floor(seq(1, nrow(MsetEx), 10000)), ]
#' set <- prepareMethylationSet(betas, pheno = pData(MsetEx))
#' methyOneVar <- DAPipeline(set, variable_names = "sex", probe_method = "ls")
#' plotEWAS(methyOneVar)
#' }
setGeneric("plotEWAS", function(object, variable = modelVariables(object)[[1]], range = NULL){
  standardGeneric("plotEWAS")
})

#' Plot RDA results
#' 
#' @name plotRDA
#' @rdname plotRDA-methods
#' @aliases plotRDA
#' @export
#' 
#' @param object \code{AnalysisRegionResults}
#' @param n_feat Numeric with the number of cpgs to be highlighted.
#' @return A plot is generated on the current graphics device.
#' @examples
#' if (require(minfiData) & require(GenomicRanges)){
#' set <- prepareMethylationSet(getBeta(MsetEx), pheno = pData(MsetEx))
#' range <- GenomicRanges::GRanges(seqnames=Rle("chrY"), 
#' ranges = IRanges(3000000, end=12300000))
#' rangeNoSNPs <- DARegionAnalysis(set, variable_names = "sex", range = range)
#' plotRDA(rangeNoSNPs)
#' }
setGeneric("plotRDA", function(object, n_feat = 5){
  standardGeneric("plotRDA")
})

#' Plot of the region
#' 
#' Plot of the beta values againts their position. Data is taken from probe analysis.
#' Cpgs with a p-value smaller than 0.05 (without adjusting) are blue and points
#' with a p-value greater than 0.05 are red.
#' 
#' @name plotRegion
#' @rdname plotRegion-methods
#' @aliases plotRegion
#' @export
#' 
#' @param object \code{AnalysisResults} or \code{AnalysisRegionResults} 
#' @param variable Character with the variable name used to obtain the probe results.
#' Note: model name should be used. Original variable name might not be valid.  
#' @param range GenomicRange whose cpgs will be shown (only for \code{AnalysisResults}
#' objects)
#' @return A plot is generated on the current graphics device.
#' @examples
#' if (require(minfiData) & require(GenomicRanges)){
#' set <- prepareMethylationSet(getBeta(MsetEx), pheno = pData(MsetEx))
#' range <- GenomicRanges::GRanges(seqnames=Rle("chrY"), 
#' ranges = IRanges(3000000, end=12300000))
#' rangeNoSNPs <- DARegionAnalysis(set, variable_names = "sex", range = range)
#' plotRegion(rangeNoSNPs)
#' }
setGeneric("plotRegion", function(object, variable = modelVariables(object)[[1]], range = NULL){
  standardGeneric("plotRegion")
})

#' Plot R2 region values
#' 
#' @name plotRegionR2
#' @rdname plotRegionR2-methods
#' @aliases plotRegionR2
#' @export
#' 
#' @param object \code{MethylationRegionResults}
#' @param feat Numeric with the index of the  feature or character with its name.
#' @param ... Further arguments passed to plotLM
#' @return A plot is generated on the current graphics device.
setGeneric("plotRegionR2", function(object, feat, ...){
  standardGeneric("plotRegionR2")
})

#' QQ plot of probe analysis
#' 
#' Generate a QQ plot using probe results.
#' @name plotQQ
#' @rdname plotQQ-methods
#' @aliases plotQQ 
#' @export
#' 
#' @param object \code{AnalysisResults} or \code{AnalysisRegionResults} 
#' @param variable Character with the variable name used to obtain the probe results.
#' Note: model name should be used. Original variable name might not be valid.
#' @return A plot is generated on the current graphics device.
#' @examples
#' if (require(minfiData)){
#' betas <- getBeta(MsetEx)[floor(seq(1, nrow(MsetEx), 10000)), ]
#' set <- prepareMethylationSet(betas, pheno = pData(MsetEx))
#' methyOneVar <- DAPipeline(set, variable_names = "sex", probe_method = "ls")
#' plotQQ(methyOneVar)
#' }
setGeneric("plotQQ", function(object, variable = modelVariables(object)[[1]]){
  standardGeneric("plotQQ")
})

#' Make a Volcano plot with the probe results
#' 
#' Plot log p-value versus the change in expression/methylation.
#' 
#' @name plotVolcano
#' @rdname plotVolcano-methods
#' @aliases plotVolcano
#' @export
#' 
#' @param object \code{MethylationResults} or \code{MethylationRegionResults} 
#' @param variable Character with the variable name used to obtain the probe results.
#' Note: model name should be used. Original variable name might not be valid.  
#' @return A plot is generated on the current graphics device.
#' @examples
#' if (require(minfiData)){
#' betas <- getBeta(MsetEx)[floor(seq(1, nrow(MsetEx), 10000)), ]
#' set <- prepareMethylationSet(betas, pheno = pData(MsetEx))
#' methyOneVar <- DAPipeline(set, variable_names = "sex", probe_method = "ls")
#' plotEWAS(methyOneVar)
#' }
setGeneric("plotVolcano", function(object, variable = modelVariables(object)[[1]]){
  standardGeneric("plotVolcano")
})

#' @export
setGeneric("probeResults", function(object){
  standardGeneric("probeResults")
})

#' @export 
setGeneric("regionLM", function(object){
  standardGeneric("regionLM")
})

#' @export 
setGeneric("regionPval", function(object){
  standardGeneric("regionPval")
})

#' @export 
setGeneric("regionR2", function(object){
  standardGeneric("regionR2")
})

#' @export
setGeneric("regionResults", function(object){
  standardGeneric("regionResults")
})

#' @export 
setGeneric("snps", function(object) {
  standardGeneric("snps")
})

#' @export 
setGeneric("snpsPvals", function(object){
  standardGeneric("snpsPvals")
})

#' @export 
setGeneric("snpsVar", function(object){
  standardGeneric("snpsVar")
})

#' Get the top features associated with the RDA
#' 
#' Get a list of the features significantly associated to the first two RDA components
#' 
#' @name topRDAhits
#' @rdname topRDAhits-methods
#' @aliases topRDAhits
#' @export
#' 
#' @param object \code{AnalysisRegionResults} 
#' @param pval numeric with the p-value threshold. Only features with a p-values below this threshold will be shown. 
#' @return data.frame with the features, the component, the correlation and the p-value
#' @examples
#' if (require(minfiData) & require(GenomicRanges)){
#' set <- prepareMethylationSet(getBeta(MsetEx), pheno = pData(MsetEx))
#' range <- GenomicRanges::GRanges(seqnames=Rle("chrY"), 
#' ranges = IRanges(3000000, end=12300000))
#' rangeNoSNPs <- DARegionAnalysis(set, variable_names = "sex", range = range)
#' topRDAhits(rangeNoSNPs)
#' }
setGeneric("topRDAhits", function(object, pval = 0.05){
  standardGeneric("topRDAhits")
})


#' @export
setGeneric("variableNames", function(object){
  standardGeneric("variableNames")
})