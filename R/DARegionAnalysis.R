#' Analyse methylation or expression in a specific range
#' 
#' Methylation analysis in a genomic range.
#' 
#' @details Set is filtered to the range specified. Probe analysis and DMR detection
#' are run using the filtering set. Finally, RDA test of the region is performed, 
#' returning the R2 between the variables and the beta matrix and a p-value of 
#' this R2.   
#' 
#' @export DARegionAnalysis
#' @param set \code{MethylationSet}, \code{ExpressionSet} or \code{MultiDataSet}.
#' @param omicset In a \code{MultiDataSet} allows to choose between methylation and
#' expression (valid values are: "methylation" or "expression").
#' @param range \code{GenomicRanges} with the desired range.
#' @param variable_names Character vector with the names of the variables that
#' will be returned as result.
#' @param variable_types Character vector with the types of the variables. By
#' default, variables type won't be changed. 
#' @param covariable_names Character vector with the names of the variables that
#' will be used to adjust the model. 
#' @param covariable_types Character vector with the types of the covariables. By
#' default, variables type won't be changed. 
#' @param equation String containing the formula to be used to create 
#' the model.
#' @param num_var Numeric with the number of variables in the matrix for which the
#' analysis will be performed. Compulsory if equation is not null. 
#' @param labels Character vector with the labels of the variables.
#' @param sva Logical indicating if Surrogate Variable Analysis should be applied.
#' @param region_methods Character vector with the methods used in \code{DARegion}. If
#' "none", region analysis is not performed.
#' @param shrinkVar Logical indicating if shrinkage of variance should be applied 
#' in probe analysis.
#' @param probe_method Character with the type of linear regression applied in 
#' probe analysis ("ls" or "robust")
#' @param max_iterations Numeric with the maximum of iterations in the robust
#' regression.
#' @param num_cores Numeric with the number of cores to be used.
#' @param verbose Logical value. If TRUE, it writes out some messages indicating progress. 
#' If FALSE nothing should be printed.
#' @param nperm Numeric with the number of permutations used to compute RDA p-values.
#' @param ... Further arguments passsed to \code{DAPipeline} function.
#' @return \code{AnalysisRegionResult} object 
#' @seealso \code{\link{preparePhenotype}}, \code{\link{DAPipeline}}
#' @examples
#' if (require(minfiData)){
#'  set <- prepareMethylationSet(getBeta(MsetEx)[1:1000, ], 
#'  pheno = data.frame(pData(MsetEx)))
#'  range <- GenomicRanges::GRanges(seqnames=Rle("chrX"), 
#'  ranges = IRanges(30000, end = 123000000))
#'  res <- DARegionAnalysis(set, range = range, variable_names = "Sample_Group",
#'  probe_method = "ls") 
#'  res
#' }
DARegionAnalysis <- function(set, range, omicset = "methylation", variable_names, 
                             variable_types = rep(NA, length(variable_names)), 
                             covariable_names = NULL, 
                             covariable_types = rep(NA, length(covariable_names)), 
                             equation = NULL,   num_var = NULL,  labels = NULL, 
                             sva = FALSE,  
                             region_methods = c("blockFinder", "bumphunter", "DMRcate"),
                             shrinkVar = FALSE, probe_method = "robust", max_iterations = 100, 
                             num_cores = 1, verbose = FALSE, nperm = 1e3,  ...){

    if (!(is(set, "MethylationSet") | is(set, "ExpressionSet"))){
        stop("set must be a MethylationSet or an ExpressionSet.")
    }
    
    msg <- validObject(set, test = TRUE)
    
    
    if (is(msg, "character")){
        message(paste(msg, collapse = "\n"))
        stop("checkProbes and checkSamples might solve validity issues.")
    }
    if (ncol(set) == 0 | nrow(set) == 0){
        stop("The set has no beta values")
    }
    
    if (ncol(pData(set)) == 0 | nrow(pData(set)) == 0){
        stop("Set has no phenotypes")
    }
    
    if (is(set, "ExpressionSet") & !all(c("chromosome", "start", "end") %in% fvarLabels(set))){
        stop("If set is an ExpressionSet, it must contain a fData with columns chromosome, start and end.")
    }
    
    if (!length(variable_names) | !sum(variable_names %in% colnames(pData(set)))){
        stop("variable_names is empty or is not a valid column of the phenoData of the set.")
    }
    
    if (!is(range, "GenomicRanges")){
        stop("range must be a GenomicRanges object.")
    }
    
    if (length(range) == 0){
        stop("range is empty.")
    }
    
    if (length(range) != 1){
        stop("range must be a GenomicRanges with only one range.")
    }
    
    if (verbose){
        message("Filtering the set.")
    }
    
    filt.set <- filterSet(set, range)
    if (!nrow(filt.set)){
        stop("There are no features in the range specified.")
    }  
    
    
    options("mc.cores" = num_cores)
    
    additives <- FALSE
    allvar <- c(variable_types, covariable_types)
    if (length(allvar) != sum(is.na(allvar))){
        addmask <- allvar == "additive"
        if (sum(addmask, na.rm = TRUE) > 0){
            allvar[addmask] <- "categorical"
            originalpheno <- preparePhenotype(phenotypes = pData(filt.set), 
                                              variable_names = c(variable_names, covariable_names),
                                              variable_types = allvar)
            additives <- TRUE
        }
    }
    
    if (verbose){
        message("Preparing phenotype object")
    }
    pData(filt.set) <- preparePhenotype(phenotypes = pData(filt.set), 
                                   variable_names = c(variable_names, covariable_names),
                                   variable_types = c(variable_types, covariable_types))
    if (verbose){
        message("Creating the model")
    }
    
    model <- createModel(data = pData(filt.set), equation = equation, names = labels)
    
    if (sva){
        if (verbose){
            message("Computing SVA")
        }
        if (is(set, "ExpressionSet")){
            vals <- exprs(set)
        }else{
            vals <- betas(set)
        }
        if (nrow(vals) > 1){
            if (length(covariable_names) == 0){
                model0 <- model.matrix(~ 1, pData(set))
            }else{
                model0 <- model.matrix(~ ., pData(set)[ , covariable_names, drop = FALSE])
            }
            n.sv <- sva::num.sv(vals, model)
            if (n.sv > 0){
                svobj <- sva::sva(vals, model, model0, n.sv = n.sv)
                model <- cbind(model, svobj$sv)
            }
        }
    }
    
    if (!is.null(equation) && is.null(num_var)){
        stop("If costum equation, num_var must be specified.")
    }
    if (is.null(num_var)){
        phenovar <- pData(filt.set)[, variable_names, drop = FALSE]
        continuous_vars <- sapply(variable_names, function(x) is.numeric(phenovar[,x]))
        if (sum(continuous_vars) < ncol(phenovar)){
            num_var <- sum(sapply(variable_names[!continuous_vars], 
                                  function(x) nlevels(phenovar[, x]) - 1), continuous_vars)
        }else{
            num_var <- sum(continuous_vars)
        }
    }
    coefficient <- 2:(1 + num_var)
    
    if (verbose){
        message("Probe Analysis started")
    }
    probes <- DAProbe(set = filt.set, model = model, coefficient = coefficient, 
                      shrinkVar = shrinkVar, method = probe_method, 
                      max_iterations = max_iterations)
    
    if (is(set, "MethylationSet")){
        
        if (verbose){
            message("Region Analysis started")
        }
        region <- DARegion(set = filt.set, model = model, methods = region_methods,  
                           coefficient = coefficient, verbose = verbose, 
                           num_cores = num_cores, ...)
        
    }else{
        region <- list(bumphunter = NA, blockFinder = NA, DMRcate = NA)
    }
    
    if (additives){
        pData(filt.set) <- originalpheno
    }
    
    if (verbose){
        message("Compute RDA")
    }
    varsmodel <- model[, 1:(1 + num_var)]
    
    if (ncol(varsmodel) != ncol(model)){
        covarsmodel <- model[, -c(1:(1 + num_var))]        
    }else{
        covarsmodel <- NULL
    }
    
    RDAres <- RDAset(filt.set, varsmodel, covarsmodel)
    rdapval <- vegan::anova.cca(RDAres, permutations = permute::how(nperm = nperm))[["Pr(>F)"]][1]
    rdaR2 <- vegan::RsquareAdj(RDAres)$r.squared
    
    ### Compute Global RDA vals
    if (is(set, "MethylationSet")){
        fullMat <- MultiDataSet::betas(set)
    }else{
        fullMat <- exprs(set)
    }
    
    globalRDA <- computeRDAR2(fullMat = fullMat, varsmodel = varsmodel, 
                              covarsmodel = covarsmodel, featNum = nrow(filt.set),
                              R2 = rdaR2, nperm = nperm)
    
    
    RDAres <- list(rda = RDAres, rdapval = rdapval, rdaR2 = rdaR2, 
                   globalpval = globalRDA["globalpval"], 
                   globalR2 = globalRDA["globalR2"])
    
    output1 <-  analysisResults(set = filt.set, model = model, regionResults = region, 
                                           probeResults = probes, 
                                num_vars = length(variable_names),
                                           num_feat = nrow(filt.set))
    output <- analysisRegionResults(analysisResults = output1, range = range,
                                        rdaRes = RDAres)                           
return(output)
}