#' Perform differential methylation analysis
#' 
#' Wrapper for analysing differential methylation and expression at region and probe level.
#' 
#' @details This function is the main wrapper of the package. First, it simplifies the 
#' the set to only contain the common samples between phenotype and features. In addition,
#' it allows to change the class of the variables and to apply genomic models (more 
#' information on \code{preparePhenotype}). Afterwards, analysis per probe and per 
#' region are done merging the results in an \code{AnalysisResults} object.
#' 
#' Default linear model will contain a sum of the variables and covariables. If 
#' interactions are desired, a costum formula can be specified. In that case, variables
#' and covariables must also be specified in order to assure the proper work of the
#' resulting \code{AnalysisResult}. In addition, the number of variables of the model
#' for which the calculation will be done \strong{must} be specified. 
#' 
#' @export runPipeline
#' @param set \code{GenomicRatioSet}, \code{eSet} derived object or 
#' \code{SummarizedExperiment}
#' @param variable_names Character vector with the names of the variables that
#' will be returned as result.
#' @param covariable_names Character vector with the names of the variables that
#' will be used to adjust the model. 
#' @param model Model matrix or formula to get model matrix from \code{set}. 
#' @param num_vars Numeric with the number of variables in the matrix for which the
#' analysis will be performed. Compulsory if equation is not null. 
#' @param sva Logical. Should Surrogate Variable Analysis be applied? 
#' (Default: FALSE) 
#' @param betas If \code{set} is a \code{GenomicRatioSet}, should beta values be
#' used? (Default: TRUE)
#' @param range \code{GenomicRanges} with the region used for RDA
#' @param region_methods Character vector with the methods used in \code{runRegionAnalysis}. If
#' "none", region analysis is not performed.
#' @param DiffMean_params List with other parameter passed to \code{runBumphunter} 
#' function.
#' @param DiffVar_params List with other parameter passed to \code{runBumphunter} 
#' function.
#' @param bumphunter_params List with other parameter passed to \code{runBumphunter} 
#' function.
#' @param blockFinder_params List with other parameter passed to \code{runBlockFinder} 
#' function.
#' @param dmrcate_params List with other parameter passed to \code{runDMRcate} 
#' function.
#' @param rda_params List with other parameter passed to \code{runRDA} 
#' function.
#' @param verbose Logical value. If TRUE, it writes out some messages indicating progress. 
#' If FALSE nothing should be printed.
#' @param warnings Should warnings be displayed? (Default:TRUE) 
#' @return \code{ResultSet} object 
#' @examples
#' if (require(minfiData)){
#' set <- ratioConvert(mapToGenome(MsetEx[1:10,]))
#'  res <- runPipeline(set, variable_names = "Sample_Group") 
#'  res
#' }
runPipeline <-  function(set, variable_names, 
                         covariable_names = NULL, 
                         model = NULL, num_vars,  
                         sva = FALSE, betas = TRUE, range,
                         region_methods = c("bumphunter", "blockFinder", "DMRcate"),
                         verbose = FALSE, warnings = TRUE, 
                         DiffMean_params = NULL, DiffVar_params = list(coefficient = 1:2),
                         bumphunter_params = NULL, 
                         blockFinder_params = NULL, dmrcate_params = NULL, 
                         rda_params = NULL) {
  ### Añadir filtro sondas con NAs
  
  
  if (verbose){
    message("Creating the model")
  }
  
  ## Create linear model
  if (is.null(model)){
    ## Get number of variables of interest
    model <- formula(paste("~", paste(variable_names, collapse = " + ")))
    num_vars <- ncol(createModel(set, model, warnings))

        if ("(Intercept)" %in% colnames(num_vars)){
      num_vars <- num_vars - 1
    }
    
    if (!is.null(covariable_names)){
      model <- formula(paste("~", paste(c(variable_names, covariable_names), 
                                            collapse = " + ")))
    } 
  } 
  model <- createModel(set, model, warnings)
  
  set <- set[, rownames(model)]
  
  
  
  ## Get matrix
  if (is(set, "ExpressionSet")){
    mat <- Biobase::exprs(set)
  } else if (is(set, "GenomicRatioSet")){
    mat <- minfi::getBeta(set)
    if (!betas) {
      mat[mat == 0] <- 1e-3
      mat[mat == 1] <- 1 - 1e-3
      mat <- minfi::logit2(mat)
    }
  } else if (is(set, "SummarizedExperiment")){
    mat <- SummarizedExperiment::assay(set)
  } else {
    stop("set must be an ExpressionSet, MethylationSet, GenomicRatioSet or SummarizedExperiment.")
  }
  
  ## Get Annotation
  fFun <- getFeatureDataFun(set)
  
  if (sva){
    if (verbose){
      message("Computing SVA")
    }
    model <- runSVA(mat, model, num_vars = num_vars)
  }
  
  # Run Differential Mean analysis
  if (verbose){
    message("Probe Analysis started")
  }
  diffmean <- tryCatch(do.call(runDiffMeanAnalysis, 
                      c(list(set = mat, model = model, resultSet = FALSE, 
                             warnings = warnings), DiffMean_params)), 
                      error = function(e) e)
  
  diffvar <- tryCatch(do.call(runDiffVarAnalysis, 
                     c(list(set = mat, model = model, resultSet = FALSE, 
                            warnings = warnings), DiffVar_params)), 
                     error = function(e) e)
  
  resList <- list(DiffMean = list(result = diffmean, 
                                  error = ifelse(is(diffmean, "try-error"), diffmean, NA)),
                  DiffVar = list(result = diffvar, 
                                 error = ifelse(is(diffvar, "try-error"), diffvar, NA)))

  if (class(set) == "GenomicRatioSet") {
    if (verbose){
      message("Region Analysis started")
    }
    region <- runRegionAnalysis(set = set, model = model, methods = region_methods,  
                             coefficient = 2, verbose = verbose, resultSet = FALSE, 
                             bumphunter_params = bumphunter_params, 
                             blockFinder_params = blockFinder_params, 
                             dmrcate_params = dmrcate_params)
    
    resList <- c(resList, region)
  }
  if (!missing(range)){
    rda <- do.call(runRDA, c(list(set = set, model = model, num_vars = num_vars, 
                                  resultSet = FALSE, range = range), rda_params))
    resList$RDA <- list(result = rda, error = NA)
  } 
  
  res <- create_resultset(fOrigin = "runPipeline", 
                          lResults = resList,
                          fData = list(main = fFun(set)),
                          lOptions = list(sva = sva))
  
  return(res)
}

runSVA <- function (mat, model, num_vars){
  
  df <- data.frame(model)
  model0 <- model[, -c(2:(1+num_vars)), drop = FALSE]
  
  Y.r <- t(resid(lm(t(mat) ~ ., df)))
  n.sv <- isva::EstDimRMT(Y.r, FALSE)$dim + 1
  if (n.sv > 0){
    svobj <- SmartSVA::smartsva.cpp(mat, model, model0, n.sv = n.sv)
    model <- cbind(model, svobj$sv)
  }
  model
}
