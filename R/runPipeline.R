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
#' @param weights weights used in the lmFit model (default NULL)
#' @param num_vars Numeric with the number of variables in the matrix for which the
#' analysis will be performed. Compulsory if equation is not null. 
#' @param sva Logical. Should Surrogate Variable Analysis be applied? 
#' (Default: FALSE) 
#' @param betas If \code{set} is a \code{GenomicRatioSet}, should beta values be
#' used? (Default: TRUE)
#' @param range \code{GenomicRanges} with the region used for RDA
#' @param analyses Vector with the names of the analysis to be run (DiffMean and/or DiffVar).
#' @param DiffMean_params List with other parameter passed to \code{runBumphunter} 
#' function.
#' @param DiffVar_params List with other parameter passed to \code{runBumphunter} 
#' function.
#' @param rda_params List with other parameter passed to \code{runRDA} 
#' function.
#' @param method String indicating the method used in the regression: "ls" or 
#' "robust". (Default: "ls")
#' @param big Logical value indicating whether SmartSVA should be used instead of SVA 
#' (TRUE recommended for methylation or when having large number of samples). 
#' Default is FALSE. 
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
                         model = NULL, weights = NULL, num_vars,  
                         sva = FALSE, betas = TRUE, range,
                         analyses = c("DiffMean"),
                         verbose = FALSE, warnings = TRUE, 
                         DiffMean_params = NULL, DiffVar_params = list(coefficient = 1:2),
                         rda_params = NULL, method = "ls",
                         big = FALSE) {
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
    model <- runSVA(mat, model, num_vars = num_vars,
                    big = big)
  }
  
  # Run Differential Mean analysis
  if (verbose){
    message("Probe Analysis started")
  }

  resList <- list()
  
  if ("DiffMean" %in% analyses){
    diffmean <- tryCatch(do.call(runDiffMeanAnalysis, 
                                 c(list(set = mat, model = model, weights = weights,
                                        resultSet = FALSE, 
                                        warnings = warnings, method = method), DiffMean_params)), 
                         error = function(e) e)
    resList$DiffMean <- list(result = diffmean, error = ifelse(is(diffmean, "try-error"), diffmean, NA))
  }
  if ("DiffVar" %in% analyses){
    
  diffvar <- tryCatch(do.call(runDiffVarAnalysis, 
                     c(list(set = mat, model = model, resultSet = FALSE, 
                            warnings = warnings), DiffVar_params)), 
                     error = function(e) e)
  resList$DiffVar <- list(result = diffvar, error = ifelse(is(diffvar, "try-error"), diffvar, NA))
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

runSVA <- function (mat, model, num_vars, big){
  
  model <- data.frame(model)
  model0 <- model[, -c(2:(1+num_vars)), drop = FALSE]
  
  if (big){
   Y.r <- t(resid(lm(t(mat) ~ ., df)))
   n.sv <- isva::EstDimRMT(Y.r, FALSE)$dim + 1
   if (n.sv > 0){
     svobj <- SmartSVA::smartsva.cpp(mat, model, model0, n.sv = n.sv)
     model <- cbind(model, svobj$sv)
   }
   else{
     IQRs <- apply(mat, 1, IQR)
     sv <- sva::sva(mat[IQRs > quantile(IQRs, prob=0.9), ],
                    mod=model, mod0=model0)
     cnames <- c(colnames(model), paste0("SV", 1:sv$n))
     model <- cbind(model, sv$sv)
     colnames(model) <- cnames
   }
  }
  return(model)
}
