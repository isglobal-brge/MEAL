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
#' @export DAPipeline
#' @param set \code{MethylationSet} or \code{ExpressionSet}
#' @param variable_names Character vector with the names of the variables that
#' will be returned as result.
#' @param variable_types Character vector with the types of the variables. As
#' default, variables type won't be changed. 
#' @param covariable_names Character vector with the names of the variables that
#' will be used to adjust the model. 
#' @param covariable_types Character vector with the types of the covariables. As
#' default, variables type won't be changed. 
#' @param equation Character containing the formula to be used to create the model.
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
#' @param num_feat Numeric with the minimum number of cpg beta values to be included in the 
#' results. 
#' @param num_cores Numeric with the number of cores to be used.
#' @param verbose Logical value. If TRUE, it writes out some messages indicating progress. 
#' If FALSE nothing should be printed.
#' @param ... Further arguments passsed to \code{DARegion} function.
#' @return \code{MethylationResult} object 
#' @seealso \code{\link{preparePhenotype}}
#' @examples
#' if (require(minfiData)){
#'  set <- prepareMethylationSet(matrix = getBeta(MsetEx)[1:10, ], pheno = pData(MsetEx))
#'  res <- DAPipeline(set, variable_names = "Sample_Group", probe_method = "ls") 
#'  res
#' }
DAPipeline <-  function(set, variable_names, 
                            variable_types = rep(NA, length(variable_names)), 
                            covariable_names = NULL, 
                            covariable_types = rep(NA, length(covariable_names)), 
                            equation = NULL, num_var = NULL, labels = NULL, sva = FALSE,
                            region_methods = c("bumphunter", "DMRcate"),
                            shrinkVar = FALSE, probe_method = "robust", max_iterations = 100, 
                            num_feat = 50, num_cores = 1, verbose = FALSE,  ...) {
  msg <- validObject(set, test = TRUE)
  
  options("mc.cores" = num_cores)
  
  if (is(msg, "character")){
    message(paste(msg, collapse = "\n"))
    stop("checkProbes and checkSamples might solve validity issues.")
  }
  
  if (ncol(set) == 0 | nrow(set) == 0){
    stop("The set has no beta values")
  }
  
  if (!(is(set, "MethylationSet") | is(set, "ExpressionSet"))){
    stop("set must be a MethylationSet or an ExpressionSet")
  }
  
  if (ncol(pData(set)) == 0 | nrow(pData(set)) == 0){
    stop("Set has no phenotypes")
  }
  
  if (!length(variable_names) | !sum(variable_names %in% colnames(pData(set)))){
    stop("variable_names is empty or is not a valid column of the phenoData of the set.")
  }
  
  if (is(set, "ExpressionSet") & !all(c("chromosome", "start", "end") %in% fvarLabels(set))){
    stop("If set is an ExpressionSet, it must contain a fData with columns chromosome, start and end.")
  }
  
  additives <- FALSE
  allvar <- c(variable_types, covariable_types)
  if (length(allvar) != sum(is.na(allvar))){
    addmask <- allvar == "additive"
    if (sum(addmask, na.rm = TRUE) > 0){
      allvar[addmask] <- "categorical"
      originalpheno <- preparePhenotype(phenotypes = pData(set), 
                                        variable_names = c(variable_names, covariable_names),
                                        variable_types = allvar)
      additives <- TRUE
    }
  }
  
  if (verbose){
    message("Preparing phenotype object")
  }
  pData(set) <- preparePhenotype(phenotypes = pData(set), 
                                 variable_names = c(variable_names, covariable_names),
                                 variable_types = c(variable_types, covariable_types))
  if (verbose){
    message("Creating the model")
  }
  
  model <- createModel(data = pData(set), equation = equation, names = labels)
  
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
    phenovar <- pData(set)[, variable_names, drop = FALSE]
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
  probes <- DAProbe(set = set, model = model, coefficient = coefficient, 
                        shrinkVar = shrinkVar, method = probe_method, 
                        max_iterations = max_iterations)
  
  if (is(set, "MethylationSet")){
    
    if (verbose){
      message("Region Analysis started")
    }
    region <- DARegion(set = set, model = model, methods = region_methods,  
                           coefficient = coefficient, proberes = probes, 
                           verbose = verbose, num_cores = num_cores, ...)
    
  }else{
    region <- list(bumphunter = NA, blockFinder = NA, DMRcate = NA)
  }
  
  if (additives){
    pData(set) <- originalpheno
  }
  results <- analysisResults(set = set, model = model, regionResults = region, 
                             probeResults = probes, num_vars = length(variable_names),
                             num_feat = num_feat)
  return(results)
}

