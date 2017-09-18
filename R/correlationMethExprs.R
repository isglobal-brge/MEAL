#' Computes the correlation between methylation and expression
#' 
#' Estimates the correlation between methylation and expression. When there are known 
#' variables that affect methylation and/or expression, their effect can be substracted
#' using a linear model and then the residuals are used. 
#' 
#' For each cpg, a range is defined by the position of the cpg plus the flank parameter
#' (upstream and downstream). Only those expression probes that are entirely in 
#' this range will be selected. For these reason, it is required that the \code{ExpressionSet} 
#' contains a featureData with the chromosome and the starting and ending positions 
#' of the probes.  
#' 
#' @export correlationMethExprs
#' 
#' @param multiset \code{MultiDataSet} containing a \code{methylation} and an 
#' \code{expression} slots.
#' @param meth_set_name Character vector with the name of the \code{MultiDataSet}'s slot containing methylation
#' data.
#' @param exprs_set_name Character vector with the name of the \code{MultiDataSet}'s slot containing expression
#' data.
#' @param vars_meth Character vector with the names of the variables that will be
#' used to obtain the methylation residuals. By default, none is used and residuals 
#' are not computed.
#' @param vars_exprs Character vector with the names of the variables that will
#' be used to obtain the expression residuals. By default, none is used and
#' residuals are not computed.
#' @param sel_cpgs Character vector with the name of the CpGs used in the analysis. If empty, all the CpGs of the 
#' methylation set will be used. 
#' @param flank Numeric with the number of pair bases used to define the cpg-expression 
#' probe pairs.
#' @param betas If \code{set} is a \code{GenomicRatioSet}, should beta values be
#' used? (Default: TRUE)
#' @param num_cores Numeric with the number of cores to be used.
#' @param verbose Logical value. If TRUE, it writes out some messages indicating progress. 
#' If FALSE nothing should be printed.
#' @return Data.frame with the results of the linear regression:
#' \itemize{
#'  \item cpg: Name of the cpg
#'  \item exprs: Name of the expression probe
#'  \item beta: coefficient of the methylation change
#'  \item se: standard error of the beta 
#'  \item P.Value: p-value of the beta coefficient
#'  \item adj.P.Val: q-value computed using B&H
#' }
correlationMethExprs <- function(multiset, 
                                 meth_set_name = NULL, exprs_set_name = NULL,
                                 vars_meth = NULL, vars_exprs = NULL, 
                                 sel_cpgs, flank = 250000, betas = TRUE,
                                 num_cores = 1, verbose = TRUE){
  
  # To DO
  ## Añadir check a vars_meth y vars_exprs
  
  ######################################################################################
  ## Data Checking
  # Check object is a MultiDataSet
  if (!is(multiset, "MultiDataSet")){
    stop("multiset must be a MultiDataSet")
  }
  
  # Check meth_set_name and exprs_set_name are characters
  if (!(is.character(meth_set_name) | is.null(meth_set_name))){
    stop("meth_set_name must be a character.")
  } 
  if (!(is.character(exprs_set_name) | is.null(exprs_set_name))){
    stop("exprs_set_name must be a character.")
  } 
  
  # Add the dataset type to the names
  meth_set_name <- paste(c("methylation", meth_set_name), collapse = "+")
  exprs_set_name <- paste(c("expression", exprs_set_name), collapse = "+")
  
  
  # Check our object has the right sets
  if (!all(c(meth_set_name, exprs_set_name) %in% names(multiset))){
    stop("multiset must contain meth_set_name and exprs_set_name.")
  }
  
  if (!all(c(meth_set_name, exprs_set_name) %in% rowRangesElements(multiset))){
    stop("multiset must contain meth_set_name and exprs_set_name.")
  }
  
  # Check that flank is numeric, positive and is only one number
  if (!is(flank, "numeric") || length(flank) > 1 || flank < 0){
    stop("flank must be a positive integer")
  }
  #####################################################################################
  
  #####################################################################################
  ## Preparation of data
  #############################
  ## Preparation of data 1
  # Select only our sets
  multiset <- multiset[, c(meth_set_name, exprs_set_name)]
  
  ## Select common samples
  multiset <- commonSamples(multiset)
  #############################
  ## Preparation of data 2
  
  mset <- multiset[[meth_set_name]]
  eset <- multiset[[exprs_set_name]]
  
  #############################
  ## Preparation of data 3
  if (!missing(sel_cpgs)){
    if (!is.character(sel_cpgs) | length(sel_cpgs) == 0){
      stop("sel_cpgs must be a character vector with the name of the CpGs to be used.")
    } else{
      mset <- mset[sel_cpgs, ]
    }
  }
  
  #############################
  ## Preparation of data 4
  # Compute Methylation-Expression pairs
  rangesMeth <- rowRanges(multiset)[[meth_set_name]]
  rangesMeth <- rangesMeth[featureNames(mset)]
  rangesExprs <- rowRanges(multiset)[[exprs_set_name]]
  #####################################################################################
  
  #####################################################################################
  ## Implementation of the algorithm
  #############################
  ## Implementation of the algorithm 1
  start(rangesMeth) <- start(rangesMeth) - flank
  end(rangesMeth) <- end(rangesMeth) + flank
  pairs <- GenomicRanges::findOverlaps(rangesExprs, rangesMeth,  type = "within")
  pairs <- data.frame(cpg = rownames(mset)[S4Vectors::subjectHits(pairs)], 
                      exprs = rownames(eset)[S4Vectors::queryHits(pairs)], 
                      stringsAsFactors = FALSE)
  if (nrow(pairs) == 0){
    warning("There are no expression probes in the range of the cpgs. An empty data.frame will be returned.")
    return(data.frame(cpg = character(0), exprs = character(0), Beta = integer(0), 
                      se = integer(0), P.Value = integer(0), adj.P.val = integer(0)))
  }
  #############################
  ## Filter sets to only features in the pairs
  mset <- mset[unique(pairs[ , 1]), ]
  eset <- eset[unique(pairs[ , 2]), ]
  
  
  if (verbose){
    message("Computing residuals")
  }
  
  #############################
  ## Implementation of the algorithm 2
  # Computing residuals
  methres <- setResidues(mset, variable_names = vars_meth, betas = betas)
  exprsres <- setResidues(eset, variable_names = vars_exprs)
  #############################
  if (verbose){
    message("Computing correlation Methylation-Expression")
  }
  
  #############################
  ## Implementation of the algorithm 3
  residualsCorr <- function (methy_res, exprs_res){
    fit <- lm(exprs_res ~ methy_res)
    return(c(summary(fit)$coef[2, 1], summary(fit)$coef[2,2], summary(fit)$coef[2,4]))
  }
  regvals <- mclapply(1:nrow(pairs), 
                      function(x) residualsCorr(methres[pairs[x, 1], ], exprsres[pairs[x, 2], ]), 
                      mc.cores = num_cores)
  #####################################################################################
  
  #####################################################################################
  ## Formatting the results
  res <- data.frame(pairs, t(data.frame(regvals)))
  colnames(res) <- c("cpg", "exprs", "Beta", "se", "P.Value")
  res$adj.P.Val <- p.adjust(res$P.Value, "BH")
  res <- res[order(res$adj.P.Val), ]
  rownames(res) <- NULL
  res
  #####################################################################################
}

setResidues <- function(set, variable_names, betas = TRUE){
  
  if (is(set, "ExpressionSet")){
    res <- Biobase::exprs(set)
  } else if (is(set, "GenomicRatioSet")){
    res <- minfi::getBeta(set)
    if (!betas) {
      res[res == 0] <- 1e-3
      res[res == 1] <- 1 - 1e-3
      res <- minfi::logit2(res)
    }
  } else if (is(set, "SummarizedExperiment")){
    res <- Biobase::assays(set)
  } 
  
  if (!is.null(variable_names)){
    model <- formula(paste("~", paste(variable_names, collapse = " + ")))
    ## Get number of variables of interest
    model <- createModel(set, model)
    
    res <- residuals(limma::lmFit(res, model), res)
  }
 
  return(res)
}