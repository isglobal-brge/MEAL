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
#' @param vars_meth Character vector with the names of the variables that will be
#' used to obtain the methylation residuals. By default, none is used and residuals 
#' are not computed.
#' @param vars_meth_types Character vector with the types of the methylation variables. 
#' By default, variables type won't be changed. 
#' @param vars_exprs Character vector with the names of the variables that will
#' be used to obtain the expression residuals. By default, none is used and
#' residuals are not computed.
#' @param vars_exprs_types Character vector with the types of the expression variables. 
#' By default, variables type won't be changed.
#' @param flank Numeric with the number of pair bases used to define the cpg-expression 
#' probe pairs.
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
correlationMethExprs <- function(multiset, vars_meth = NULL, vars_exprs = NULL, 
                                 vars_meth_types = rep(NA, length(vars_meth)), 
                                 vars_exprs_types = rep(NA, length(vars_exprs)), 
                                 flank = 250000, num_cores = 1, verbose = TRUE){
  
  options("mc.cores" = num_cores)
  
  if (!is(multiset, "MultiDataSet")){
    stop("multiset must be a MultiDataSet")
  }
  
  if (!all(c("methylation", "expression") %in% names(multiset))){
    stop("multiset must contain methylation and expression data.")
  }
  
  mset <- multiset[["methylation"]]
  eset <- multiset[["expression"]]
  
  if (ncol(mset) == 0 | nrow(mset) == 0){
    stop("The mset has no beta values")
  }
  
  if (ncol(eset) == 0 | nrow(eset) == 0){
    stop("The eset has no expression values")
  }
  
  if (!all(c("chromosome", "start", "end") %in% fvarLabels(eset))){
    stop("eset must contain a featureData with columns chromosome, start and end")
  }
  
  if (!is(flank, "numeric") || length(flank) > 1 || flank < 0){
    stop("flank must be a positive integer")
  }
  
  if (verbose){
    message("Calculating cpg-expression probe pairs")
  }
  
  pairs <- pairsExprsMeth(fData(mset), fData(eset), flank)
  
  if (nrow(pairs) == 0){
    warning("There are no expression probes in the range of the cpgs. An empty data.frame will be returned.")
    return(data.frame(cpg = character(0), exprs = character(0), Beta = integer(0), 
                      se = integer(0), P.Value = integer(0), adj.P.val = integer(0)))
  }
  
  mset <- mset[unique(pairs[ , 1]), ]
  eset <- eset[unique(pairs[ , 2]), ]
  
  
  if (verbose){
    message("Computing residuals")
  }
  
  methres <- setResidues(mset, vars_names = vars_meth, vars_types = vars_meth_types)
  exprsres <- setResidues(eset, vars_names = vars_exprs, vars_types = vars_exprs_types)
    
  if (verbose){
    message("Computing correlation Methylation-Expression")
  }
  
  regvals <- mclapply(1:nrow(pairs), 
                    function(x) residualsCorr(methres[pairs[x, 1], ], exprsres[pairs[x, 2], ]))
  
  res <- data.frame(pairs, t(data.frame(regvals)))
  colnames(res) <- c("cpg", "exprs", "Beta", "se", "P.Value")
  res$adj.P.Val <- p.adjust(res$P.Value, "BH")
  res <- res[order(res$adj.P.Val), ]
  rownames(res) <- NULL
  res
}


pairsExprsMeth <- function(meth, exprs, flank){
  meth$start <- meth$position - flank
  meth[ meth < 0] <- 0
  meth$end <- meth$position + flank
  methGR <- GenomicRanges::makeGRangesFromDataFrame(meth, seqnames.field = "chromosome", 
                                                    ignore.strand = TRUE)
  exprsGR <- GenomicRanges::makeGRangesFromDataFrame(exprs, seqnames.field = "chromosome",
                                                     ignore.strand = TRUE)
  pairs <- GenomicRanges::findOverlaps(exprsGR, methGR,  type = "within")
  pairs <- data.frame(cpg = rownames(meth)[GenomicRanges::subjectHits(pairs)], 
                      exprs = rownames(exprs)[GenomicRanges::queryHits(pairs)], 
                      stringsAsFactors = FALSE)
  pairs
}

residualsCorr <- function (methy_res, exprs_res){
  fit <- lm(exprs_res ~ methy_res)
  return(c(summary(fit)$coef[2, 1], summary(fit)$coef[2,2], summary(fit)$coef[2,4]))
}

setResidues <- function(set, vars_names, vars_types){
  if (length(vars_names) != 0){
    if (ncol(pData(set)) == 0 | nrow(pData(set)) == 0){
      warning("set has no phenotypes. No residues will be computed")
    } else if (!sum(vars_names %in% colnames(pData(set)))){
      if (is(set, "ExpressionSet")){
        warning("vars_exprs is/are not a valid column of the eset phenoData. No residues will be computed")
      } else{
        warning("vars_meth is/are not a valid column of the mset phenoData. No residues will be computed")
      }
    }else{
      pData(set) <- preparePhenotype(phenotypes = pData(set), 
                                     variable_names = vars_names,
                                     variable_types = vars_types)
      model <- createModel(data = pData(set))
      if (is(set, "ExpressionSet")){
        vals <- exprs(set)
      } else{
        vals <- betas(set)
        vals <- minfi::logit2(vals)
      }
      res <- residuals(limma::lmFit(vals, model), vals)
      if (is(set, "MethylationSet")){
        res <- minfi::ilogit2(res)
      }
      return(res)
    }
  }
  if (is(set, "ExpressionSet")){
    return(exprs(set))
  } else{
    return(betas(set))
  }
}