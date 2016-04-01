#' Perform per probe analysis
#' 
#' Compute statistics (t estimate and p-value) for methylation or expression data
#' using linear or robust linear regression.
#' 
#' @export DAProbe
#' 
#' @param set \code{MethylationSet}, matrix of M-values or \code{ExpressionSet}. 
#' @param model Matrix with the model
#' @param coefficient Numeric with the index of the model matrix used to perform
#'  the analysis. If a vector is supplied, a list will be returned.
#' @param shrinkVar Logical indicating if shrinkange of variance should be 
#' performed.
#' @param method String indicating the method used in the regression ("ls" or 
#' "robust")
#' @param max_iterations Numeric indicating the maximum number of iterations
#' done in the robust method. 
#' @return Data.frame or list of data.frames containing intercept and slope values. 
#' If the set is a MethylationSet, probe's position, chromosome and the nearest
#'  gene are also returned.
#' @examples
#' if (require(minfiData)){
#'  mvalues <- getM(MsetEx)[1:100, ]
#'  model <- model.matrix(~ Sample_Group, data = pData(MsetEx)) 
#'  res <- DAProbe(mvalues, model, method = "ls")
#'  head(res)
#' }
DAProbe <- function(set, model, coefficient = 2, shrinkVar = FALSE, method = "robust",
                        max_iterations = 100)
{
  if (ncol(set) == 0 | nrow(set) == 0){
    stop("The set is empty.")
  }
  if (ncol(model) == 0 | nrow(model) == 0){
    stop("The model matrix is empty.")
  }
  if (ncol(set) != nrow(model)){
    stop("The number of samples is different in the set and in the model")
  }
  if (length(method) != 1){
    stop("method should be one of \"ls\", \"robust\"")
  }
  if (!method %in% c("ls", "robust")){
    stop("method should be one of \"ls\", \"robust\"")
  }
  if (length(max_iterations) == 0 || !is(max_iterations, "numeric") || max_iterations < 1){
    warning("max_iterations must be a numeric greater than 1. The default value will be used.")
    max_iterations <- 100
  }
  if (length(coefficient) > 1){
    results <- lapply(coefficient,
                      function(x) DAProbe(set = set, model = model,
                                              coefficient = x, shrinkVar = shrinkVar, 
                                              method = method, 
                                              max_iterations = max_iterations))
    names(results) <- colnames(model)[coefficient]
    return(results)                  
  }
  
  if (is(set, "MethylationSet"))
  {
    msg <- validObject(set, test = TRUE)
    if (is(msg, "character")){
      message(paste(msg, collapse = "\n"))
      stop("checkProbes and checkSamples might solve validity issues.")
    }
    M <- getMs(set)
  } 
  else if (is(set, "matrix")){
    M <- set
  } else if (is(set, "ExpressionSet")){
    M <- exprs(set)
  } else{
    stop("set must be a MethylationSet or a matrix.")
  }
  
  fit <- limma::lmFit(M, model, method = method, maxit = max_iterations)
  if (shrinkVar) {
    fit <- limma::contrasts.fit(fit, contrasts(model))
  }
  fit <- limma::eBayes(fit)
  toptable <- limma::topTable(fit, coef = coefficient, sort.by = "none", 
                              number = nrow(fit$coefficients), confint = TRUE) 
  CI <- qt(0.975, fit$df.total[1])

  if (is(set, "ExpressionSet")){
    tab <- data.frame(intercept = fit$coefficients[ , 1],
                      fit$coefficients[ , 2], sd = (toptable[ , 1] - toptable[ , 2])/CI, 
                      toptable[ , 5:7])
  }
  else {
    tab <- data.frame(intercept = minfi::ilogit2(fit$coefficients[ , 1]),
                      minfi::ilogit2(fit$coefficients[ , 2] + fit$coefficients[ , 1]) -
                        minfi::ilogit2(fit$coefficients[ , 1]), sd = (toptable[ , 1] - toptable[ , 2])/CI, 
                      toptable[ , 5:7])
  }
  rownames(tab) <- rownames(M)
  
  o <- order(tab$P.Value)
  tab <- tab[o, , drop = FALSE]
  colnames(tab)[2] <- colnames(model)[coefficient]
  if (is(set, "MethylationSet")){
    positions <- fData(set)[ , c("chromosome", "position", "genes", "group")]
  } else if (is(set, "ExpressionSet")){
    positions <- fData(set)
  }
  if (is(set, "MethylationSet") | is(set, "ExpressionSet")){
    tab <- cbind(tab, positions[rownames(tab), , drop = FALSE])
  }
  tab <- tab[order(tab$P.Value, na.last = TRUE),]
  tab
}