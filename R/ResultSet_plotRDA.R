#' Plot RDA results
#' 
#' @name plotRDA
#' @aliases plotRDA
#' @export
#' 
#' @param object \code{ResultSet}
#' @param pheno data.frame with the variables used to color the samples. 
#' @param n_feat Numeric with the number of cpgs to be highlighted.
#' @param main Character with the plot title.
#' @return A plot is generated on the current graphics device.
#' @examples
#' if (require(minfiData)){
#' set <- ratioConvert(mapToGenome(MsetEx[1:10,]))
#' model <- model.matrix(~set$sex)
#' rda <- runRDA(set, model)
#' plotRDA(rda, pheno = data.frame(factor(set$sex)))
#' }
plotRDA <- function(object, pheno, n_feat = 5, main = "RDA plot"){

  ## TO DO
  ## Check passing a pheno with character
  ## Check pass a vector in pheno
  
  stopifnot("RDA" %in% names(object))
  ans <- getAssociation(object, rid = "RDA")
  
  classes <- sapply(colnames(pheno), function(x) class(pheno[, x]))
  factormatrix <- pheno[, classes == "factor", drop = FALSE]
  phenocont <- pheno[, classes == "numeric", drop = FALSE]
  factor <- as.vector(sapply(colnames(factormatrix), function(x) levels(factormatrix[, x])))
  
  if (ncol(factormatrix) > 1){
    phenofactor <- as.factor(sapply(rownames(factormatrix),
                               function(x) paste0(unlist(factormatrix[x, ]), collapse = "")))
  }else{
    phenofactor <- factormatrix[, 1]
  }
  
  r2 <- ans$rdaR2
  pval <- ans$pval
  if (pval < 1e-3){
    pval <- format(signif(pval, 3), scientific = TRUE)
  } else{
    pval <- round(pval, 3)
  }
  
  if (!is.null(factor)){
    ans$CCA$centroids <- t(data.matrix(data.frame(
      lapply(colnames(factormatrix), function(x) 
        getCentroids(ans, factormatrix[ , x])),
      check.names = FALSE)))
  }
  
  ## Remove from biplot variables already present in centroids
  pattern <- paste(colnames(factormatrix), collapse = "|")
  ans$CCA$biplot[grep(pattern, rownames(ans$CCA$biplot)), ] <- 0
  
  temp <- ans$CCA$v
  
  if (ncol(temp) == 1){
    temp <- cbind(temp, ans$CA$v[ , 1])
  }
  o1 <- order(abs(temp[,1]), decreasing=TRUE)
  o2 <- order(abs(temp[,2]), decreasing=TRUE)
  
  plot(ans, display=c("sp", "cn", "wa"), type="n", main = main, scaling = 3)
  
  filter <- union(o1[1:n_feat], o2[1:n_feat])
    
  text(ans, display = "species", select = filter, cex=0.6, scaling = 3)
  
  points(ans, display = "wa", col = ggplot2::alpha(as.numeric(phenofactor), 0.5), 
         pch=19, scaling = 3)
  
  if (!is.null(factor)){
    text(ans, display = "cn", col = "blue", label = factor, scaling = 3)
  }
  
  legend("topleft", c(as.expression(bquote(R^2 ~ "=" ~ .(round(r2, 3)))), paste("p-value:", pval), levels(phenofactor)), 
         col = c("white", "white", 1:nlevels(phenofactor)), cex = 0.8, bty = "n", pch = 19) 
}


getCentroids <- function(rda, factor){
  posRDA <- data.frame(rda$CCA$u)
  splitpos <- split(posRDA, factor)
  if (is.null(dim(splitpos[[1]]))){
    splitpos <- lapply(splitpos, function(RDA1) data.frame(RDA1))
  }
  pos <- data.frame(lapply(splitpos, colMeans))
  colnames(pos) <- names(splitpos)
  colnames(pos) <- paste0(colnames(factor), colnames(pos))
  return(pos)
}