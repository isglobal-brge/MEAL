#' Plot RDA results
#' 
#' @name plotRDA
#' @aliases plotRDA
#' @export
#' 
#' @param object \code{ResultSet}
#' @param pheno data.frame with the variables used to color the samples. 
#' @param n_feat Numeric with the number of cpgs to be highlighted. Default: 5.
#' @param main Character with the plot title.
#' @param alpha Numeric with the alpha level for colour transparance. Default: 1; no
#' transparency.
#' @return A plot is generated on the current graphics device.
#' @examples
#' if (require(minfiData)){
#' set <- ratioConvert(mapToGenome(MsetEx[1:10,]))
#' model <- model.matrix(~set$sex)
#' rda <- runRDA(set, model)
#' plotRDA(rda, pheno = data.frame(factor(set$sex)))
#' }
plotRDA <- function(object, pheno = data.frame(), n_feat = 5, main = "RDA plot", 
                    alpha = 1){

  stopifnot("RDA" %in% names(object))
  ans <- getAssociation(object, rid = "RDA")
  
  if (is.null(dim(pheno))){
    pheno <- data.frame(pheno)
  }
  
  classes <- vapply(pheno, class, character(1))
  pheno[, classes == "character"] <- lapply(pheno[, classes == "character", drop = FALSE], function(x) factor(x))
  
  classes <- vapply(pheno, class, character(1))
  factormatrix <- pheno[, classes == "factor", drop = FALSE]
  phenocont <- pheno[, classes %in% c("numeric", "integer"), drop = FALSE]
  factor <- as.vector(sapply(colnames(factormatrix), function(x) levels(factormatrix[, x])))
  
  
  
  if (ncol(factormatrix) > 1){
    factormatrix[] <- lapply(factormatrix, as.character)
    
    phenofactor <- as.factor(sapply(rownames(factormatrix),
                               function(x) paste0(unlist(factormatrix[x, ]), collapse = "")))
  } else if (ncol(factormatrix) == 1) {
    phenofactor <- factormatrix[, 1]
  } else {
    phenofactor <- 1
  }
  
  r2 <- ans$rdaR2
  pval <- ans$pval
  if (pval < 1e-3){
    pval <- format(signif(pval, 3), scientific = TRUE)
  } else{
    pval <- round(pval, 3)
  }
  
  ans$CCA$centroids <- NULL
  if (length(factor)){
    ans$CCA$centroids <- t(data.matrix(data.frame(
      lapply(colnames(factormatrix), function(x) 
        getCentroids(ans, factormatrix[ , x])),
      check.names = FALSE)))
  }
  
  ## Remove from biplot variables already present in centroids
  ans$CCA$biplot[] <- 0
  
  if (ncol(phenocont)) {
    xx <- data.matrix(phenocont)
    ans$CCA$biplot <- (1/sqrt(colSums(xx^2))) * crossprod(xx, ans$CCA$u)
  }

  temp <- ans$CCA$v
  
  if (ncol(temp) == 1){
    temp <- cbind(temp, ans$CA$v[ , 1])
  }
  o1 <- order(abs(temp[,1]), decreasing=TRUE)
  o2 <- order(abs(temp[,2]), decreasing=TRUE)
  
  plot(ans, display=c("sp", "cn", "wa"), type="n", main = main, scaling = 3)
  
  filter <- union(o1[1:n_feat], o2[1:n_feat])
    
  text(ans, display = "species", select = filter, cex=0.6, scaling = 3)
  
  ## Set points
  pch = (15:18)[as.numeric(phenofactor)]
  points(ans, display = "wa", col = ggplot2::alpha(as.numeric(phenofactor), alpha), 
         pch = pch, scaling = 3)
  
  if (length(factor)){
    vegan::ordilabel(ans, display = "cn", col = "blue", label = factor)
  }
  
  if (ncol(phenocont)){
    text(ans, display = "bp", col = "blue", scaling = 3)
  }
  
  legend("topleft", c(as.expression(bquote(R^2 ~ "=" ~ .(round(r2, 3)))), paste("p-value:", pval), levels(phenofactor)), 
         col = c("white", "white", ggplot2::alpha(1:nlevels(phenofactor), alpha)), 
         cex = 0.8, bty = "n", pch = c(1, 1, c(15:18)[seq_len(nlevels(phenofactor))]))
}


getCentroids <- function(rda, factor){
  posRDA <- data.frame(rda$CCA$wa)
  splitpos <- split(posRDA, factor)
  if (is.null(dim(splitpos[[1]]))){
    splitpos <- lapply(splitpos, function(RDA1) data.frame(RDA1))
  }
  pos <- data.frame(lapply(splitpos, colMeans))
  colnames(pos) <- names(splitpos)
  colnames(pos) <- paste0(colnames(factor), colnames(pos))
  return(pos)
}