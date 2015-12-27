#' @describeIn AnalysisRegionResults  Plot RDA results
#' @aliases AnalysisRegionResults-methods
#' @param n_feat Numeric with the number of features to be highlighted.
setMethod(
 f = "plotRDA",
 signature = "AnalysisRegionResults",
 definition =  function(object, n_feat = 5){
  if (n_feat > length(feats(object))){
    warning("n_feat is greater than the total number of cpgs in the range.")
    n_feat <- length(feats(object))
  }
  if (n_feat < 1){
    stop("n_feat must be greater than 1.")
  }
  
  phenos <- pData(object)
  phenos <- phenos[ , variableNames(object), drop = FALSE]
  classes <- sapply(colnames(phenos), function(x) class(phenos[, x]))
  factormatrix <- phenos[, classes == "factor", drop = FALSE]
  phenocont <- phenos[, classes == "numeric", drop = FALSE]
  factor <- as.vector(sapply(colnames(factormatrix), function(x) levels(factormatrix[, x])))
  if (ncol(factormatrix) > 1){
    phenofactor <- as.factor(sapply(rownames(factormatrix),
                               function(x) paste0(unlist(factormatrix[x, ]), collapse = "")))
  }else{
    phenofactor <- unlist(factormatrix)
  }
  
  ans <- getRDA(object)
  r2 <- round(vegan::RsquareAdj(ans)$r.squared,3)
  
  if (!is.null(factor)){
    ans$CCA$centroids <- t(data.matrix(data.frame(
      lapply(colnames(factormatrix), function(x) getCentroids(ans, factormatrix[ , x, drop = FALSE])),
      check.names = FALSE)))
  }
  
  
  temp <- ans$CCA$v
  
  if (ncol(temp) == 1){
    temp <- cbind(temp, ans$CA$v[ , 1])
  }
  o1 <- order(abs(temp[,1]), decreasing=TRUE)
  o2 <- order(abs(temp[,2]), decreasing=TRUE)
  
  plot(ans, display=c("sp", "cn", "wa"), type="n")
  
  filter <- union(o1[1:n_feat], o2[1:n_feat])
    
  points(ans, display="sp", select = !(1:nrow(temp) %in% filter), col="red",
         pch=3, cex=0.4)
  text(ans, display = "species", select = filter, cex=0.6)
  
  points(ans, display = "wa", col = ggplot2::alpha(as.numeric(phenofactor), 0.5), pch=19)
  
  if (!is.null(factor)){
    text(ans, display = "cn", col = "blue", label = factor)
  }
  
  legend("topleft", c(as.expression(bquote(R^2 ~ "=" ~ .(r2))), levels(phenofactor)), 
         col = c("white", 1:nlevels(phenofactor)), cex = 0.8, bty = "n", pch = 19) 
}
)

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