#' Plots a rda object from multiCorrMethExprs
#' 
#' Function to plot the rda object obtained in \code{multiCorrMethExprs} function. Samples are displayed with a point
#' whose color can correspond to categorical variable. The 5 methylation probes most related to the first component and the top 
#' 5 most related to the second component will be displayed with a blue arrow. The 5 expression probes most related 
#' to the first component and the top 5 most related to the second component will be labelled. 
#' 
#' @export plotRDAMulti
#' 
#' @param object A rda object like the obtained from multiCorrMethExprs
#' @param color Vector with the colors of the samples. By default, all samples will be displayed in black. 
#' @return A plot in the graphical device
#' @seealso \code{\link{multiCorrMethExprs}}
#' 
#' @details The legend's R-square correspond to the first two RDA components. When a factor is passed to 
#' color argument, the levels will also be displayed in the legend. In addition, an analysis of the variance will be
#' performed to evaluate if the factor can classify the samples. The p-value of this computation is also displayed
#' on the legend. The smaller the p-value, the better the factor can separate the samples. 
plotRDAMulti <- function(object, color = "black"){
  if (!all(class(object) == c("rda", "cca"))){
    stop("Object must be a rda object.")
  }
  
  eigrda <- object$CCA$eig
  eigpca <- object$CA$eig
  r2 <- round(sum(eigrda[1:2])/sum(c(eigrda, eigpca)), 3)
  
  exps <- object$CCA$v
  if (ncol(exps) == 1){
    exp <- cbind(exps, object$CA$v[ , 1])
  }
  e1 <- order(abs(exps[,1]), decreasing=TRUE)
  e2 <- order(abs(exps[,2]), decreasing=TRUE)
  exps.filt <- union(e1[1:5], e2[1:5])
  
  meth <- object$CCA$biplot
  if (ncol(meth) == 1){
    meth <- cbind(meth, object$CA$v[ , 1])
  }
  m1 <- order(abs(meth[,1]), decreasing=TRUE)
  m2 <- order(abs(meth[,2]), decreasing=TRUE)
  m.filt <- union(m1[1:5], m2[1:5])
  
  plot(object, display=c("sp", "cn", "wa"), type = "n")
  points(object, display="sp", select = exps.filt, col = "red", pch=3, cex=0.4)
  text(object, display = "species", select = exps.filt, cex=0.6)
  
  if (length(color) != nrow(object$CCA$u) && color != "black"){
    warning("Length of color is different from the number of samples. All samples will be coloured in black.")
    color <- "black"
  }
  points(object, display = "wa", pch=19, col = color)
  
  text(object, display = "bp", select = m.filt, cex=0.6, col = "blue")
  
  if (is.factor(color)){
    pval <- vegan::adonis(object$CCA$u[, 1:2] ~ color, method = "euclidean")$aov.tab$`Pr(>F)`
    legend("topleft", c(as.expression(bquote(R^2 ~ "=" ~ .(r2))), levels(color), 
                        as.expression(bquote(p.val ~ "=" ~ .(pval)))), 
           col = c("white", 1:nlevels(color)), cex = 0.8, bty = "n", pch = 19)
  } else{
    legend("topleft", as.expression(bquote(R^2 ~ "=" ~ .(r2))), cex = 0.8, bty = "n", pch = 19) 
  }
}
