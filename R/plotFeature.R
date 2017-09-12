#' Plot values of a feature
#' 
#' Plot values of a feature splitted by one or two variables.
#' 
#' @export
#' @param set \code{ExpressionSet}, \code{GenomicRatioSet} or \code{SummarizedExperiment.}
#' @param feat Numeric with the index of the feature or character with its name.
#' @param variables Character vector with the names of the variables to be used 
#' in the splitting. Two variables is the maximum allowed.
#' @param betas If \code{set} is a \code{GenomicRatioSet}, should beta values be
#' used? (Default: TRUE)
#' @return A plot is generated on the current graphics device.
#' @examples
#' if (require(minfiData)){
#'  set <- prepareMethylationSet(getBeta(MsetEx)[1:1000, ], 
#'  pheno = data.frame(pData(MsetEx)))
#'  plotFeature(set, 1, variables = "Sample_Group")
#'  }

plotFeature <- function(set, feat, variables = colnames(pheno)[1], betas = TRUE){
  if (length(feat) != 1){
    stop("feat must contain only one value.")
  }
  methplot <- FALSE
  ## Get matrix
  if (is(set, "ExpressionSet")){
    vals <- Biobase::exprs(set)[feat, ]
  } else if (is(set, "GenomicRatioSet")){
    if (betas) {
      vals <- minfi::getBeta(set)[feat, ]
      methplot <- TRUE
    } else {
      vals <- minfi::getM(set)[feat, ]
    }
  } else if (is(set, "SummarizedExperiment")){
    vals <- SummarizedExperiment::assay(set)[feat, ]
  } else {
    stop("set must be an ExpressionSet, GenomicRatioSet or SummarizedExperiment.")
  }
  
  ## Change getters depending on set class
  if (is(set, "eSet")){
    pheno <- Biobase::pData(set)
  }
  else if (is(set, "SummarizedExperiment")){
    pheno <- SummarizedExperiment::colData(set)
  } 
  
  if (length(variables) < 1 || length(variables) > 2){
    stop("variables must have one or two values.")
  }  
  if (!all(variables %in% colnames(pheno))){
    stop("Not all variables are present in set phenodata.")
  }
  
  values <- as.vector(vals)
  
  datamatrix <- data.frame(Values = values, factor = pheno[ , variables[1]])
  if (length(variables) == 2){
    datamatrix <- cbind(datamatrix, factor2 = pheno[ , variables[2]] )
  }
  if (is.numeric(datamatrix$factor)){
    p <- ggplot2::ggplot(datamatrix, ggplot2::aes_string(x = "factor", y = "Values")) + 
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = lm, se = FALSE) + 
      ggplot2::scale_x_continuous(name = "") + 
      ggplot2::ggtitle(feat)
    if (length(variables) == 2){
      p <- p +  ggplot2::facet_grid(. ~ factor2) + 
        ggplot2::scale_x_continuous(name = variables[2])
    }
  }else{  
    groups <- split(datamatrix, datamatrix$factor)
    means <- data.frame(means = sapply(groups, function(x) mean(x[ , 1])), 
                        factor = names(groups))
    p <- ggplot2::ggplot(datamatrix, 
                         ggplot2::aes_string(x = "factor", y = "Values", fill = "factor")) + 
      ggplot2::geom_violin(alpha = 1/2) +
      ggplot2::geom_point(alpha = 1/2, position = ggplot2::position_jitter(width = 0.1)) + 
      ggplot2::geom_errorbar(data = means, 
                             ggplot2::aes_string(x = "factor", y = "means", ymax = "means", 
                                                 ymin = "means"),
                             width = 0.4, size = 1) +
      ggplot2::scale_x_discrete(name = "") + ggplot2::ggtitle(feat) + 
      ggplot2::scale_fill_discrete(name = variables[1])
    if (length(variables) == 2){
      p <- p +  ggplot2::facet_grid(. ~ factor2) + 
        ggplot2::scale_x_discrete(name = variables[2])
    }
  }
  if (methplot){
    p <- p + ggplot2::scale_y_continuous("Methylation (Beta)", limits = c(0, 1))
  }
    
 p
}