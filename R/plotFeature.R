#' Plot values of a feature
#' 
#' Plot values of a feature splitted by one or two variables.
#' 
#' @export
#' @param set \code{AnalysisResults}, \code{AnalysisRegionResults}, \code{ExpressionSet} or 
#' \code{MethylationSet}
#' @param feat Numeric with the index of the feature or character with its name.
#' @param variables Character vector with the names of the variables to be used 
#' in the splitting. Two variables is the maximum allowed. Note: default values 
#' are only valid for \code{MethylationResults} objects. 
#' @return A plot is generated on the current graphics device.
#' @examples
#' if (require(minfiData)){
#'  set <- prepareMethylationSet(getBeta(MsetEx)[1:1000, ], 
#'  pheno = pData(MsetEx))
#'  plotFeature(set, 1, variables = "Sample_Group")
#'  }

plotFeature <- function(set, feat, variables = variableNames(set)[1]){
  if (length(feat) != 1){
    stop("feat must contain only one value.")
  }
  if (is.numeric(feat)) {
    if (is(set, "AnalysisResults")){
      if (feat < 1 || feat > length(feats(set))){
        stop("feat index must be greater than 0 and smaller than the number of cpgs.")
      }
      feat <- feats(set)[feat]
    }else if (is(set, "MethylationSet") || is(set, "ExpressionSet")){
      if (feat < 1 || feat > nrow(set)){
        stop("feat index must be greater than 0 and smaller than the number of cpgs.")
      }
      feat<- rownames(set)[feat]
    }else{
      stop("set must be of class AnalysisResults, MethylationSet or ExpressionSet")
    }
  }
  if (is(set, "AnalysisResults")) {
    if (!feat %in% feats(set)){
      stop("feat name is not present in the set.")
    }
    values <- featvals(set)[feat, ]
  } else if (is(set, "MethylationSet")){
    if (!feat %in% rownames(betas(set))){
      stop("feat name is not present in the set.")
    }
    values <- betas(set)[feat, ]
  } else if (is(set, "ExpressionSet")){
    if (!feat %in% rownames(exprs(set))){
      stop("feat name is not present in the set.")
    }
    values <- exprs(set)[feat, ]
  } else{
    stop("set must be of class MethylationResults, MethylationSet or ExpressionSet")
  }
  
  if (length(variables) < 1 || length(variables) > 2){
    stop("variables must have one or two values.")
  }  
  if (!all(variables %in% colnames(pData(set)))){
    stop("Not all variables are present in set phenodata.")
  }
  
  values <- as.vector(values)
  
  datamatrix <- data.frame(vals = values, factor = pData(set)[ , variables[1]])
  if (length(variables) == 2){
    datamatrix <- cbind(datamatrix, factor2 = pData(set)[ , variables[2]] )
  }
  if (is.numeric(datamatrix$factor)){
    p <- ggplot2::ggplot(datamatrix, ggplot2::aes_string(x = "factor", y = "vals")) + ggplot2::geom_point() +
      ggplot2::geom_smooth(method = lm, se = FALSE) + ggplot2::scale_x_continuous(name = "") + 
      ggplot2::ggtitle(feat)
    if (length(variables) == 2){
      p <- p +  ggplot2::facet_grid(. ~ factor2) + ggplot2::scale_x_continuous(name = variables[2])
    }
  }else{  
    groups <- split(datamatrix, datamatrix$factor)
    means <- data.frame(means = sapply(groups, function(x) mean(x[ , 1])), factor = names(groups))
    p <- ggplot2::ggplot(datamatrix, ggplot2::aes_string(x = "factor", y = "vals")) + 
      ggplot2::geom_point(alpha = 1/2, position = ggplot2::position_jitter(width = 0.1)) + 
      ggplot2::geom_errorbar(data = means, ggplot2::aes_string(x = "factor", y = "means", ymax = "means", ymin = "means"),
                             width = 0.4, size = 1) +
      ggplot2::scale_x_discrete(name = "") + ggplot2::ggtitle(feat)
    if (length(variables) == 2){
      p <- p +  ggplot2::facet_grid(. ~ factor2) + ggplot2::scale_x_discrete(name = variables[2])
    }
  }
  if (max(datamatrix$vals) > 1){
    p <- p + ggplot2::scale_y_continuous("Expression")
  }else{
    p <- p + ggplot2::scale_y_continuous("Methylation (Beta)", limits = c(0, 1))
  }
    
  print(p)
  
}