#' Plot a vector of R2
#' 
#' Plot a vector of R2 where the first value is the main variable and the last one,
#' if named \emph{covariates} is treated as covariates.
#' 
#' @export plotLM
#' 
#' @param Rsquares Numerical vector of R2
#' @param title Character with the plot title
#' @param feat_name Name of the feature used in default title.
#' @param variable_name Character for the first column name
#' @param max_columns Numerical with the maximum number of columns to be plotted.
plotLM <- function(Rsquares, title = paste("Variance Explained in", feat_name), feat_name = NULL,
                   variable_name = names(Rsquares)[1], max_columns = 6){
  names(Rsquares)[1] <- variable_name
  if (length(Rsquares) == 1){
    warning("There is no influence of snps for this cpg")  
  }
  if (length(Rsquares) > max_columns){
    if (names(Rsquares[length(Rsquares)]) == "covariates"){
      vars <- sort(Rsquares[2:(length(Rsquares) - 1)], decreasing = TRUE)[1:4]
      Rsquares <- c(Rsquares[1], vars, Rsquares[length(Rsquares)])
    }else{
      vars <- sort(Rsquares[2:(length(Rsquares))], decreasing = TRUE)[1:5]
    }
  }else{}
  if (names(Rsquares[length(Rsquares)]) == "covariates"){
    if (length(Rsquares) == 2){
      colors <- 1:2
    }else{
      colors <- c(1, rep(4, length(Rsquares) - 1))
      colors[length(colors)] <- 2
      colors[which(Rsquares == max(Rsquares[2:(length(Rsquares) - 1)]))] <- 3
    }
    colorcodes <- c("#444444","blue",  "#CC3333", "#009900") 
  }else{
    colors <- c(1, rep(3, length(Rsquares) - 1))
    colors[which(Rsquares == max(Rsquares[2:length(Rsquares)]))] <- 2   
    colorcodes <- c("#444444","#CC3333", "#009900") 
  }
  colors <- as.factor(colors)
  Rsquares[Rsquares < 0] <- 0
  Rsquaresdf <- data.frame(variables = factor(names(Rsquares), levels = names(Rsquares)),
                           values = Rsquares, colors = colors)
  p <- ggplot2::ggplot(Rsquaresdf, ggplot2::aes_string(x = "variables", y = "values", fill = "colors")) + 
    ggplot2::geom_bar(stat = "identity") + ggplot2::scale_fill_manual(values = colorcodes) +
    ggplot2::theme(legend.position = "none") + ggplot2::scale_y_continuous("R squared", limits = c(0,1)) + 
    ggplot2::scale_x_discrete(name = "") + ggplot2::ggtitle(title)
  print(p)
}