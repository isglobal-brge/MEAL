#' @describeIn AnalysisResults  QQ plot of probe analysis
#' @aliases AnalysisResults-methods 
setMethod(
  f = "plotQQ", 
  signature = "AnalysisResults",
  definition = function(object, variable = modelVariables(object)[1]){
  if (is(object, "AnalysisRegionResults")){
    warning("QQplot is not recommended in range analyis.")
  }
  if (length(variable) != 1){
    stop("variable must have one value.")
  }  
  if (!variable %in% modelVariables(object)){
    stop("Variable is not present in phenodata.")
  }
  qqplotBand(probeResults(object)[[variable]]$P.Value)  
}
) 

qqplotBand <- function(x, main = NULL, alpha = 0.1)
{
  o = -log10(sort(x, decreasing = FALSE))
  e = -log10(1:length(o)/length(o))
  df <- data.frame(o = o, e = e)
  names(o) <- names(x)
  names(e) <- names(x)
  
  N <-length(x)
  conf <- 1 - (alpha/2)
  c95 <- mapply(qbeta, conf, 1:N, N:1)
  c05 <- mapply(qbeta, 1 - conf, 1:N, N:1)
  band <- cbind(c05, c95)
  pol <- data.frame(x = c(e, e[N:1]), y = c(-log10(band[ , 1]), -log10(band[N:1, 2])))
  p <- ggplot2::ggplot(df, ggplot2::aes_string(x = "e", y = "o")) + 
    ggplot2::geom_polygon(data = pol, ggplot2::aes_string(x = "x", y = "y"), fill = "skyblue")
  
  p <-  p + ggplot2::geom_point() + 
    ggplot2::scale_x_continuous(expression(Expected ~~-log[10](italic(p)))) + 
    ggplot2::scale_y_continuous(expression(Observed ~~-log[10](italic(p)))) + 
    ggplot2::ggtitle(main) + ggplot2::geom_smooth(ggplot2::aes_string(y = "e"), method = lm)
  
  print(p)

}