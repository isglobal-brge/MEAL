#' @describeIn AnalysisResults  Plot of the region
#' @aliases AnalysisResults-methods 
setMethod(
  f = "plotRegion",
  signature = "AnalysisResults",
  definition = function(object, variable = modelVariables(object)[[1]], range = NULL){
    if (!"position" %in% colnames(probeResults(object)[[1]])){
      stop("Results must have a column called position to perform this plot.")
    }
    if (length(variable) != 1){
      stop("variable must have one value.")
    }  
    if (!variable %in% modelVariables(object)){
      stop(paste("Variable is not present in modelVariables.\nValid variables are", 
                 paste(modelVariables(object), collapse = ", ")))
    }
    
    results <- probeResults(object)[[variable]]
    
    if (!is(object, "AnalysisRegionResults")){
      if (is.null(range)){
        stop("range must be present to use plotRegion with a AnalysisResults")
      }
      results <- filterResults(results, range)
    }
    
    results <- data.frame(beta = results[ , 2], positions = results[ , 8], P.Value = results[ , 5])
    p <- ggplot2::ggplot(results, ggplot2::aes_string(x = "positions", y = "beta", color = "P.Value" < 0.05)) +
      ggplot2::geom_point(alpha = 1/2) + ggplot2::geom_line(ggplot2::aes(y = 0.05, color = "green")) +
      ggplot2::geom_line(ggplot2::aes(y = -0.05, color = "green")) +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_y_continuous(expression(paste("Change in methylation (", Delta, "beta)")))
    print(p)
  }
)