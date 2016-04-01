#' @describeIn AnalysisResults  Make a Volcano plot with the probe results
#' @param mindiff Numeric with the threshold to consider a difference in methylation or expression significant.
#' @aliases AnalysisResults-methods
setMethod(
  f = "plotVolcano",
  signature = "AnalysisResults",
  definition = function(object, variable = modelVariables(object)[1], mindiff = NULL){
    if (length(variable) != 1){
      stop("variable must have one value.")
    }  
    if (!is.character(variable)){
      stop("variable must be a character vector.")
    }
    if (!variable %in% modelVariables(object)){
      stop("Variable is not present in the model.")
    }
    
    results <- probeResults(object, drop = FALSE)[[variable]]
    results <- data.frame(beta = results[ , 2], logP = -log10(results[ , 5]), 
                          adj.p = results[ , 6], names = rownames(results))
    bonflevel <- -log10(0.05/nrow(results))
    if (is.null(mindiff)){
      if(object@originalclass == "ExpressionSet"){
        mindiff <- 1
      }      else{
        mindiff <- 0.1
      }
    }
    results$col <- ifelse(results$logP > bonflevel & abs(results$beta) > mindiff, 1, 2)
    
    p <- ggplot2::ggplot(results, ggplot2::aes_string(x = "beta", y = "logP", color = "col", label = "names")) +
      ggplot2::geom_point() + ggplot2::geom_vline(xintercept = mindiff, linetype = "dashed", col = 'red') +
      ggplot2::geom_vline(xintercept = -mindiff, linetype = "dashed", col = 'red') +
      ggplot2::scale_y_continuous(expression(~~-log[10](italic(p))), limits = c(0, NA)) + 
      ggplot2::ggtitle(paste("Volcano plot of", variable, "results")) +
      ggplot2::geom_hline(yintercept = bonflevel, linetype = "dashed", col = 'red') + 
      ggplot2::theme(legend.position = "none") 
    
    if (object@originalclass == "ExpressionSet"){
      p <- p + ggplot2::scale_x_continuous("Log Fold Change") +
        ggplot2::geom_text(ggplot2::aes(label = ifelse(col == 1, as.character(names),'')))
    }else{
      p <- p + ggplot2::scale_x_continuous("Change in methylation") +
        ggplot2::geom_text(ggplot2::aes(label = ifelse(col == 1, as.character(names),'')))
    }
    print(p) 
  }
)

