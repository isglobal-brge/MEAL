#' @describeIn AnalysisResults  Plot a Manhattan plot with the probe results
#' @aliases AnalysisResults-methods
#' @param variable Character with the variable name used to obtain the probe results.
#' Note: model name should be used. Original variable name might not be valid.  
#' @param range GenomicRange whose probes will be highlighted
setMethod(
  f = "plotEWAS",
  signature = "AnalysisResults",
  definition = function(object, variable = modelVariables(object)[1], range = NULL){
    if (is(object, "AnalysisRegionResults")){
      warning("For a better visualization of the region, use plotRegion.")
      range <- NULL
    }
    if (length(variable) != 1){
      stop("variable must have one value.")
    }  
    if (!variable %in% modelVariables(object)){
      stop("Variable is not present in the model.")
    }
    
    dmp <- probeResults(object)[[variable]]
    
    dmp$CHR <- gsub("chr", "", dmp$chromosome)
    dmp$CHR[dmp$CHR %in% c("X","x")] <- 23
    dmp$CHR[dmp$CHR %in% c("Y","y")] <- 24
    dmp$CHR[dmp$CHR == "XY"] <- 25
    dmp$CHR[dmp$CHR == "MT"] <- 26
    dmp$CHR <- as.numeric(dmp$CHR)
    if ("start" %in% colnames(dmp)){
      dmp$position <- dmp$start
    }
    dmp <- dmp[order(dmp$CHR, dmp$position), ]
    dmp$idx <- seq_len(nrow(dmp))
    dmp$logP <- -log10(dmp$P.Value)
    dmp$col <- dmp$CHR %% 2 + 1
    if (!is.null(range)){
      selectedcpgs <- rownames(filterResults(dmp, range))
      dmp[selectedcpgs, "col"] <- 3
    }
    dmp$col <- as.factor(dmp$col)
    breaks <- vapply(split(dmp, dmp$CHR), 
                     function(x) round((max(x$idx) - min(x$idx))/2) + min(x$idx), numeric(1))
    cbPalette <- c("#999999", "#000000", "#33FF33")
    bonflevel <- -log10(0.05/nrow(dmp))
    p <- ggplot2::ggplot(dmp, ggplot2::aes_string(x = "idx", y = "logP", color = "col")) + 
      ggplot2::geom_point(aes(fill = col)) + 
      ggplot2::scale_x_continuous("Chromosome", labels = unique(dmp$chr), 
                                  breaks = breaks, limits = c(1, max(dmp$idx))) +
      ggplot2::scale_y_continuous(expression(~~-log[10](italic(p)))) +
      ggplot2::theme(legend.position = "none",  axis.text.x  = ggplot2::element_text(angle = 45, size = 8)) + 
      ggplot2::geom_hline(yintercept = bonflevel, linetype = 1, col = 'red') +
      ggplot2::scale_colour_manual(values = cbPalette)
    print(p)
  }
)

