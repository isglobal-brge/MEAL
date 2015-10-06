#' @aliases MethylationSet-methods
#' @describeIn MethylationSet Filter probes with annotation
setMethod(
  f = "checkProbes",
  signature = "MethylationSet",
  definition = function(object) {
    bothprobes <- intersect(rownames(betas(object)), rownames(fData(object)))
    if(!length(bothprobes)){
      stop("There are no probes with annotation.")
    }      
    
    object <- object[bothprobes, ]
    annot <- fData(object)
    annot <- annot[bothprobes, , drop = FALSE]
    if (is(annot[, "chromosome"], "numeric")){
      annot[, "chromosome"] <- chrNumToChar(annot[, "chromosome"])
    }
    if (!is(annot[, "position"], "numeric")){
      annot[, "position"] <- as.numeric(annot[, "position"])
    }
    orderVector <- order(annot[, "chromosome"], annot[, "position"])
    object <- object[orderVector, ]
    annot <- annot[orderVector, ]
    fData(object) <-  annot
    featureNames(object) <- bothprobes[orderVector]  
    return(object)
  }
)
