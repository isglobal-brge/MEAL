setMethod(
  f = "show",
  signature = "AnalysisResults",
  definition = function(object) {
    scat <- function(fmt, vals = character(), exdent = 2, ...) {
      vals <- ifelse(nzchar(vals), vals, "''")
      lbls <- paste(selectSome(vals), collapse = " ")
      txt <- sprintf(fmt, length(vals), lbls)
      cat(strwrap(txt, exdent = exdent, ...), sep = "\n")
    }
    cat("class:", class(object), "\n")
    cat("original class:", object@originalclass, "\n")
    features <- feats(object)
    scat("features(%d): %s\n", features)
    samples <- sampleNames(object)
    scat("samples(%d): %s\n", samples)
    cat("variables: ", paste(variableNames(object), collapse=", "), "\n")
    cat("model variables names: ", paste(modelVariables(object), collapse=", "), "\n")
    covars <- covariableNames(object)
    if (!length(covars)){
      covars <- "none"
    }
    cat("covariables:", paste(covars, collapse=", "), "\n")
    
    bumps <- bumps(object)
    cat("Number of bumps:", paste(sapply(bumps, nrow), collapse = ", "), "\n")
    
    blocks <- blocks(object)
    cat("Number of blocks:", paste(sapply(blocks, nrow), collapse = ", "), "\n")
    
    cate <- dmrCate(object)
    cat("Number of regions:", paste(sapply(cate, nrow), collapse = ", "), "\n")
    
    probe <- probeResults(object, drop = FALSE)
    if (object@originalclass == "MethylationSet"){
      text <- "methylated"
    }else{
      text <- "expressed"
    }
    cat("Number of differentially", text,"probes:", 
        paste(sapply(probe, function(x) sum(x$adj.P.Val < 0.05)), collapse = ", "), "\n")
  }
)
