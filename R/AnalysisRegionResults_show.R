setMethod(
  f = "show",
  signature = "AnalysisRegionResults",
  definition = function(object) {
    callNextMethod()
    scat <- function(fmt, vals = character(), exdent = 2, ...) {
      vals <- ifelse(nzchar(vals), vals, "''")
      lbls <- paste(selectSome(vals), collapse = " ")
      txt <- sprintf(fmt, length(vals), lbls)
      cat(strwrap(txt, exdent = exdent, ...), sep = "\n")
    }
    range <- getRange(object)
    cat("Range:\n\tChr: ", as.character(S4Vectors::runValue(GenomicRanges::seqnames(range))),"\tstart: ", 
        GenomicRanges::start(range),"\tend: ", GenomicRanges::end(range),"\n")
    cat("R2: ", round(regionR2(object), 3), "\n")
    pval <- RDAPval(object)
    if (pval < 1e-3){
      cat("RDA P-value: ", format(signif(pval, 3), scientific = TRUE),"\n")
    }else{
      cat("RDA P-value: ", signif(pval, 3), "\n")
    }
    globPval <- globalPval(object)
    if (globPval < 1e-3){
      cat("Region P-value: ", format(signif(globPval, 3), scientific = TRUE),"\n")
    }else{
      cat("Region P-value: ", signif(globPval, 3), "\n")
    }
    cat("Global R2: ", round(globalR2(object), 3), "\n")
    
  }
)
