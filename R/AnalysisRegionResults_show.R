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
    relsnps <- snps(object)
    scat("Relevant snps(%d): %s\n", relsnps)
    cat("Snps Variance: ", snpsVar(object), "\n")
    range <- getRange(object)
    cat("Range:\n\tChr: ", as.character(S4Vectors::runValue(GenomicRanges::seqnames(range))),"\tstart: ", 
        GenomicRanges::start(range),"\tend: ", GenomicRanges::end(range),"\n")
    cat("Rsquared: ", regionR2(object),"\n")
    cat("P-value: ", regionPval(object),"\n")
  }
)
