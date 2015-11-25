#' @describeIn AnalysisRegionResults  Get the top features associated with the RDA
#' @aliases AnalysisRegionResults-methods
#' @param pval numeric with the p-value threshold. Only features with a p-values below this threshold will be shown.
setMethod(
  f = "topRDAhits",
  signature = "AnalysisRegionResults",
  definition =  function(object, pval = 0.05){
    
    rdavals <- getRDA(object)$CCA$u
    if (ncol(rdavals) > 2){
      rdavals <- rdavals[, 1:2]
    }
    vals <- t(featvals(object))
    resmat <- vapply(1:ncol(rdavals), function(x) vapply(1:ncol(vals), function(y){
      res <- cor.test(rdavals[,x], vals[,y])
      res <- c(res$estimate, res$p.value)
      names(res) <- c("cor", "pval")
      res
    }, numeric(2)), numeric(ncol(vals)*2))
    colnames(resmat) <- colnames(rdavals)
    cors <- resmat[seq(1, nrow(resmat), 2), , drop = FALSE]
    pvals <- resmat[seq(2, nrow(resmat), 2), , drop = FALSE]
    rownames(cors) <- colnames(vals)
    rownames(pvals) <- colnames(vals)
    adjpvals <- vapply(1:ncol(pvals), function(x) p.adjust(pvals[,x]), numeric(nrow(pvals)))
    signif.feats <- which(adjpvals < pval, arr.ind = TRUE)
    res <- data.frame(feat = rownames(cors)[signif.feats[ , 1]], RDA = colnames(cors)[signif.feats[ , 2]],
                      cor = cors[signif.feats], P.Value = pvals[signif.feats], adj.P.Val = adjpvals[signif.feats])
    res <- res[order(res$P.Value), ]
    res
  }
)