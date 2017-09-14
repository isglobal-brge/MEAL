#' Get the top features associated with the RDA
#' 
#' Get a list of the features significantly associated to the first two RDA components
#' 
#' @name topRDAhits
#' @export
#' 
#' @param object \code{ResultSet} 
#' @param tPV numeric with the p-value threshold. Only features with a p-values below this threshold will be shown. 
#' @return data.frame with the features, the component, the correlation and the p-value
#' @examples
#' if (require(minfiData) & require(GenomicRanges)){
#' set <- ratioConvert(mapToGenome(MsetEx[1:10,]))
#' model <- model.matrix(~set$sex)
#' rda <- runRDA(set, model)
#' topRDAhits(rda)
#' }
topRDAhits <- function(object, tPV = 0.05){
    
    stopifnot("RDA" %in% names(object))
    rda <- getAssociation(object, rid = "RDA")
    rdavals <- rda$CCA$u
    
    if (ncol(rdavals) > 2){
      rdavals <- rdavals[, 1:2]
    }
    vals <- vegan::ordiYbar(rda, "partial")
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
    signif.feats <- which(adjpvals < tPV, arr.ind = TRUE)
    res <- data.frame(feat = rownames(cors)[signif.feats[ , 1]], 
                      RDA = colnames(cors)[signif.feats[ , 2]],
                      cor = cors[signif.feats], 
                      P.Value = pvals[signif.feats], 
                      adj.P.Val = adjpvals[signif.feats])
    res <- res[order(res$P.Value), ]
    res
  }