#' Calculate RDA for a set
#' 
#' Perform RDA calculation for a \code{AnalysisRegionResults}. Feature values will
#' be considered the matrix X and phenotypes the matrix Y. Adjusting for covariates 
#' is done using a model matrix passed in covarsmodel.
#'  
#' @export RDAset
#' @param set \code{MethylationSet}, \code{ExpressionSet} or \code{matrix}
#' @param varsmodel Matrix with the model
#' @param covarsmodel Matrix with the covariables model
#' @return Object of class \code{rda}
#' @seealso \code{\link[vegan]{rda}}
#' @examples
#' if (require(minfiData)){
#' set <- prepareMethylationSet(getBeta(MsetEx)[1:50,], pheno = data.frame(pData(MsetEx)))
#' model <- model.matrix(~set$age)
#' rda <- RDAset(set, model)
#' rda
#' }
RDAset <- function(set, varsmodel = NULL, covarsmodel = NULL){
        
        if (ncol(set) == 0 | nrow(set) == 0){
                stop("The set is empty.")
        }
        if (ncol(varsmodel) == 0 | nrow(varsmodel) == 0){
                stop("The model matrix is empty.")
        }
        if (ncol(set) != nrow(varsmodel)){
                stop("The number of samples is different in the set and in the model")
        }
        
        if (is(set, "MethylationSet")){
                msg <- validObject(set, test = TRUE)
                if (is(msg, "character")){
                        message(paste(msg, collapse = "\n"))
                        stop("checkProbes and checkSamples might solve validity issues.")
                }
                vals <- t(MultiDataSet::betas(set))
        } 
        else if (is(set, "matrix")){
                vals <- t(set)
        } else if (is(set, "ExpressionSet")){
                vals <- t(exprs(set))
        } else{
                stop("set must be a MethylationSet, an ExpressionSet or a matrix.")
        }
        
        
        rdaRes <- vegan::rda(vals, varsmodel, covarsmodel)
        
        return(rdaRes)
}