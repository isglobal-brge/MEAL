#' Calculate RDA for a set
#' 
#' Perform RDA calculation for a \code{AnalysisRegionResults}. Feature values will
#' be considered the matrix X and phenotypes the matrix Y. Adjusting for covariates 
#' is done using covariable_names stored in the object.
#'  
#' @export RDAset
#' @param set \code{AnalysisResults}
#' @param varsmodel Matrix with the model
#' @param covarsmodel Matrix with the covariables model
#' @return Object of class \code{rda}
#' @seealso \code{\link[vegan]{rda}}
#' @examples
#' if (require(minfiData)){
#' set <- prepareMethylationSet(getBeta(MsetEx)[1:50,], pheno = pData(MsetEx))
#' methyOneVar <- DAPipeline(set, variable_names = "sex", probe_method = "ls")
#' rda <- RDAset(methyOneVar)
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