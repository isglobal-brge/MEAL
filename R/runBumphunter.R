#' Run bumphunter
#' 
#' Run bumphunter to a methylation dataset. This function contains all steps of
#' bumphunter analysis, from model.matrix creation to running the analysis. 
#' 
#' @details runBumphunter is a wrapper for minfi \code{bumphunter}. This function 
#' runs all the steps required prior running bumphunter from the methylation set
#' and the formula of the model. This implementation allows running bumphunter to other 
#' objects than \code{GenomicRatioSet}. The result can be encapsulated in a 
#' \code{ResultSet} to take adavantege of its plotting capabilities.
#' 
#' If the user wants to run permutations to calculate p-values, this implementation
#' can filter the bumps to avoid doing a very high number of permutations and to
#' reduce computation time. To do so, we can set the maximum number of bumps that 
#' we want to permute with the \code{bumps_max} parameter. runBumphunter increases
#' \code{bumphunter_cutoff} value until the number of bumps is lower than 
#' \code{bumps_max}.
#' 
#' @export
#' 
#' @param set \code{GenomicRatioSet}, \code{eSet} derived object or 
#' \code{SummarizedExperiment}
#' @param model Model matrix or formula to get model matrix from \code{set}. 
#' @param coefficient Numeric with the column of model matrix used in the analysis.
#' (Default: 2)
#' @param bumphunter_cutoff Numeric with the minimum cutoff to include a probe
#' in a block. (Default: 0.1)
#' @param num_permutations Numeric with the number of permutations run to compute 
#' the bumps p-value. (Default: 0)
#' @param betas If \code{set} is a \code{GenomicRatioSet}, should beta values be
#' used? (Default: TRUE)
#' @param check_perms Logical. Should we check that there are less bumps than 
#' \code{bumps_max}? This parameter only applies when \code{num_permutations} is 
#' greater than 0. (Default: TRUE)
#' @param bumps_max Numeric with the maximum number of bumps used in the permutation.
#' This parameter only applies when \code{num_permutations} is greater than 0. 
#' (Default: 30000)
#' @param verbose Logical value. Should the function be verbose? (Default: FALSE)
#' @param resultSet Should results be encapsulated in a \code{resultSet}? (Default: TRUE)
#' @param ... Further arguments passed to \code{bumphunter}.
#' @return data.frame or \code{resultSet} with the result of \code{bumphunter}
#' @seealso \code{\link[minfi]{bumphunter}}
runBumphunter <- function(set, model, coefficient = 2, bumphunter_cutoff = 0.1,
                          num_permutations = 0, bumps_max = 30000, betas = TRUE, 
                          check_perms = TRUE, verbose = FALSE, 
                          resultSet = FALSE, ...){
  #Activate parallelization 
  ### Change to BiocParallel !!!!!!!
  # if (num_cores > 1){
  #   doParallel::registerDoParallel(cores = num_cores)
  # }
  # 
  ## Create model matrix from formula
  if (is(model, "formula")){
    model <- createModel(set, model)
    set <- set[, rownames(model)]
  }
  

  ## Get matrix
  if (is(set, "GenomicRatioSet")){
    mat <- minfi::getBeta(set)
    if (!betas) {
      mat[mat == 0] <- 1e-3
      mat[mat == 1] <- 1 - 1e-3
      mat <- minfi::logit2(mat)
    }
  } else if (is(set, "SummarizedExperiment")){
    mat <- Biobase::assays(set)
  } else {
    stop("set must be a MethylationSet, GenomicRatioSet or SummarizedExperiment.")
  }
  
  if (is(set, "eSet")){
    fFun <- function(set) {
      df <- Biobase::fData(set)
      ## Change position by start to harmonize fDatas
      if ("position" %in% colnames(df)){
        colnames(df)[colnames(df) == "position"] <- "start"
      } 
      df
    }
  } else if (is(set, "RangedSummarizedExperiment")){
    fFun <- function(set) { 
      df <- as.data.frame(SummarizedExperiment::rowRanges(set))
      colnames(df)[1] <- "chromosome"
      df
    }
  } else if (is(set, "SummarizedExperiment")){
    fFun <- SummarizedExperiment::rowData
  } else if (is.matrix(set)){
    fFun <- function(set) data.frame(matrix(vector(), nrow(set), 0))
  }
  annot <- fFun(set)
  
  if (check_perms == FALSE | num_permutations == 0){
    res <-  minfi::bumphunter(object = mat, design = model, coef = coefficient,
                                     chr = annot[, "chromosome"], pos = annot[, "start"],
                                     cutoff = bumphunter_cutoff, B = num_permutations,
                                     nullMethod = "bootstrap", verbose = verbose, ...)$table
  } else {
  # Ensure that permutations are applied to a reasonable number of bumps
    i <- 1
    res <-  minfi::bumphunter(object = mat, design = model, coef = coefficient,
                              chr = annot[ , "chromosome"], pos = annot[, "start"],
                              cutoff = bumphunter_cutoff, B = 0, 
                              nullMethod = "bootstrap", verbose = verbose, ...)$table
    if (verbose){
      message(paste("Iteration",i,"Num bumps:", nrow(res), 
                    "cutoff:", bumphunter_cutoff))  
    }
      ## Increase cut off until getting a reasonable number of bumps
      while(nrow(res) > bumps_max){
        i <- i +1
        bumphunter_cutoff <- bumphunter_cutoff + 0.05
        res <-  minfi::bumphunter(object = mat, design = model, coef = coefficient,
                                         chr = annot[ , "chromosome"], pos = annot[, "start"],
                                         cutoff = bumphunter_cutoff, B = 0, 
                                         nullMethod = "bootstrap", verbose = verbose, ...)$table
        if (verbose){
          message(paste("Iteration",i,"Num bumps:", nrow(res), 
                        "cutoff:", bumphunter_cutoff))
        }
      }
   } 
    if (length(res) == 1){
      res <- data.frame()
    }
  if (resultSet)
  {
    fFun <- getFeatureDataFun(set)
    
    res <- create_resultset("runBumphunter", lResults = 
                              list(bumphunter = list(result = res, error = NA)),  
                            fData = list(main = fFun(set)), lOptions = list())
  }
  res
}