#' Run bumphunter
#' 
#' @details This function has been deprecated and will be defunct in the new version.
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
                          check_perms = FALSE, verbose = FALSE, 
                          resultSet = FALSE, ...){

  .Deprecated()
}