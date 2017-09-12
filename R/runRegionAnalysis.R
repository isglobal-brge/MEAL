#' Run different DMR detection methods
#' 
#' This function is a wrapper of two known region differentially methylated detection 
#' methods: \emph{Bumphunter}, \code{blockFinder} and \emph{DMRcate}.
#' 
#' @export runRegionAnalysis
#' @details \code{runRegionAnalysis} performs a methylation region analysis using 
#' \emph{bumphunter}, \code{blockFinder} and \emph{DMRcate}. Bumphunter allows the modification of several 
#' parameters that should be properly used. 
#'  
#' Cutoff will determine the number of bumps that will be detected. The smaller the cutoff, the higher the
#' number of positions above the limits, so there will be more regions and they
#' will be greater. \code{Bumphunter} can pick a cutoff using the null distribution,
#' i.e. permutating the samples. There is no standard cutoff and it will depend
#' on the features of the experiment. Permutations are used to estimate p-values and, 
#' if needed, can be used to pick a cutoff. The advised number of permutation is 1000. 
#' The number of permutations will define the maximum number of bumps that will be considered
#' for analysing. The more bumps, the longer permutation time. As before,
#' there is not an accepted limit but \code{minfi} tutorial recommends not to exceed
#' 30000 bumps. Finally, if supported, it is very advisable to use parallelization
#' to perform the permutations.
#'  
#' Due to \code{minfi} design, \emph{BlockFinder} can only be run using own minfi
#' annotation. This annotation is based on hg19 and Illumina 450k chipset. Cpg sites 
#' not named like in this annotation package will not be included. As a result,
#' the use of \emph{BlockFinder} is not recommended.
#' 
#' \emph{DMRcate} uses a first step where linear regression is performed in order
#' to estimate coefficients of the variable of interest. This first step is equal
#' to the calculation performed in \code{DAProbe}, but using in this situation
#' linear regression and not robust linear regression. 
#'  
#' 
#' @param set \code{GenomicRatioSet}, \code{eSet} derived object or 
#' \code{SummarizedExperiment}
#' @param model Model matrix representing a linear model.
#' @param methods Character vector with the names of the methods used to estimate 
#' the regions. Valid names are: "blockFinder", "bumphunter" and "DMRcate".
#' @param coefficient Numeric with the index of the model matrix used to perform
#' the analysis.
#' @param bumphunter_params List with other parameter passed to \code{runBumphunter} 
#' function.
#' @param blockFinder_params List with other parameter passed to \code{runBlockFinder} 
#' function.
#' @param dmrcate_params List with other parameter passed to \code{runDMRcate} 
#' function.
#' @param verbose Logical value. Should the function be verbose? (Default: FALSE)
#' @param resultSet Should results be encapsulated in a \code{resultSet}? (Default: TRUE)
#' @return List or \code{resultSet} with the result of the DMR detection methods.
#' @seealso \code{\link[minfi]{bumphunter}}, \code{\link[minfi]{blockFinder}}, 
#' \code{\link[DMRcate]{dmrcate}} 
#' @examples
#' if (require(minfiData)){
#'  set <- prepareMethylationSet(minfi::getBeta(MsetEx)[1:10, ], pheno = data.frame(pData(MsetEx)))
#'  model <- model.matrix(~Sample_Group, data = pData(MsetEx))
#'  res <- DARegion(set, model) 
#'  res
#' }
runRegionAnalysis <- function (set, model, methods = c("blockFinder", "bumphunter", "DMRcate"), 
                            coefficient = 2, bumphunter_params = NULL, 
                            blockFinder_params = NULL, dmrcate_params = NULL,
                            verbose = FALSE, resultSet = TRUE)
{
  methods <- match.arg(methods, choices = c("blockFinder", "bumphunter", "DMRcate"),
                       several.ok = TRUE)
  
  dmrs <- list(bumphunter = list(result = data.frame(), error = NA), 
               blockFinder = list(result = data.frame(), error = NA), 
               dmrcate = list(result = data.frame(), error = NA))
  if("bumphunter" %in% methods) {
    bumphunter <- tryCatch(do.call(runBumphunter, c(list(set = set, model = model, 
                           coefficient = coefficient, verbose = verbose), bumphunter_params)), 
                           error = function(e) as.character(e))
    if (is.data.frame(bumphunter)){
      dmrs$bumphunter$result <- bumphunter
    } else {
      dmrs$bumphunter$error <- bumphunter
    }
  }
  if("blockFinder" %in% methods) {
    blockFinder <- tryCatch(do.call(runBlockFinder, c(list(set = set, model = model, 
                            coefficient = coefficient, verbose = verbose), blockFinder_params)), 
                            error = function(e) as.character(e))
    if (is.data.frame(blockFinder)){
      dmrs$blockFinder$result <- blockFinder
    } else {
      dmrs$blockFinder$error <- blockFinder
    }
  }
  if("DMRcate" %in% methods) {
    dmrcate <- tryCatch(do.call(runDMRcate,  c(list(set = set, model = model, 
                              coefficient = coefficient, verbose = verbose), dmrcate_params)),
                        error = function(e) as.character(e))
    if (is.data.frame(dmrcate)){
      dmrs$dmrcate$result <- dmrcate
    } else {
      dmrs$dmrcate$error <- dmrcate
    }
  } 
  if (resultSet){
    dmrs <- create_resultset(fOrigin = "RegionAnalysis", lResults = dmrs, 
                             fData = list(main = data.frame()))
  }
  dmrs
}