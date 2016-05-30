#' Detect regions differentially methylated
#' 
#' This function is a wrapper of two known region differentially methylated detection 
#' methods: \emph{Bumphunter} and \emph{DMRcate}. \code{blockFinder} implementation 
#' present in \code{minfi} package is also available.
#' 
#' @export DARegion
#' @param set \code{MethylationSet}. 
#' @param model Model matrix representing a linear model.
#' @param proberes Data.frame or list of data.frames with the results of 
#' \code{DAProbe}
#' @param methods Character vector with the names of the methods used to estimate 
#' the regions. Valid names are: "blockFinder", "bumphunter" and "DMRcate".
#' @param coefficient Numeric with the index of the model matrix used to perform
#' the analysis.
#' @param num_permutations Numeric with the number of permutations used to
#' calculate p-values in \code{bumphunter} and \code{blockFinder}
#' @param bumphunter_cutoff Numeric with the threshold to consider a probe
#' significant. If one number is supplied, the lower limit is minus the upper one.
#' If two values are given, they will be upper and lower limits.
#' @param bumps_max Numeric with the maximum number of bumps allowed.
#' @param num_cores Numeric with the number of cores used to perform the permutation.
#' @param verbose Logical value. If TRUE, it writes out some messages indicating progress. 
#' If FALSE nothing should be printed.
#' @param ... Further arguments passsed to \code{bumphunter} function.
#' @return List with the main results of the three methods. If a method is not
#' chosen, NA is returned in this position. 
#'  
#' @details \code{DARegion} performs a methylation region analysis using 
#' \emph{bumphunter} and \emph{DMRcate}. Bumphunter allows the modification of several 
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
#' linear regression and not robust linear regression. The results of \code{DAProbe}
#' can be supplied in proberes argument, skipping this first step.
#'  
#' \code{DARegion} supports multiple variable analyses. If coefficient is a vector,
#' a list of lists will be returned. Each member will be named after the name of the
#' column of the model matrix.
#'  
#' @seealso \code{\link[minfi]{bumphunter}}, \code{\link[minfi]{blockFinder}}, 
#' \code{\link[DMRcate]{dmrcate}} 
#' @examples
#' if (require(minfiData)){
#'  set <- prepareMethylationSet(minfi::getBeta(MsetEx)[1:10, ], pheno = pData(MsetEx))
#'  model <- model.matrix(~Sample_Group, data = pData(MsetEx)) 
#'  res <- DARegion(set, model) 
#'  res
#' }

DARegion <- function (set, model, proberes, methods = c("blockFinder", "bumphunter", "DMRcate"), 
                      coefficient = 2, num_permutations = 0, bumphunter_cutoff = 0.05, 
                      bumps_max = 30000, num_cores = 1, verbose = FALSE, ...)
{
  if (sum(methods %in% c("blockFinder", "bumphunter", "DMRcate", "none")) == 0){
    stop("Method variable is empty or none of the methods introduced is valid.")
  }
  if (nrow(set) == 0 | ncol(set) == 0){
    stop("The set is empty.")
  }
  if (nrow(model) == 0 | ncol(set) == 0){
    stop("The model matrix is empty.")
  }
  if (ncol(set) != nrow(model)){
    stop("The number of samples is different in the set and in the model")
  }
  if (!is(set, "MethylationSet")){
    stop("set must be a MethylationSet.")
  }
  if (missing(proberes)){
    proberes <- rep(NA, length(coefficient))
  }
  if (length(coefficient) > 1){
    regions <- lapply(coefficient, 
                      function(x) DARegion(set = set, model = model, methods = methods, 
                                           coefficient = x, proberes = proberes[[which(coefficient == x)]],
                                           num_permutations = num_permutations, 
                                           bumphunter_cutoff = bumphunter_cutoff,
                                           bumps_max = bumps_max, num_cores = num_cores, 
                                           verbose = verbose, ...))
    names(regions) <- colnames(model)[coefficient]
    return(regions)
  }
  #Activate parallelization
  if (num_cores > 1){
    doParallel::registerDoParallel(cores = num_cores)
  }
  dmrs <- list()
  if("bumphunter" %in% methods){
    M <- getMs(set)
    annot <- fData(set)
    Bumphunter <-  minfi::bumphunter(object = M, design = model, coef = coefficient,
                                     chr = annot[, "chromosome"], pos = annot[, "position"],
                                     cutoff = bumphunter_cutoff, 
                                     nullMethod = "bootstrap", verbose = verbose, ...)$table
    
    #Assure that permutations are applied to a reasonable number of bumps
    if (num_permutations > 0){
      i <- 1
      while(nrow(Bumphunter) > bumps_max){
        if (verbose){
          message(paste("Iteration",i,"Num bumps:", nrow(Bumphunter), 
                        "cutoff:", bumphunter_cutoff))
        }
        bumphunter_cutoff <- bumphunter_cutoff + 0.05
        Bumphunter <-  minfi::bumphunter(object = M, design = model, coef = coefficient,
                                         chr = annot[ , "chromosome"], pos = annot[, "position"],
                                         cutoff = bumphunter_cutoff, B = 0, 
                                         nullMethod = "bootstrap", verbose = verbose, ...)$table
        i <- i +1
      }
      if (verbose){
        message(paste("Iteration",i,"Num bumps:", nrow(Bumphunter), 
                      "cutoff:", bumphunter_cutoff))
      }
      Bumphunter <-  minfi::bumphunter(object = M, design = model, coef = coefficient,
                                       chr = annot[, "chromosome"], pos = annot[, "position"],
                                       cutoff = bumphunter_cutoff, B = num_permutations, 
                                       nullMethod = "bootstrap", verbose = verbose, ...)$table
    } 
    if (length(Bumphunter) == 1){
      Bumphunter <- data.frame()
    }
    dmrs[["bumphunter"]] <- Bumphunter
  }else{
    dmrs[["bumphunter"]] <- NA
  }
  if("blockFinder" %in% methods){
    #Create a GenomicRatioSet from minfi to run cpgCollapse
    granges <- fData(set)[, c("chromosome", "position", "position")]
    granges <- createRanges(granges)
    beta <- betas(set)
    minfiset <- minfi::GenomicRatioSet(gr = granges, Beta = beta, pData = pData(set),
                                       CN = matrix(1, ncol = ncol(beta), nrow = nrow(beta)),
                                       annotation = c(array = "IlluminaHumanMethylation450k",
                                                      annotation = "ilmn12.hg19"))
    blockset <- minfi::cpgCollapse(minfiset, returnBlockInfo = FALSE, verbose = verbose )
    blockfinder <- tryCatch(minfi::blockFinder(blockset, design = model, coef = coefficient,
                                               cutoff = bumphunter_cutoff, nullMethod = "bootstrap",
                                               verbose = verbose, ...)$table, 
                            error = function(e) NULL)
    if (is.null(blockfinder)){
      blockfinder <- data.frame()
      dmrs[["blockFinder"]] <- blockfinder
    }else{
      
      #Assure that permutations are applied to a reasonable number of bumps
      if (num_permutations){
        i <- 1
        while(nrow(blockfinder) > bumps_max){
          message(paste("Iteration",i,"Num bumps:", nrow(blockfinder), 
                        "cutoff:", bumphunter_cutoff))
          bumphunter_cutoff <- bumphunter_cutoff + 0.05
          blockfinder <- minfi::blockFinder(blockset, design = model, coef = coefficient, 
                                            cutoff = bumphunter_cutoff, B = 0, 
                                            verbose = verbose, ...)$table
          i <- i +1
        }
        if (verbose){
          message(paste("Iteration",i,"Num bumps:", nrow(blockfinder), 
                        "cutoff:", bumphunter_cutoff))
        }
        blockfinder  <- minfi::blockFinder(blockset, design = model, coef = coefficient,
                                           cutoff = bumphunter_cutoff, B = num_permutations, 
                                           nullMethod = "bootstrap",
                                           verbose = verbose, ...)$table
      } 
      if (length(blockfinder) == 1){
        blockfinder <- data.frame()
      }
      dmrs[["blockFinder"]] <- blockfinder
    }
  }else{
    dmrs[["blockFinder"]] <- NA
  }
  if("DMRcate" %in% methods){
    if (!is(proberes, "data.frame")){
      proberes <- DAProbe(set = set, model = model, coefficient = coefficient,
                          method = "ls")
    }
    myannotation <- DMRcate::cpg.annotate(datatype = "array", object = MultiDataSet::getMs(set), 
                                          design = model, coef = coefficient)
    dmrcoutput <- tryCatch(DMRcate::dmrcate(myannotation, lambda = 1000, C = 2), 
                           error = function(e) NULL)
    if (is.null(dmrcoutput)){
      dmrs[["DMRcate"]] <-  data.frame()
    }else {
      dmrs[["DMRcate"]] <- dmrcoutput$results
    }
  }else{
    dmrs[["DMRcate"]] <- NA
  }
  dmrs
}

dmrcateCreator <- function(proberesults){
  dmrcate <- list()
  proberesults <- proberesults[order(proberesults$chr, proberesults$pos), ]
  dmrcate[["ID"]] <- rownames(proberesults)
  dmrcate[["stat"]] <- proberesults$t
  dmrcate[["CHR"]] <- proberesults$chromosome
  dmrcate[["pos"]] <- proberesults$position
  dmrcate[["gene"]] <- proberesults$genes
  dmrcate[["group"]] <- proberesults$group
  dmrcate[["betafc"]] <- proberesults[, 6]
  dmrcate[["indfdr"]] <- proberesults$adj.P.Val
  dmrcate[["is.sig"]] <- proberesults$adj.P.Val < 0.05
  class(dmrcate) <- "annot"
  return(dmrcate)
}