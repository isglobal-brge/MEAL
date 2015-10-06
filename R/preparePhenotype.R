#' Process a table of phenotypes
#' 
#' Given a data.frame containing phenotypic variables, select the desired columns
#' and transform them to the desired types. 
#' 
#' @export preparePhenotype
#' 
#' @param phenotypes Data.frame with the phenotypic features
#' @param variable_names Vector with the names or the positions of the desired
#' variables.
#' @param variable_types Vector with the types of the variables. 
#' @details \code{preparePhenotype} supports five types of variables. Categorical 
#' and continuous correspond to factor and numerical types in R. The other three
#' are genomic models as defined in \code{SNPassoc}: dominant, recessive and additive.
#' In order to use these types, only two alleles can be present and genotypes should
#' be specified in the form \emph{a/b}. 
#' 
#' If transformation of variables is not needed, the variable_types can be passed
#' as a vector of \code{NA}.
#' @return Data.frame with the columns selected and with the types desired.
#' @examples
#' pheno <- data.frame(a = sample(letters[1:2], 5, replace = TRUE), b = runif(5), 
#' c = sample(c("a/a","a/b", "b/b"), 5, replace = TRUE))
#' pheno <- preparePhenotype(pheno, variable_names = c("a", "c"), 
#' variable_types = c("categorical", "dominant"))
#' pheno

preparePhenotype <- function(phenotypes, variable_names, variable_types = rep(NA, length(variable_names))){
  if (is.null(dim(phenotypes))){ 
    stop("Phenotypes is not a valid data.frame")
  }
  if (ncol(phenotypes) == 0 || nrow(phenotypes) == 0){
    stop("Phenotypes data.frame is empty.")
  }
  if (length(variable_names) == 0){
    stop("variable_names is empty.")
  }
  if (length(variable_types) == 0){
    stop("variable_types is empty.")
  }
  if (!is.null(dim(variable_names))){
    stop("variable_names must be a vector.")
  }
  if (!is.null(dim(variable_types))){
    stop("variable_types must be a vector.")
  }
  if (is(variable_names, "character")){
    logical_names <- variable_names %in% colnames(phenotypes)
    if (!all(logical_names)){
      stop(paste(paste(variable_names[!logical_names], collapse = ", "), "are not columns of the phenotypes table."))
    }
  }else if (is(variable_names, "numeric")){
    logical_names <- variable_names >= 1 & variable_names <= ncol(phenotypes)
    if (!all(logical_names)){
      stop(paste("variable_names indexes must be between 1 and", ncol(phenotypes)))
    }
  }
  
  #Remove rows with NAs from phenotypes matrix
  phenotypes <- phenotypes[, variable_names, drop = FALSE]
  phenotypes <- phenotypes[apply(!is.na(phenotypes), 1, all), ,drop=FALSE]
  if (length(variable_names) > length(variable_types)){
    warning("The number of types of variables is smaller than the number of variables. Variables without type won't be changed.")
    variable_types = c(variable_types, 
                       rep(NA, length(variable_names) - length(variable_types)))
  }
  temp <- do.call(data.frame, lapply(1:ncol(phenotypes), 
                                     function(x) changeVectorType(phenotypes[,x], variable_types[x])))
  colnames(temp) <- colnames(phenotypes)
  rownames(temp) <- rownames(phenotypes)
  temp
}