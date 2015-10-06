#' Calculate R2 for different variables
#' 
#' Using a data.frame as input, calculates the R2 between a dependent variable and
#' some independent variables. Base adjusting by covariates can also be used.
#' 
#' @details \code{explainedVariance} computes R2 via linear models. The first column
#' is considered to be the dependent variable. Therefore, a lineal model will be
#' constructed for each of the remaining variables. In case that covariates 
#' were included, they will be included in all the models and, in addition, a model
#' containing only the covariates will be returned.
#' 
#' Some variables can be grouped in the models to assess their effect together.
#' 
#' @export
#' @param data Data.frame containing the dependent variable in the first column.
#' @param num_mainvar Numerical with the number of variables that should be grouped.
#' They should be at the beggining.
#' @param num_covariates Numerical with the number of variables that should be
#' considered as covariates. Covariates variables must be at the end.
#' @param variable_label Character with the name of the main variable in the results.
#' @return Numeric vector with the R2 explained by each of the variables.
#' @examples
#' data(mtcars)
#' R2 <- explainedVariance(mtcars)
#' R2
explainedVariance <- function(data, num_mainvar = 1, num_covariates = 0, 
                              variable_label = NULL){
  if (!is(data, "data.frame") && !is(data, "matrix")){
    stop("data must be a data.frame or a matrix.")
  }
  if (ncol(data) < 2){
    stop("data must contain at lest two columns.")
  }
  if (num_mainvar > (ncol(data) - 1)){
    stop("Maximum number of variables is number of columns minus 1.")
  }
  if (num_mainvar < 0){
    stop("Number of main variables must be a positive numerical.")
  }
  if (num_covariates > (ncol(data) - num_mainvar - 1)){
    stop("Maximum number of covariates is number of columns minus number of main variables minus 1.")
  }
  if (num_covariates < 0){
    stop("Number of covariates must be a positive numerical.")
  }
  dependent <- colnames(data)[1]
  mainvariable <- colnames(data)[2:(1 + num_mainvar)]
  variables <- NULL
  covariates <- NULL
  if (ncol(data) > 2){
    if(ncol(data) > 1 + num_mainvar + num_covariates){
      variables <- colnames(data)[(2 + num_mainvar):(ncol(data) - num_covariates)]
    }else{}
    if(num_covariates){
      covariates <- colnames(data)[(ncol(data) - num_covariates + 1):ncol(data)]
    }
  }else{}
  variablemodel <- summary(lm(
    formula(paste(dependent, "~", paste(c(mainvariable, covariates), collapse = "+ "))),
    data = data))$adj.r.squared
  result <- c(variable = variablemodel)
  if(!is.null(variables)){
    variablesmodels <- sapply(variables, function(x) summary(lm(
      formula(paste(dependent, "~", paste(c(x, covariates), collapse = "+"))),
      data = data))$adj.r.squared)
    result <- c(result, variablesmodels)
  }
  if(!is.null(covariates)){
    covariatesmodel <- summary(lm(
      formula(paste(dependent, "~", paste(covariates, collapse = "+ "))),
      data = data))$adj.r.squared
    result <- c(result, covariates = covariatesmodel)
  }
  if (!is.null(variable_label)){
    names(result)[1] <- variable_label
  }
  return(result)
}
