# Convert a data.frame to a model matrix
# 
# Given a data.frame, transform it to a model matrix. If a equation is not specified, 
# the model will contain all the variables summed.
# 
# @details \code{createModel} is used to generate a model matrix. Equation can 
# contain a formula compatible with \code{formula} function. Therefore, it is possible
# to include interactions in the matrix model. It should be noticed that all the 
# variables present in data will be included additively in the final model, although 
# not being included in the equation.
# 
# If names vector length is different than number of columns of the model, the
# original names will be preserved. 
#  
# @param data Data.frame with variables.
# @param equation Character string containing the formula to be used to create 
# the model.
# @param names Character vector with the names of the columns of the model matrix. 
# @return Model.matrix
createModel <- function(data, equation = NULL, names = NULL){
  if (is.null(dim(data))){ 
    stop("data is not a valid data.frame")
  }
  if (ncol(data) == 0 || nrow(data) == 0){
    stop("data is empty.")
  }
  if (is.null(equation)){
    model <- model.matrix(~ ., data = data)
  }else{
    equation <- paste(equation, "+ .")
    equation <- tryCatch(formula(equation), 
                         error = function(e) stop("Equation is an invalid formula."))
    model <- tryCatch(model.matrix(equation, data = data), 
                      error = function(e) stop("A model can't be created with this equation."))
  }
  if (!is.null(names)){
    if (length(names) == 0){
      warning("names variable is empty.")
    }else if (length(names) != ncol(model)){
      warning("Length of names is different than column names of model. Names won't be changed.")
    }else{
      colnames(model) <- names
    }
  }
  
  model
}