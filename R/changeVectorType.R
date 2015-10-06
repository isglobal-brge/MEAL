# Changes vector's class
# 
# vector: vector to be transformed
# type: string with the destiny type. Available types are: categorical,
# continuous, dominant, recessive and additive. If NA, no change is performed. 
# return: The same vector with the type changed. If a genetic model is used, values 
# will be changed.  

changeVectorType <- function(vector, type){
  if (is.na(type)){
    vector <- vector
  }else if (type == "categorical"){
    vector <- as.factor(vector)
  }
  else if (type == "continuous"){
    vector <- as.numeric(vector)
  }
  else if (type == "dominant"){
    vector <- SNPassoc::dominant(vector)
  }
  else if (type == "recessive"){
    vector <- SNPassoc::recessive(vector)
  }
  else if (type == "additive"){
    vector <- SNPassoc::additive(vector)
  }else{
    warning("Invalid type. No change will be performed")    
  }
  return(vector)
}