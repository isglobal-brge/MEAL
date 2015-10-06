setValidity("MethylationSet", function(object){
  msg <- validMsg(NULL, assayDataValidMembers(assayData(object), "meth"))
  if (is.null(msg)){
    TRUE
  } else{
    msg
  }
})