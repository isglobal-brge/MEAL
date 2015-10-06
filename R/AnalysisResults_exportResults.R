#' @describeIn AnalysisResults  Exports results data.frames to csv files.
#' @aliases AnalysisResults-methods
#' @param dir Character with the path to export.
#' @param prefix Character with a prefix to be added to all file names.
#' @param vars Character vector with the names of the variables to be exported. Note: 
#' names should be that of the model. 
setMethod(
  f = "exportResults",
  signature = "AnalysisResults",
  definition = function(object, dir = "./", prefix = NULL, vars = modelVariables(object)) {
    if (substr(dir, nchar(dir), nchar(dir)) != "/"){
      dir <- paste0(dir, "/")
    }
    if (sum(vars %in% modelVariables(object)) == 0){
      stop("Vars introduced are not present in the results object.")
    }
    if (!file.exists(dir)){
      dir.create(dir)
    }
    if (length(vars) == 1){
      temp <- probeResults(object)[[1]]
      write.csv2(temp, file = paste0(dir, prefix, "probeResults.csv"))
      temp <- dmrCate(object)[[1]]
      write.csv2(temp, file = paste0(dir, prefix, "dmrCateResults.csv"))
      temp <- bumps(object)[[1]]
      write.csv2(temp, file = paste0(dir, prefix, "bumphunterResults.csv"))
      temp <- blocks(object)[[1]]
      write.csv2(temp, file = paste0(dir, prefix, "blockFinderResults.csv")) 
      
    } else {
      for (var in vars){
        tempdir <- paste0(dir, var, "/")
        if (!file.exists(tempdir)){
          dir.create(tempdir)
        }
        temp <- probeResults(object)[[var]]
        write.csv2(temp, file = paste0(tempdir, prefix, "probeResults.csv"))
        temp <- dmrCate(object)[[var]]
        write.csv2(temp, file = paste0(tempdir, prefix, "dmrCateResults.csv"))
        temp <- bumps(object)[[var]]
        write.csv2(temp, file = paste0(tempdir, prefix, "bumphunterResults.csv"))
        temp <- blocks(object)[[var]]
        write.csv2(temp, file = paste0(tempdir, prefix, "blockFinderResults.csv")) 
      }    
    }
  }
)