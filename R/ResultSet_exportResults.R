#' Exports results data.frames to csv files.
#' 
#' Exports results to csv files. If more than one variable is present, subfolders
#' with the name of the variable are created. For each variable, four files will 
#' be generated: probeResults.csv, dmrCateResults.csv, bumphunterResults.csv
#' and blockFinderResults.csv
#' 
#' @name exportResults
#' @aliases exportResults 
#' 
#' @export 
#' 
#' @param object \code{ResultSet}
#' @param dir Character with the path to export.
#' @param prefix Character with a prefix to be added to all file names.
#' @param fNames Names of the columns of \code{object} fData that will be added to 
#' the results data.frame.
#' @return Files are saved into the given folder.
#' @examples
#' if (require(minfiData)){
#' set <- ratioConvert(mapToGenome(MsetEx[1:10,]))
#' methyOneVar <- runPipeline(set, variable_names = "sex")
#' exportResults(methyOneVar)
#' }
exportResults <- function(object, dir = "./", prefix = NULL, 
                          fNames = c("chromosome", "start")) {
  stopifnot(is(object, "ResultSet"))
  if (substr(dir, nchar(dir), nchar(dir)) != "/"){
    dir <- paste0(dir, "/")
  }
  if (!file.exists(dir)){
    dir.create(dir)
  }
  if ("DiffMean" %in% names(object)) {
    temp <- getProbeResults(object, rid = "DiffMean", fNames = fNames)
    write.csv2(temp, file = paste0(dir, prefix, "DiffMeanResults.csv"))
  }
  if ("DiffVar" %in% names(object)) {
    temp <- getProbeResults(object, rid = "DiffVar", fNames = fNames)
    write.csv2(temp, file = paste0(dir, prefix, "DiffVarResults.csv"))
  }
  if ("dmrcate" %in% names(object)) {
    temp <- getAssociation(object, rid = "dmrcate")
    write.csv2(temp, file = paste0(dir, prefix, "dmrCateResults.csv"))
  }
  if ("bumphunter" %in% names(object)) {
    temp <- getAssociation(object, rid = "bumphunter")
    write.csv2(temp, file = paste0(dir, prefix, "bumphunterResults.csv"))
  }
  if ("blockFinder" %in% names(object)) {
    temp <- getAssociation(object, rid = "blockFinder")
    write.csv2(temp, file = paste0(dir, prefix, "blockFinderResults.csv")) 
  }
}