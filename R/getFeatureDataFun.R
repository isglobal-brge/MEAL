# Return the function to get feature data
#
# Returns the function required to get feature data depending on the incoming 
# object. This function ensures the standarization of \code{resultSet} fData.
# @param set \code{eSet}, \code{SummarizedExperiment} or \code{RangedSummarizedExperiment}
getFeatureDataFun <- function(set){
  if (is(set, "eSet")){
    fFun <- function(set) {
      df <- Biobase::fData(set)
      ## Change position by start to harmonize fDatas
      if ("position" %in% colnames(df)){
        colnames(df)[colnames(df) == "position"] <- "start"
        df$end <- df$start
      }
      if ("chr" %in% colnames(df)){
        colnames(df)[colnames(df) == "chr"] <- "chromosome"
      } 
      df
    }
  } else if (is(set, "RangedSummarizedExperiment")){
    fFun <- function(set) { 
      df <- as.data.frame(SummarizedExperiment::rowRanges(set), 
                          stringsAsFactors = FALSE)
      colnames(df)[1] <- "chromosome"
      df[, 1] <- as.character(df[, 1])
      df
    } 
  } else if (is(set, "SummarizedExperiment")){
    fFun <- SummarizedExperiment::rowData
  } else if (is.matrix(set)){
    fFun <- function(set) data.frame(matrix(vector(), nrow(set), 0))
  }
  return(fFun)
}