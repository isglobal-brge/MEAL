# Filter the data.frame obtained from probe analysis
# 
# @param results Data.frame with the results of probe analysis
# @param range \code{GenomicRanges} with the desired range. 
# @param position Character with the name of the column containing the positions
# @param chr Character with the name of the column containing the chromosome
# @return Data.frame with the results of the probes of the range

filterResults <- function(results, range, position = "position", chr = "chromosome"){
  if(ncol(results) == 0 | nrow(results) == 0){
    stop("results is empty")
  }
    
  if (!is(range, "GenomicRanges")){
    stop("range must be a GenomicRanges object.")
  }
  if (length(range) == 0){
    stop("range is empty.")
  }
  
  if (length(range) != 1){
    stop("range must be a GenomicRanges with only one range.")
  }
      
  maskResults <- as.vector(results[ , chr] == as.character(S4Vectors::runValue(GenomicRanges::seqnames(range)))) & 
    results[ , position] >= start(range) &  results[ , position] <= end(range)
  results <- results[maskResults, ]
  if (nrow(results) == 0){
    warning("There are no cpgs in the range. An empty data.frame will be returned.")
  }
  return(results)
}