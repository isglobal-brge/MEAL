#' Filter a \code{MethylationSet}, an \code{ExpressionSet} or a \code{SnpSet}
#' 
#' @name filterSet
#' @rdname filterSet-methods
#' @export 
#' @param set \code{MethylationSet}, \code{ExpressionSet} or a \code{SnpSet}
#' @param range \code{GenomicRanges} with the desired range. 
#' @return \code{MethylationSet}, \code{ExpressionSet} or a \code{SnpSet} with only 
#' the features of the range.
#' @examples
#' if (require(minfiData) & require(GenomicRanges)){
#' range <- GRanges(seqnames=Rle("chrY"), 
#' ranges = IRanges(3000000, end=12300000))
#' set <- prepareMethylationSet(MsetEx[1:100, ], pData(MsetEx))
#' set
#' filteredset <- filterSet(set, range)
#' filteredset
#' }
filterSet <- function(set, range){
  if (!(is(set, "MethylationSet") | is(set, "ExpressionSet") | is(set, "SnpSet"))){
    stop("set must be a MethylationSet, an ExpressionSet or a SnpSet")
  }
  if(ncol(set) == 0 | nrow(set) == 0){
    stop("set is empty")
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
  
  if (is(set, "MethylationSet") | is(set, "SnpSet")){
    positions <- fData(set)[ , c("chromosome", "position")]
    positions <- positions[apply(!is.na(positions), 1, all), ]
    mask <- as.vector(positions[ , 1] == as.character(S4Vectors::runValue(GenomicRanges::seqnames(range)))) & 
      positions[ , 2] >= start(range) & positions[ , 2] <= end(range)
  }else{
    if (!all(c("chromosome", "start", "end") %in% fvarLabels(set))){
      stop("Expression needs to contain feature data with chromosome, start and end columns.")
    }
    positions <- fData(set)[ , c("chromosome", "start", "end")]
    positions <- positions[apply(!is.na(positions), 1, all), ]
    mask <- as.vector(positions[ , 1] == as.character(S4Vectors::runValue(GenomicRanges::seqnames(range)))) & 
      ((positions[ , 2] <= start(range) & positions[ , 3] >= start(range)) | 
      (positions[ , 2] >= start(range) & positions[ , 3] <= end(range)))
  }

  set <- set[mask, ]
  
  if (sum(mask) == 0){
    warning("There are no features in the range. An empty set will be returned.")
  }
  
  set
}
