#' Create \code{GenomicRanges} from data.frame
#' 
#' Convert a data.frame with chromosomes in the first column, starting positions in
#' the second one and ending position in the third one to \code{GenomicRanges}. Names 
#' of the data.frame are preserved in the output \code{GenomicRanges}.
#' 
#' @param ranges Data.frame or matrix
#' @return \code{GenomicRanges}
#' 
#' @export createRanges
#' @examples
#' dfranges <- data.frame(chr = c("chr1", "chr2", "chr1"), start = c(1290, 1250, 4758),
#' end = c(64389, 632409, 16430), stringsAsFactors = FALSE)
#' names(dfranges) <- c("range1", "range2", "range3")
#' ranges <- createRanges(dfranges)
#' ranges
createRanges <- function(ranges){
  if (is(ranges, "data.frame") | is(ranges, "matrix")){
    Range <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(ranges[ , 1]), 
                     ranges = IRanges::IRanges(ranges[ , 2], end = ranges[ , 3]))
    names(Range) <- rownames(ranges)
    Range
  }
}