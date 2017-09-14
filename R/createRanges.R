#' Create \code{GenomicRanges} from data.frame
#' 
#' @export createRanges
#' @return Deprecated
#' @examples
#' createRanges()
createRanges <- function(){
  .Deprecated("makeGRangesFromDataFrame", 
              msg = "createRanges is deprecated. use makeGRangesFromDataFrame from GRanges package.")
}
