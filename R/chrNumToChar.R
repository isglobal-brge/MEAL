#' Convert chr numbers to chr strings
#' 
#' Given a vector of number representing the chromosomes, convert them to string 
#' (e.g 1 to chr1). 23 is consider chrX, 24 is chrY, 25 is chrXY (probes shared between 
#' chromosomes X and Y) and 26 is chrMT.
#' 
#' @export
#' @param vector The vector with the chromosome numbers
#' @return A vector with the chromosomes in string format.
#' @examples
#' chromosomes <- c(1, 3, 4, 23, 15)
#' stringChrs <- chrNumToChar(chromosomes)
#' stringChrs

chrNumToChar <- function(vector){
  vector <- as.character(vector)
  vector[vector == "23"] <- "X"
  vector[vector == "24"] <- "Y"
  vector[vector == "25"] <- "XY"
  vector[vector == "26"] <- "MT"
  vector <- paste0("chr", vector)
}