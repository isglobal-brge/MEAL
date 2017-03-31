#' Generating a \code{MethylationSet}
#' 
#' This function creates a \code{MethylationSet} using from a matrix of beta values and
#' a data.frame of phenotypes.
#' 
#' @details prepareMethylationSet is a useful wrapper to create \code{MethylationSet}.
#' Rigth now, prepareMethylationSet supports two entry points: a \code{minfi}
#' object and a matrix of betas.
#' 
#' Phenotypes are compulsory and can be supplied as data.frame or AnnotatedDataFrame.
#' 
#' By default, annotation is taken from \code{minfi} package and 
#' \code{IlluminaHumanMethylation450kanno.ilmn12.hg19} package is used, being the 
#' default arguments adapted to use this annotation. To use this annotation, 
#' \code{IlluminaHumanMethylation450kanno.ilmn12.hg19} must be installed and 
#' methylation sites must be named like in Illumina 450k chip. Use of this annotation
#' ensures correct results in all the analysis.  
#' 
#' If custom annotation is desired, there are two compulsory features: chromosomes and
#' positions. Chromosomes should be supplied in the character form (e.g. chr1). Two
#' additional features will be used during the presentation of results but not during
#' the analyses: genes and group. Genes are the gene names of the genes around the cpg
#' site and group defines the groups of the genes. Both columns will appear in the 
#' results but they are not used through the workflow. It should be noticed that \code{BlockFinder}
#' only supports \code{minfi} annotation, so it is not advised to be used with custom
#' annotation.  
#' 
#' @export prepareMethylationSet
#' @param matrix Data.frame or a matrix with samples on the columns and cpgs
#' on the rows. A \code{minfi} object can be used to. 
#' @param phenotypes Data.frame or vector with the phenotypic features of the samples.
#' Samples will be in the rows and variables in the columns. If matrix is a \code{minfi}
#' object, phenotypes can be taken from it. 
#' @param annotation Character with the name of the annotation package or data.frame
#' or AnnotationDataFrame with the annotation.
#' @param chromosome Character with the column containing chromosome name in the
#' annotation data. 
#' @param position chromosome Character with the column containing position coordinate in the
#' annotation data. 
#' @param genes Character with the column containing gene names related to the
#' methylation site in the annotation data. (Optional)
#' @param group Character with the column containing the position of the probe
#' related to the gene named in gene column. (Optional)
#' @param filterNA_threshold Numeric with the maximum percentage of NA allowed for 
#' each of the probes. If 1, there will be no filtering, if 0 all probes containing at 
#' least a NA will be filtered. 
#' @param verbose Logical value. If TRUE, it writes out some messages indicating progress. 
#' If FALSE nothing should be printed.
#' @return \code{MethylationSet} with phenotypes and annotation.
#' @examples
#' if (require(minfiData)){
#'  betas <- getBeta(MsetEx)[1:1000, ]
#'  pheno <- pData(MsetEx)
#'  set <- prepareMethylationSet(betas, pheno)
#'  }
prepareMethylationSet <- function(matrix, phenotypes, 
                                  annotation = "IlluminaHumanMethylation450kanno.ilmn12.hg19", 
                                  chromosome = "chr", position = "pos", 
                                  genes = "UCSC_RefGene_Name",
                                  group = "UCSC_RefGene_Group",
                                  filterNA_threshold = 0.05,
                                  verbose = FALSE){
    if (verbose){
        message("Creating the object...")
    }
    
    if (length(matrix) == 0){
        stop("Matrix is empty.")
    }
    
    #If matrix is a data.frame or a matrix, assure that only contains numbers, 
    #filter NAs and convert to RatioSet
    if (is(matrix, "data.frame")){
        matrix <- data.matrix(matrix)
    }
    if (is(matrix, "matrix")){
        matrix <- matrix[rowMeans(is.na(matrix)) <= filterNA_threshold, , drop = FALSE]
    } else if (class(matrix) %in% c("GenomicMethylSet", "MethylSet", "RatioSet", "GenomicRatioSet")){
        if (missing(phenotypes)){
            phenotypes <- colData(matrix)
        }
        matrix <- minfi::getBeta(matrix)
    }else{
        stop("matrix is not a minfi class nor a data.frame nor a matrix")
    }
    
    if (is.null(rownames(matrix))){
        stop("Rownames of matrix must contain probe names.")
    }
    if (is.null(colnames(matrix))){
        stop("Colnames of matrix must contain sample names.")
    }
    
    if (length(phenotypes) == 0){
        stop("Phenotypes is empty.")
    }
    if (is.null(dim(phenotypes))){
      stop("Phenotypes must be a data.frame, a DataFrame or an AnnotatedDataFrame.")
    }

    if (is(annotation, "character")){
        if (annotation == "IlluminaHumanMethylation450kanno.ilmn12.hg19"){
            annoChar <- annotation
            annotation <- data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19::Locations, 
                                     IlluminaHumanMethylation450kanno.ilmn12.hg19::Other)
        } else if (annotation == "IlluminaHumanMethylationEPICanno.ilm10b2.hg19"){
            annoChar <- annotation
            annotation <- data.frame(IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Locations, 
                                     IlluminaHumanMethylationEPICanno.ilm10b2.hg19::Other)
        }
        else{
            stop(sprintf("Annotation package %s can not be found.", annotation))
        }
        
    }else{
        annoChar <- "custom"
    }
    
    if (!chromosome %in% colnames(annotation)){
        stop("Annotation data.frame must contain a column name equal to \"chromosome\" argument.")
    }else{
        colnames(annotation)[colnames(annotation) == chromosome] <- "chromosome"
    }
    
    if (!position %in% colnames(annotation)){
        stop("Annotation data.frame must contain a column name equal to \"position\" argument.")
    }else{
        colnames(annotation)[colnames(annotation) == position] <- "position"
    }
    
    if (!genes %in% colnames(annotation)){
        warning("Annotation data.frame does not contain a column name equal to \"gene\" parameter. Results may be incomplete.")
        if (is(annotation, "AnnotatedDataFrame")){
            pData(annotation)[ , "genes"] <- rep("", nrow(annotation))
        } else{
            annotation[ , "genes"] <- rep("", nrow(annotation))
        }
    }else{
        colnames(annotation)[colnames(annotation) == genes] <- "genes"
    }
    
    if (!group %in% colnames(annotation)){
        warning("Annotation data.frame does not contain a column name equal to \"group\" parameter. Results may be incomplete.")
        if (is(annotation, "AnnotatedDataFrame")){
            pData(annotation)[ , "group"] <- rep("", nrow(annotation))
        } else{
            annotation[ , "group"] <- rep("", nrow(annotation))
        }
    }else{
        colnames(annotation)[colnames(annotation) == group] <- "group"
    }
    
    set <- methylationSet(betas = matrix, phenotypes = phenotypes[colnames(matrix), , drop = FALSE],
                          annotationDataFrame = annotation[rownames(matrix), , drop = FALSE], 
                          annoString = annoChar)
    if (verbose){
        message("Checking the object...")
    }
    set <- checkSamples(set)
    set <- checkProbes(set)
    validObject(set)
    set
}