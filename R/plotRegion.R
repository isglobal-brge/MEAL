#' Plot results in a genomic region
#' 
#' Plot the results from the different analyses of a \code{ResultSet} in a specific
#' genomic region. It can plot all the results from \code{runPipeline}. 
#' 
#' @details This plot allows to have a quick summary of the methylation or gene
#' expression analyses in a given region. If we use a \code{ResultSet} obtained 
#' from methylation data, transcripts annotation is obtained from archive. If we
#' use a \code{ResultSet} obtained from gene expression data, transcripts annotation
#' is taken from fData. 
#' 
#' This plot can be used to plot the results of one dataset (methylation or gene
#' expression) or to represent the association between methylation and gene 
#' expression data. If only one dataset is used, the p-values and the coefficients 
#' of DiffMean and DiffVar analyses are plotted. If we pass two \code{ResultSet}s, 
#' \code{rset} should contain methylation results and a \code{rset2} the gene expression 
#' results.
#' 
#' @param rset \code{ResultSet}
#' @param range \code{GenomicRanges} with the region coordinates
#' @param results Character with the analyses that will be included in the plot.
#' By default, all analyses available are included. 
#' @param genome String with the genome used to retrieve transcripts annotation:
#' hg19, hg38, mm10. (Default: "hg19")
#' @param tPV Threshold for P-Value
#' @param rset2 Additional \code{ResultSet}
#' @return Regional plot
#' @export
plotRegion <- function(rset, range, results = names(rset),
                       genome = "hg19", rset2, tPV = 5,
                       fNames = c("chromosome", "start", "end"),
                       fNames2 = c("chromosome", "start", "end")){

  ### Afegir linea per separar blocs. 

  ## Create Basic tracks
  ## Add basic tracks tracks superiores (Idiogram más trancripción)
  basetracks <- list(Gviz::IdeogramTrack(genome = genome, 
                                         chromosome = as.character(seqnames(range))), 
                     Gviz::GenomeAxisTrack())
  
  rset@results <- rset@results[results]
  
  if (missing(rset2)){
    
    
    
    if ("end" %in% colnames(fData(rset)[[1]])){
      endf <- "end"
    } else {
      endf <- "start"
    }
    annot <- GenomicRanges::makeGRangesFromDataFrame(fData(rset)[[1]], 
                                                     seqnames.field = c("chromosome", "seqnames"), 
                                                     start.field = c("position", "start"), 
                                                     end.field = endf)
    if (width(annot[1]) == 1){
      data(dmrcatedata, envir = environment(), package = "DMRcatedata")
      switch(genome, hg19 = {
        tx = tx.hg19
      }, hg38 = {
        tx = tx.hg38
      }, mm10 = {
        tx = tx.mm10
      })
      annotTrack <- c(
        Gviz::GeneRegionTrack(subsetByOverlaps(tx, range), name = "Transcripts", 
                              symbol = subsetByOverlaps(tx, range)$gene_name, 
                              fill = "lightblue", 
                              gene = subsetByOverlaps(tx, range)$gene_name,
                              showId = TRUE, geneSymbol = TRUE, cex.title = 0.6,
                              shape = "arrow", transcriptAnnotation = "symbol",
                              collapseTranscripts = TRUE, rotation.title = 0), 
        Gviz::GeneRegionTrack(subsetByOverlaps(annot, range), name = "CpGs", 
                              fill = "black", stacking = "dense",rotation.title = 0)) 
    } else {
      annotTrack <- Gviz::GeneRegionTrack(subsetByOverlaps(annot, range), name = "Transcripts", 
                                          symbol = names(subsetByOverlaps(annot, range)), 
                                          fill = "lightblue", 
                                          gene = names(subsetByOverlaps(annot, range)),
                                          showId = TRUE, geneSymbol = TRUE,
                                          shape = "arrow", transcriptAnnotation = "symbol",
                                          collapseTranscripts = TRUE) 
    }
  } else {
    rset2@results <- rset2@results[names(rset2) %in% results]
    
    annot <- GenomicRanges::makeGRangesFromDataFrame(fData(rset)[[1]], 
                                                     seqnames.field = c("chromosome", "seqnames"), 
                                                     start.field = c("position", "start"), 
                                                     end.field = "start")
    annot2 <- GenomicRanges::makeGRangesFromDataFrame(fData(rset2)[[1]], 
                                                      seqnames.field = c("chromosome", "seqnames"), 
                                                      start.field = c("position", "start"), 
                                                      end.field = "end")
    annotTrack <- c(
      Gviz::GeneRegionTrack(subsetByOverlaps(annot2, range), name = "Transcripts", 
                            symbol = names(subsetByOverlaps(annot2, range)), 
                            fill = "lightblue", 
                            gene = names(subsetByOverlaps(annot2, range)),
                            showId = TRUE, geneSymbol = TRUE,
                            shape = "arrow", transcriptAnnotation = "symbol",
                            collapseTranscripts = TRUE), 
      Gviz::GeneRegionTrack(subsetByOverlaps(annot, range), name = "CpGs", 
                            fill = "black", stacking = "dense",rotation.title = 0)) 
  }
  resultTrack <- list()
  
  if (missing(rset2)){
    
    if ("DiffMean" %in% names(rset)){
      diffMean <- getProbeResults(rset, rid = "DiffMean", fNames = c("chromosome", "start", "end"))
      diffMeanGR <- GenomicRanges::makeGRangesFromDataFrame(
        diffMean, seqnames.field = c("chromosome", "seqnames"), 
        start.field = "start", 
        end.field = "end", keep.extra.columns = TRUE)
      diffMeanGR <- subsetByOverlaps(diffMeanGR, range)
      resultTrack <- c(resultTrack, Gviz::DataTrack(range = diffMeanGR, 
                                                    data = diffMeanGR$logFC, 
                                                    genome = genome, name = "DiffMean", type = "h", 
                                                    baseline = 0),
                       Gviz::DataTrack(range = diffMeanGR, data = -log10(diffMeanGR$P.Value), 
                                       genome = genome, name = "DiffMean -log10 p-value", 
                                       type = "p", col = "darkblue",
                                       baseline = tPV))
    }
    
    
    if ("DiffVar" %in% names(rset)){
      DiffVar <- getProbeResults(rset, rid = "DiffVar", 
                                 fNames = c("chromosome", "start", "end"))
      DiffVarGR <- GenomicRanges::makeGRangesFromDataFrame(
        DiffVar, seqnames.field = c("chromosome", "seqnames"), 
        start.field = "start", 
        end.field = "end", keep.extra.columns = TRUE)
      DiffVarGR <- subsetByOverlaps(DiffVarGR, range)
      resultTrack <- c(resultTrack, Gviz::DataTrack(
        range = DiffVarGR, data = DiffVarGR$DiffLevene, 
        genome = genome, name = "DiffVar", type = "h", baseline = 0, col = "darkgreen"),
        Gviz::DataTrack(range = DiffVarGR, data = -log10(DiffVarGR$P.Value), 
                        baseline = tPV, genome = genome, name = "DiffVar -log10 p-value", 
                        type = "p", col = "darkolivegreen"))
    }
  } else {
    if ("DiffMean" %in% names(rset)){
      diffMean <- getProbeResults(rset, rid = "DiffMean", fNames = c("chromosome", "start", "end"))
      diffMeanGR <- GenomicRanges::makeGRangesFromDataFrame(
        diffMean, seqnames.field = c("chromosome", "seqnames"), 
        start.field = "start", 
        end.field = "end", keep.extra.columns = TRUE)
      diffMeanGR <- subsetByOverlaps(diffMeanGR, range)
      resultTrack <- c(resultTrack, Gviz::DataTrack(range = diffMeanGR, 
                                                    data = diffMeanGR$logFC, 
                                                    genome = genome, 
                                                    name = "DiffMean Res1", type = "h", 
                                                    baseline = 0, 
                                                    col = "black"))
    }
    if ("DiffMean" %in% names(rset2)){
      diffMean <- getProbeResults(rset2, rid = "DiffMean", fNames = fNames2)
      diffMeanGR <- GenomicRanges::makeGRangesFromDataFrame(
        diffMean, seqnames.field = c("chromosome", "seqnames"), 
        start.field = fNames2[2], 
        end.field = fNames2[3], keep.extra.columns = TRUE)
      diffMeanGR <- subsetByOverlaps(diffMeanGR, range)
      resultTrack <- c(resultTrack, Gviz::DataTrack(range = diffMeanGR, 
                                                    data = diffMeanGR$logFC, 
                                                    genome = genome, name = "DiffMean Res2", type = "h", 
                                                    baseline = 0))
    }
    
    if ("DiffVar" %in% names(rset)){
      DiffVar <- getProbeResults(rset, rid = "DiffVar", fNames = fNames)
      DiffVarGR <- GenomicRanges::makeGRangesFromDataFrame(
        DiffVar, seqnames.field = c("chromosome", "seqnames"), 
        start.field = fNames[2], 
        end.field = fNames[3], keep.extra.columns = TRUE)
      DiffVarGR <- subsetByOverlaps(DiffVarGR, range)
      resultTrack <- c(resultTrack, Gviz::DataTrack(
        range = DiffVarGR, data = DiffVarGR$DiffLevene, col = "black",
        genome = genome, name = "DiffVar Res1", type = "h", baseline = 0))
    }
    if ("DiffVar" %in% names(rset2)){
      DiffVar <- getProbeResults(rset2, rid = "DiffVar", fNames = fNames)
      DiffVarGR <- GenomicRanges::makeGRangesFromDataFrame(
        DiffVar, seqnames.field = c("chromosome", "seqnames"), 
        start.field = fNames[2], 
        end.field = fNames[3], keep.extra.columns = TRUE)
      DiffVarGR <- subsetByOverlaps(DiffVarGR, range)
      resultTrack <- c(resultTrack, Gviz::DataTrack(
        range = DiffVarGR, data = DiffVarGR$DiffLevene, 
        genome = genome, name = "DiffVar Res2", type = "h", baseline = 0))
    }
  }
  
  dmrTracks <- list()
  # Create tracks from Region Analysis
  if ("bumphunter" %in% names(rset)){
    bumphunter <- MultiDataSet::getAssociation(rset, rid = "bumphunter")
    if (nrow(bumphunter) != 0){
      bumphunterGR <- GenomicRanges::makeGRangesFromDataFrame(
        bumphunter, seqnames.field = "chr", start.field = "start",
        end.field = "end", keep.extra.columns = TRUE)
      bumphunterGR$id <- names(bumphunterGR)
      if (!"p.value" %in% colnames(bumphunter)){
        bumphunterGR$col <- "null"
      } else {
        bumphunterGR$col <- ifelse(bumphunter$p.value < 0.05, "sig", "no.sig")
      }
      bumphunterGR <- subsetByOverlaps(bumphunterGR, range)
      bumptrack <- Gviz::AnnotationTrack(bumphunterGR,
                                         name = "Bumphunter", shape = "box",
                                         group = bumphunterGR$id, 
                                         feature = bumphunterGR$col,
                                         rotation.title = 0, 
                                         cex.title = 0.6)
      dmrTracks <- c(dmrTracks, bumptrack)
    }
  }
  
  if ("blockFinder" %in% names(rset)){
    blockFinder <- MultiDataSet::getAssociation(rset, rid = "blockFinder")
    if (nrow(blockFinder) != 0){
      blockFinderGR <- GenomicRanges::makeGRangesFromDataFrame(
        blockFinder, seqnames.field = "chr", start.field = "start", end.field = "end",
        keep.extra.columns = FALSE)
      blockFinderGR$id <- names(blockFinderGR)
      if (!"p.value" %in% colnames(blockFinder)){
        blockFinderGR$col <- "null"
      } else {
        blockFinderGR$col <- ifelse(blockFinder$p.value < 0.05, "sig", "no.sig")
      }
      blockFinderGR <- subsetByOverlaps(blockFinderGR, range)
      blockTrack <- Gviz::AnnotationTrack(blockFinderGR,
                                          name = "blockFinder", shape = "box",
                                          group = blockFinderGR$id, 
                                          feature = blockFinderGR$col,
                                          rotation.title = 0, 
                                          cex.title = 0.6)
      dmrTracks <- c(dmrTracks, blockTrack)
    }
  }
  
  if ("dmrcate" %in% names(rset)){
    dmrcate <- MultiDataSet::getAssociation(rset, rid = "dmrcate")
    if (nrow(dmrcate) != 0){
      colnames(dmrcate)[colnames(dmrcate) == "Stouffer"] <- "p.val"
      dmrcateGR <- GRanges(dmrcate$coord)
      mcols(dmrcateGR) <- dmrcate[, -1]
      dmrcateGR$id <- 1:nrow(dmrcate)
      dmrcateGR$col <- ifelse(dmrcateGR$p.val < 0.05, "sig", "no.sig")
      dmrcateGR <- subsetByOverlaps(dmrcateGR, range)
      dmrcateTrack <- Gviz::AnnotationTrack(dmrcateGR,
                                            name = "DMRcate", shape = "box", genome = genome,
                                            group = dmrcateGR$id,  
                                            cex.title = 0.6, 
                                            feature = dmrcateGR$col,
                                            rotation.title = 0, showTitle = TRUE)
      dmrTracks <- c(dmrTracks, dmrcateTrack)
    }
  }
  Gviz::plotTracks(c(basetracks, annotTrack, dmrTracks, resultTrack), 
                   from = start(range), to = end(range), title.width = 1.5,
                   groupAnnotation = "group", just.group = "above", 
                   null = "white", sig = "forestgreen", no.sig = "grey")
}














