#' @title Compute differential expression analysis for transposable elements 
#'
#' @param metafile 
#' @param folder 
#' @param useCtrlGenes 
#' @param shrinklog2FC 
#'
#' @returns Table of results
#' @export
#'
TE_DEA <- function(metafile, folder, useCtrlGenes=FALSE, 
                         shrinklog2FC=FALSE) {
  
  # Read metadata
  metadata <- read_metadata(metafile)
  
  # Read counts into a count data frame
  countData <- readTEexpress(metadata, folder)
  
  # DESeq2 analysis
  dds <- call_deseq2(countData, metadata, useCtrlGenes)
  res <- results_deseq2(dds, shrinklog2FC)
  
  # Save normalized counts
  norm.counts <- norm_counts(dds)
  
  res
}