#' @import stringr
#' @import DESeq2
#' @import dplyr
#' Performs differential analysis expression using DESeq2
#'
#' @param countData Data frame with matrix counts obtained from readTEcounts
#' @param metadata  Metadata read with read_metadata
#' @param useCtrlGenes Boolean for using only genes in estimating DESeq2 size factors (default: FALSE)
#'
#' @returns DESeq2 dds object
#' @export
#'
#' @examples
call_deseq2 <- function(countData, metadata, useCtrlGenes=FALSE) {
  
  # Check correct factors in Condition.
  expected_levels <- c("Control", "Treat")
  if (!metadata$Condition %in% expected_levels) {
    stop(paste0("Metadata file 'Condition' columns has values ",
                "different to 'Control' and 'Treat'"))
  }
  
  # Create coldata from metadata
  colData <- metadata[, -1] %>% 
           dplyr::mutate(condition = factor(Condition, levels = expected_levels))
  
  # DESeq2 analysis
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = count.df,
                                colData = colData,
                                design = ~ condition)
  
  # Use genes+TEs or only genes to estimate size factors.
  # We assume TE annotation includes ":"
  if(useCtrlGenes) {
    ctrlGenes <- !stringr::str_detect(rownames(dds), ":")
    dds <- DESeq2::estimateSizeFactors(dds, controlGenes=ctrlGenes)
  }
  
  # DESeq2 analysis
  dds <- DESeq2::DESeq(dds)
  dds
}