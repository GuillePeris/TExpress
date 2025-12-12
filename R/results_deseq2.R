#' Title
#'
#' @param dds Objecto from call_deseq2 function
#' @param shrinklog2FC Boolean for applying log2FC shrinking in DESEq2 (default: FALSE)
#'
#' @returns Data frame with DESeq2 results
#' @export
#'
#' @examples
results_deseq2 <- function(dds, shrinklog2FC=FALSE) {
  
  # Shrink log2FC if requires
  if(shrinklog2FC) {
    res <- DESeq2::lfcShrink(dds, coef=paste0("condition_", expected_levels[2], 
                                              "_vs_", expected_levels[1]),  
                             type="apeglm")
  } else {
    res <- DESeq2::results(dds)
  } 
  
  # Remove NA rows and sort by padj
  res_complete <- res[complete.cases(res), ]
  res_sorted <- res_complete[order(res_complete$padj), ]
  
  as.data.frame(res)
}