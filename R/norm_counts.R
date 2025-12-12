#' @import DESeq2
#' @title Get normalized counts from DESeq2 object
#'
#' @param dds DEseq2 object
#'
#' @returns Data frame with normalized counts by samples
#' @export
#'
norm_counts <- function(dds) {
  DESeq2::counts(dds, normalized=TRUE)
}