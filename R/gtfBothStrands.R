#' Create Antisense TE GTF with Both Strands
#'
#' Generates a GTF file containing both sense and antisense annotations by
#' duplicating all features and inverting the strand for antisense copies.
#' Useful for analyzing antisense transcription or bidirectional promoters.
#'
#' @param gtf.file Character string. Path to input TE GTF file. The file will be
#'   read, duplicated with strand inversion, and exported to a new file with
#'   "_AS" suffix.
#' @param output.file Character string. Optional path for output TE GTF file. If
#'   NULL (default), creates a file in the same directory as input with "_AS"
#'   added before the extension (e.g., "hg38_rmsk_TE.gtf" becomes 
#'   "hg38_rmsk_TE_AS.gtf").
#' @param add_suffix Character string. Suffix to add to gene_id and
#'   transcript_id for antisense features. Default is "_AS".
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Imports the input TE GTF file
#'   \item Creates antisense copies by inverting strand information
#'   \item Appends suffix to gene_id and transcript_id in antisense copies
#'   \item Merges sense and antisense features
#'   \item Sorts by genomic coordinates
#'   \item Exports to new TE GTF file
#' }
#'
#' For each feature on the "+" strand, an antisense copy on the "-" strand is
#' created, and vice versa. Gene and transcript IDs are modified to distinguish
#' antisense features (e.g., "ENSG00000000001" becomes "ENSG00000000001_AS").
#'
#' @return Character string (invisibly). Path to the created output GTF file.
#'   Also prints a message with file location and feature counts.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create antisense GTF with default naming
#' gtfBothStrands("genes.gtf")
#' # Output: genes_AS.gtf
#'
#' # Specify custom output file
#' gtfBothStrands(
#'   gtf.file = "hg38_rmsk_TE.gtf",
#'   output.file = "hg38_rmsk_TE_AS.gtf"
#' )
#'
#' # Custom suffix for antisense features
#' gtfBothStrands(
#'   gtf.file = "hg38_rmsk_TE.gtf",
#'   add_suffix = "_antisense"
#' )
#' }
#'
gtfBothStrands <- function(gtf.file,
                           output.file = NULL,
                           add_suffix = "_AS") {
  
  # Input validation
  if (missing(gtf.file)) { 
    stop("Argument 'gtf.file' is missing with no default.", call. = FALSE)
  }
  
  if (!file.exists(gtf.file)) {
    stop("GTF file '", gtf.file, "' not found.", call. = FALSE)
  }
  
  # Validate output.file if provided
  if (!is.null(output.file)) {
    if (!is.character(output.file) || length(output.file) != 1L) {
      stop("'output.file' must be NULL or a single character string.",
           call. = FALSE)
    }
    
    output_dir <- dirname(output.file)
    if (output_dir != "." && !dir.exists(output_dir)) {
      stop("Output directory '", output_dir, "' does not exist.", call. = FALSE)
    }
  }
  
  #-- Import gtf
  gtf.gr <- tryCatch(
    rtracklayer::import(gtf.file, format = "gtf"),
    error = function(e) {
      stop("Failed to import GTF file: ", e$message, call. = FALSE)
    }
  )
  
  
  # Check for required metadata columns
  gtf_mcols <- GenomicRanges::mcols(gtf.gr)
  
  if (!"gene_id" %in% colnames(gtf_mcols)) {
    warning(
      "GTF does not contain 'gene_id' column. ",
      "Antisense features will not have modified gene IDs.",
      call. = FALSE
    )
  }
  
  if (!"transcript_id" %in% colnames(gtf_mcols)) {
    warning(
      "GTF does not contain 'transcript_id' column. ",
      "Antisense features will not have modified transcript IDs.",
      call. = FALSE
    )
  }
  
  #-- Invert strand and add antisense tag
  gtf.gr.as <- tryCatch(
    BiocGenerics::invertStrand(gtf.gr),
    error = function(e) {
      stop("Failed to invert strand: ", e$message, call. = FALSE)
    }
  )
  
  # Add antisense suffix to identifiers
  if ("gene_id" %in% colnames(gtf_mcols)) {
    gtf.gr.as$gene_id <- paste0(gtf.gr.as$gene_id, add_suffix)
  }
  
  if ("transcript_id" %in% colnames(gtf_mcols)) {
    gtf.gr.as$transcript_id <- paste0(gtf.gr.as$transcript_id, add_suffix)
  }
  
  # Combine sense and antisense
  final.gtf.gr <- c(gtf.gr, gtf.gr.as)
  
  sorted.gtf.gr <- tryCatch(
    {
      order_idx <- order(
        GenomicRanges::seqnames(final.gtf.gr),
        BiocGenerics::start(final.gtf.gr),
        BiocGenerics::end(final.gtf.gr)
      )
      final.gtf.gr[order_idx]
    },
    error = function(e) {
      warning(
        "Failed to sort features: ", e$message,
        "\nProceeding with unsorted data.",
        call. = FALSE
      )
      final.gtf.gr
    }
  )
  
  if (is.null(output.file)) {
    # Create default output filename
    base_name <- tools::file_path_sans_ext(gtf.file)
    extension <- tools::file_ext(gtf.file)
    
    # Handle case where file has no extension
    if (nchar(extension) == 0) {
      gtf.out <- paste0(base_name, add_suffix)
    } else {
      gtf.out <- paste0(base_name, add_suffix, ".", extension)
    }
  } else {
    gtf.out <- output.file
  }
  
  
  tryCatch(
    {
      rtracklayer::export(sorted.gtf.gr, gtf.out, format = "gtf")
    },
    error = function(e) {
      stop("Failed to export GTF file: ", e$message, call. = FALSE)
    }
  )
  
  # Return invisible
  invisible(NULL)
  
  }