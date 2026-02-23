filter_GTF <- function(gtf.genes.file, 
                               features = c("protein_coding"),
                               format.in = "gtf",
                               format.out = "gtf",
                               suffix = "filtered") {
  
  # Validate gtf.genes.file
  if (missing(gtf.genes.file)) {
    stop("Argument 'gtf.genes.file' is missing with no default.", call. = FALSE)
  }
  
  if (!file.exists(gtf.genes.file)) {
    stop("Gene GTF file '", gtf.genes.file, "' not found.", call. = FALSE)
  }
  
  # Import gtf file
  gtf.genes <- tryCatch(
    importGTF(gtf.genes.file, format = format.in),
    error = function(e) {
      stop("Failed to import gene GTF file: ", e$message, call. = FALSE)
    }
  )
  
  # Check required columns exist
  gtf_mcols <- GenomicRanges::mcols(gtf.genes)
  required_gtf_cols <- c("type", "transcript_biotype", "gene_biotype")
  missing_gtf_cols <- setdiff(required_gtf_cols, colnames(gtf_mcols))
  
  if (length(missing_gtf_cols) > 0L) {
    stop(
      "Gene GTF is missing required column(s): ",
      paste(missing_gtf_cols, collapse = ", "),
      "\nRequired columns: ", paste(required_gtf_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  
  # Filter to protein-coding genes
  filter_features <- 
    (gtf.genes$type == "gene" & 
       gtf.genes$gene_biotype %in% features ) |
    (gtf.genes$type != "gene" &
       gtf.genes$transcript_biotype %in% features ) 
    
  gtf.genes <- gtf.genes[filter_features]
  
  # Export gtf file
  base_name <- tools::file_path_sans_ext(gtf.genes.file)
  extension <- tools::file_ext(gtf.genes.file)
  output <- paste0(base_name, "_", suffix, ".", format.out)
  tryCatch(
    rtracklayer::export(gtf.genes, output, format = format.out),
    error = function(e) {
      stop("Failed to export gene GTF file: ", e$message, call. = FALSE)
    }
  )
  
}
