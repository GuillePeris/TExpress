#' Read Metadata File with Sample Information
#'
#' Reads and validates a tab-separated metadata file containing sample information
#' for differential expression analysis. The file must contain specific columns
#' and follow formatting requirements for DESeq2 compatibility.
#'
#' @param datafile Character string. Path to a tab-separated file (.txt, .tsv, or .csv)
#'   containing metadata. Must have 4 columns with headers: File, Sample,
#'   Group, and Condition, plus an optional fifth column with header Batch.
#' @param batch Logical. If TRUE, the metadata file must contain a fifth column
#'   with header Batch, used to model batch effects in DESeq2, and every batch
#'   must contain samples of both conditions. If FALSE, a Batch column is
#'   accepted but not required. Default is FALSE.
#'
#' @details
#' The metadata file must follow these requirements:
#' \itemize{
#'   \item 4 columns with headers: File, Sample, Group, Condition, optionally
#'     followed by a fifth column with header Batch (required if
#'     \code{batch = TRUE})
#'   \item At least one row of data (beyond the header)
#'   \item The "Condition" column must contain exactly two unique values: "Control" and "Treat"
#'   \item No duplicate values in the "Sample" or "File" columns
#'   \item Tab-separated format (other delimiters not supported)
#' }
#'
#' The "Control" and "Treat" labels in the Condition column are required for
#' proper DESeq2 analysis, where "Control" will be set as the reference level.
#'
#' @return A data frame with four columns (File, Sample, Group, Condition), plus
#'   Batch if present, and validated metadata. The Condition column is converted
#'   to a factor with "Control" as the reference level.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Read metadata from file
#' metadata <- read_metadata("path/to/metadata.txt")
#'
#' # Example valid file structure:
#' # File            Sample  Group   Condition
#' # WT1.cntTable    WT1     WT      Control
#' # WT2.cntTable    WT2     WT      Control
#' # KO1.cntTable    KO1     KO      Treat
#' # KO2.cntTable    KO2     KO      Treat
#'
#' # Read metadata including a Batch column
#' metadata <- read_metadata("path/to/metadata_batch.txt", batch = TRUE)
#'
#' # Example valid file structure with batch:
#' # File            Sample  Group   Condition  Batch
#' # WT1.cntTable    WT1     WT      Control    1
#' # WT2.cntTable    WT2     WT      Control    2
#' # KO1.cntTable    KO1     KO      Treat      1
#' # KO2.cntTable    KO2     KO      Treat      2
#' }
read_metadata <- function(datafile, batch = FALSE) {

  # Input validation
  if (missing(datafile)) {
    stop("Argument 'datafile' is missing with no default.")
  }

  if (!isTRUE(batch) && !isFALSE(batch)) {
    stop("'batch' must be TRUE or FALSE.", call. = FALSE)
  }

  # Read and validate metadata file
  if (!file.exists(datafile)) {
    stop("File '", datafile, "' not found.")
  }
  
  # Read metadata file
  metadata <- tryCatch(
    utils::read.table(
      datafile,
      sep = "\t",
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      comment.char = "",
      quote = "\"",
      na.strings = c("NA", "")
    ),
    error = function(e) {
      stop("Failed to read file '", datafile, "': ", e$message, call. = FALSE)
    }
  )  
  # Validate structure. An optional fifth column must be named Batch.
  expected_cols <- c("File", "Sample", "Group", "Condition")
  has_batch <- ncol(metadata) == length(expected_cols) + 1 &&
    colnames(metadata)[length(expected_cols) + 1] == "Batch"

  if (batch && !has_batch) {
    stop("'batch = TRUE' requires a fifth column named 'Batch' in file '",
         datafile, "'.", call. = FALSE)
  }

  if (has_batch) {
    expected_cols <- c(expected_cols, "Batch")
  }

  if (ncol(metadata) != length(expected_cols)) {
    stop("File '", datafile, "' should have ", length(expected_cols),
         " columns, but has ", ncol(metadata), ".")
  }

  # Columns are interpreted by position (File, Sample, Group, Condition[, Batch]).
  colnames(metadata) <- expected_cols
  
  # Validate Condition column
  unique_conditions <- unique(metadata$Condition)
  required_conditions <- c("Control", "Treat")
  
  if (!setequal(unique_conditions, required_conditions)) {
    stop(
      "The 'Condition' column must contain only 'Control' and 'Treat'.\nFound: ",
      paste(unique_conditions, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Reject duplicate sample or file identifiers.
  for (col in c("Sample", "File")) {
    dups <- unique(metadata[[col]][duplicated(metadata[[col]])])
    if (length(dups) > 0) {
      stop(
        "The '", col, "' column must not contain duplicates.\nDuplicated: ",
        paste(dups, collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  # Set Condition as a factor with 'Control' as the reference level.
  metadata$Condition <- factor(metadata$Condition, levels = required_conditions)

  # Each batch needs samples of both conditions; otherwise the batch term
  # absorbs those samples and they no longer inform the condition effect.
  if (batch) {
    batch_table <- table(metadata$Batch, metadata$Condition)
    bad <- rownames(batch_table)[apply(batch_table, 1, min) < 1]
    if (length(bad) > 0) {
      stop("With 'batch = TRUE', every batch must contain samples of both ",
           "conditions. Batches failing: ",
           paste0(bad, " (Control=", batch_table[bad, "Control"],
                  ", Treat=", batch_table[bad, "Treat"], ")", collapse = "; "),
           ".", call. = FALSE)
    }
  }

  metadata
}