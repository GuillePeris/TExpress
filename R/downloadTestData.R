#' Download Test Data for Package Examples and Testing
#'
#' Downloads sample data files from Dropbox for testing and running package
#' examples. The data includes human genome annotations (chromosome 1 only)
#' and TE count data, suitable for demonstrating differential expression
#' analysis workflows.
#'
#' @param folder Character string. Directory where test data will be downloaded.
#'   If NULL (default), uses a temporary directory via \code{tempdir()}.
#'   The directory will be created if it doesn't exist.
#' @param timeout Numeric. Download timeout in seconds. Default is 300 (5 minutes).
#'   Increase if downloads fail due to slow connections.
#'
#' @details
#' The function downloads three datasets:
#' \enumerate{
#'   \item Gene GTF (\strong{~78 MB}): Human gene annotations for chromosome 1
#'     from Ensembl GRCh38.115
#'   \item TE GTF (\strong{~10 MB}): Transposable element annotations for chromosome 1
#'     from RepeatMasker (hg38)
#'   \item Count Data (\strong{~2 MB compressed}): Example TEtranscripts count
#'     files in a ZIP archive. Contains 6 sample count files (3 control, 3 treatment)
#'     and a csv file for TEexpress analysis
#' }
#'
#' Files are downloaded from Dropbox using direct download links. The count data
#' ZIP file is automatically extracted and the archive is removed after extraction.
#'
#' @section Downloaded Files:
#' The following files will be created in the specified folder:
#' \itemize{
#'   \item \code{Homo_sapiens.GRCh38.115.chr1.gtf} - Gene annotations
#'   \item \code{hg38_rmsk_TE_chr1.gtf} - TE annotations
#'   \item \code{*.cntTable} - Multiple count files (extracted from ZIP)
#'   \item \code{data.csv} - csv file including info on count files
#' }
#'
#' @section File Sizes:
#' Total download size is approximately 90 MB. Ensure adequate disk space and
#' internet bandwidth. On slow connections, increase the \code{timeout} parameter.
#'
#' @return A list with four elements:
#' \describe{
#'   \item{folder}{Character string with path to download directory}
#'   \item{metafile}{Full path to metadata CSV file (data.csv)}
#'   \item{gtf.gene.file}{Full path to gene GTF file}
#'   \item{gtf.TE.file}{Full path to TE GTF file}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Download to temporary directory (default)
#' my.data <- downloadTestData()
#' 
#' # Download to specific directory
#' my.data <- downloadTestData(folder = "test_data")
#' 
#' # Use downloaded data in analysis
#' results <- TE_DEA(
#'   metafile = my.data$metafile,
#'   folder = my.data$folder,
#'   gtf.TE.file = my.data$gtf.TE.file
#' )
#'}
downloadTestData <- function(folder = NULL,
                             timeout = 300) {
  
  # ============================================================
  # Input Validation
  # ============================================================
  
  # Validate folder
  if (!is.null(folder)) {
    target_dir <- folder
  } else {
    target_dir <- tempdir()
  }
  
  # ============================================================
  # Setup
  # ============================================================
  
  # Create directory if it doesn't exist
  if (!dir.exists(target_dir)) {
    tryCatch(
      dir.create(target_dir, recursive = TRUE, showWarnings = FALSE),
      error = function(e) {
        stop(
          "Failed to create directory '", target_dir, "': ",
          e$message,
          call. = FALSE
        )
      }
    )
  }
  
  # Set download timeout
  options(timeout = timeout)
  
  # ============================================================
  # Define Files to Download
  # ============================================================
  
  files_to_download <- list(
    gene_gtf = list(
      filename = "Homo_sapiens.GRCh38.115_chr22.gtf",
      url = "https://www.dropbox.com/scl/fi/ehzjxyixfag0gvjbnnjsr/Homo_sapiens.GRCh38.115_chr22.gtf?rlkey=1x484mjw3xkzz58h2ibaeav7z&st=dl1yqy3x&dl=1",
      description = "Gene annotations (Ensembl GRCh38.115, chr1)",
      size_mb = "~78 MB"
    ),
    te_gtf = list(
      filename = "hg38_rmsk_TE_chr22.gtf",
      url = "https://www.dropbox.com/scl/fi/nyb0d8bexyq8zjy93vyt2/hg38_rmsk_TE_chr22.gtf?rlkey=6hdwa59zzpjwg26iwiuwq6amz&st=fw0rurph&dl=1",
      description = "TE annotations (RepeatMasker hg38, chr1)",
      size_mb = "~10 MB"
    ),
    count_data = list(
      filename = "countData.zip",
      url = "https://www.dropbox.com/scl/fi/tcodv0hi1m934pdsg5xcj/countData.zip?rlkey=2mh3z3vrhxetnuoeecdcld67i&st=mhddcfog&dl=1",
      description = "Count data files (6 samples)",
      size_mb = "~2 MB",
      extract = TRUE,
      remove_after = TRUE
    )
  )
  
  # ============================================================
  # Download Files
  # ============================================================
 
  for (file_info in files_to_download) {
    dest_file <- file.path(target_dir, file_info$filename)
    
    download_success <- tryCatch(
      {
        downloader::download(
          file_info$url,
          dest_file,
          mode = "wb"
        )
        TRUE
      },
      error = function(e) {
        message("  [ERROR] Download failed: ", e$message)
        FALSE
      },
      warning = function(w) {
        message("  [WARNING] ", w$message)
        TRUE  # Continue despite warning
      }
    )
    
    # Extract if ZIP file
    if (!is.null(file_info$extract) && file_info$extract) {
      extract_success <- tryCatch(
        {
          utils::unzip(
            zipfile = dest_file,
            exdir = target_dir,
            overwrite = TRUE
          )
          TRUE
        },
        error = function(e) {
          message("  [ERROR] Failed to extract ZIP: ", e$message)
          FALSE
        }
      )
      
      if (extract_success) {
        # Remove ZIP file if requested
        if (!is.null(file_info$remove_after) && file_info$remove_after) {
          tryCatch(
            {
              file.remove(dest_file)
            },
            error = function(e) {
            }
          )
        }
      }
    }
    
  }
  
  # List to return
  return_list <- list(
    metafile = file.path(target_dir, "data.csv"),
    folder = target_dir,
    gtf.gene.file = file.path(target_dir, files_to_download$gene_gtf$filename),
    gtf.TE.file = file.path(target_dir,files_to_download$te_gtf$filename)
  )
 
  return_list

}
