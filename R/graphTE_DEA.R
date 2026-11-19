#' Plot Differential Expression Analysis Results for Transposable Elements
#'
#' Creates volcano and MA plots from DESeq2 differential expression results
#' and saves them in multiple file formats. 
#'
#' @param res Data frame or DESeqResults object containing differential expression
#'   results. Must include columns: \code{log2FoldChange}, \code{padj}, and
#'   \code{baseMean}.
#' @param maxpadj Numeric. Adjusted p-value threshold for significance.
#'   Features with padj < maxpadj are considered significant. Default is 0.05.
#' @param minlfc Numeric. Minimum absolute log2 fold change threshold for
#'   differential expression. Features with |log2FC| > minlfc are considered
#'   differentially expressed. Default is 1.
#' @param device Character vector. Output file format(s). Supported formats:
#'   "svg", "eps", "png", "tiff", "jpeg". Multiple formats can be
#'   specified as a vector (e.g., \code{c("tiff", "png")}). Default is "png".
#' @param output_folder Character string. Directory where plots will be saved.
#'   Must exist or be creatable. Default is current directory (".").
#' @param width Numeric. Plot width in inches. Default is 7.
#' @param height Numeric. Plot height in inches. Default is 7.
#' @param plot.title Character string. Title for both plots. If NULL (default),
#'   uses generic titles "Volcano Plot" and "MA Plot".
#' @param minBM Numeric. Only features with baseMean > minBM are shown in the
#'   volcano plot (the MA plot is not filtered). If greater than 0,
#'   "_BM_<minBM>" is appended to the volcano plot file name and the threshold
#'   is shown as its subtitle. Default is 0 (no filtering).
#'
#' @return Does not return any object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assume 'results' is a DESeq2 results object
#' library(DESeq2)
#'
#' # Basic usage - save as PNG in current directory
#' graphTE_DEA(
#'   res = results,
#'   maxpadj = 0.05,
#'   minlfc = 1
#' )
#'
#' # Save in multiple formats with custom title
#' graphTE_DEA(
#'   res = results,
#'   maxpadj = 0.01,
#'   minlfc = 1.5,
#'   device = c("png", "svg"),
#'   output_folder = "results/plots",
#'   width = 10,
#'   height = 8,
#'   plot.title = "KO vs WT"
#' )
#'
#' # Custom dimensions for publication
#' graphTE_DEA(
#'   res = results,
#'   maxpadj = 0.05,
#'   minlfc = 1,
#'   device = "jpeg",
#'   output_folder = "figures",
#'   width = 6,
#'   height = 5,
#'   plot.title = "Differential TE Expression"
#' )
#' }
#'
graphTE_DEA <- function(res,
                        maxpadj = 0.05,
                        minlfc = 1,
                        device = "png",
                        output_folder = ".",
                        width = 7,
                        height = 7,
                        plot.title = NULL,
                        minBM = 0) {
  
  # Input validation 
  if (missing(res)) {
    stop("Argument 'res' is missing with no default.", call. = FALSE)
  }
  
  if (!is.data.frame(res)) {
    if (inherits(res, "DESeqResults")) {
      res <- as.data.frame(res)
    } else {
      stop("'res' must be a data frame or DESeqResults object.", call. = FALSE)
    }
  }

  if (!is.numeric(minBM) || length(minBM) != 1 || is.na(minBM) || minBM < 0) {
    stop("'minBM' must be a single non-negative number.", call. = FALSE)
  }

  # Keep only features with baseMean over minBM for the volcano plot
  res.volcano <- res
  if (minBM > 0) {
    res.volcano <- res[res$baseMean > minBM, , drop = FALSE]
  }
  bm.suffix <- if (minBM > 0) paste0("_BM_", minBM) else ""
  bm.subtitle <- if (minBM > 0) paste0("baseMean > ", minBM) else NULL

  # Create output folder if it doesn't exist
  if (!dir.exists(output_folder)) {
    tryCatch(
      {
        dir.create(output_folder, recursive = TRUE)
      },
      error = function(e) {
        stop(
          "Failed to create output directory '", output_folder, "': ",
          e$message,
          call. = FALSE
        )
      }
    )
  }
  
  # Set default titles if not provided
  volcano_title <- if (is.null(plot.title)) "Volcano Plot" else plot.title
  ma_title <- if (is.null(plot.title)) "MA Plot" else plot.title
  
  # Create volcano plot
  p.volcano <- tryCatch(
    volcanoPlot(res.volcano, maxpadj, minlfc, volcano_title) +
      ggplot2::labs(subtitle = bm.subtitle),
    error = function(e) {
      stop(
        "Failed to create volcano plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Create MA plot
  p.maplot <- tryCatch(
    MAPlot(res, maxpadj, minlfc, ma_title),
    error = function(e) {
      stop(
        "Failed to create MA plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Save volcano plot
  filename.volcano <- file.path(output_folder, paste0("volcanoPlot", bm.suffix))
  
  tryCatch(
    {
      printDevice(p.volcano, filename.volcano, device, 
                  width = width, height = height)
    },
    error = function(e) {
      stop(
        "Failed to save volcano plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Save MA plot
  filename.maplot <- file.path(output_folder, "maPlot")
  
  tryCatch(
    {
      printDevice(p.maplot, filename.maplot, device, 
                  width = width, height = height)
    },
    error = function(e) {
      stop(
        "Failed to save MA plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Don't return anything
  invisible()
}
