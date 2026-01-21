#' Create Pie Charts for TE Expression Classification
#'
#' Generates different types of charts showing the proportion of self-expressed 
#' vs. gene-dependent transposable elements across different TE classes (LINE,
#' SINE, LTR, DNA, etc.).
#'
#' @param res Data frame containing TE classification results. 
#' @param plot.title Character string. Title for the plot. If NULL (default),
#'   no title is displayed.
#' @param device Character vector. File format(s) for output plots. Supported:
#'   "pdf", "svg", "eps", "png", "tiff", "jpeg". Default is "png".
#' @param width.pie Numeric. Plot width in inches for pie graph. Default is 14 (suitable for
#'   4-6 TE classes).
#' @param height.pie Numeric. Plot height in inches for pie graph. Default is 7.
#' @param height.bar Numeric. Plot height in inches for stack bar. Default is 7.
#' @param width.bar Numeric. Plot width in inches for stack. Default is 7 (suitable for
#'   4-6 TE classes).
#' @param output_folder Character string. Directory where plots will be saved.
#'   Default is current directory (".").
#' @param colors Named character vector. Colors for TE expression types.
#'   Default is c("dependent" = "#FF1F5B", "self" = "#009ADE"). Names must
#'   match values in \code{TE_expression} column.
#' @param labels Named character vector. Labels for TE expression types in
#'   legend. Default is c("dependent" = "Gene-dependent TEs",
#'   "self" = "Self-expressed TEs").
#' @param save Character string. Which TEs to save in output file:
#'   \itemize{
#'     \item "all": All expressed TEs (default)
#'     \item "dys": Only significantly dysregulated TEs (up or down)
#'     \item "up": Only significantly upregulated TEs
#'     \item "down": Only significantly downregulated TEs
#'   }
#'
#' @details
#' The function creates a pie chart for each TE class showing:
#' \itemize{
#'   \item **Self-expressed TEs**: TEs transcribed from their own promoter
#'   \item **Gene-dependent TEs**: TEs transcribed as part of gene runthrough
#' }
#'
#' Each pie slice is labeled with its percentage. TE classes are displayed
#' as separate facets, arranged horizontally (up to \code{max_cols} columns).
#'
#' @return Invisible NULL. Called for side effects (creating plot files).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming 'classified_TEs' has TE_expression and TE_class columns
#' TE_classify_pie(
#'   res = TE_results$res.TEs,
#'   plot.title = "TE Expression Classification",
#'   device = c("pdf", "png"),
#'   output_folder = "results/classification"
#' )
#' }
#'
#'   
TE_classify_pie <- function(res,
                            plot.title = NULL,
                            device = "png",
                            width.pie = 14,
                            height.pie = 7,
                            width.bar = 7,
                            height.bar = 7,
                            output_folder = ".",
                            colors = c("dependent" = "#FF1F5B", "self" = "#009ADE"),
                            labels = c("dependent" = "Gene-dependent TEs",
                                       "self" = "Self-expressed TEs"),
                            save = "all") {
  
  # Input validation
  if (missing(res)) {
    stop("Argument 'res' is missing with no default.", call. = FALSE)
  }
  
  if (!is.data.frame(res)) {
    stop("'res' must be a data frame.", call. = FALSE)
  }
  
  if (nrow(res) == 0L) {
    stop("'res' contains no data rows.", call. = FALSE)
  }
  
  # Check required columns
  required_cols <- c("TE_expression", "TE_class")
  missing_cols <- setdiff(required_cols, colnames(res))
  
  if (length(missing_cols) > 0L) {
    stop(
      "'res' must contain columns: ",
      paste(required_cols, collapse = ", "),
      "\nMissing: ", paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!dir.exists(output_folder)) {
    tryCatch(
      dir.create(output_folder, recursive = TRUE),
      error = function(e) {
        stop(
          "Cannot create output directory '", output_folder, "': ",
          e$message,
          call. = FALSE
        )
      }
    )
  }
  
  
  # Create pie plots
  filename <- file.path(output_folder, paste0("pie_TE_classes_", save))
  
  p <- tryCatch(
    create_pie_plot(res = res, 
                      plot.title = plot.title,
                      colors = colors,
                      labels = labels
                      ),
    error = function(e) {
      stop(
        "Failed to create pie TE classify plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Save plot 
  tryCatch(
    printDevice(p, filename, device,
                width = width.pie, height = height.pie),
    error = function(e) {
      warning(
        "Failed to save pie TE classify plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Create stack bar plots
  filename <- file.path(output_folder, paste0("stack_bar_TE_classes_", save))
  
  p <- tryCatch(
    create_grouped_stack_bar(res = res, 
                              plot.title = plot.title
                             ),
    error = function(e) {
      stop(
        "Failed to create stack bar TE classify plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  # Save plot 
  tryCatch(
    printDevice(p, filename, device,
                width = width.bar, height = height.bar),
    error = function(e) {
      warning(
        "Failed to save pie TE classify plot: ", e$message,
        call. = FALSE
      )
    }
  )
  
  invisible(NULL)
}

#' Create Pie Plot for TE Classification
#'
#' Internal function to generate the actual pie chart visualization.
#'
#' @param res Data frame with TE_expression and TE_class columns
#' @param plot.title Character string for plot title
#' @param colors Named character vector of colors
#' @param labels Named character vector of legend labels
#'
#' @return A ggplot object
#' @keywords internal
#' @noRd
#'
create_pie_plot <- function(res,
                             plot.title = NULL,
                             colors = c("dependent" = "#FF1F5B", "self" = "#009ADE"),
                             labels = c("dependent" = "Gene-dependent TEs",
                                        "self" = "Self-expressed TEs")) {

  # Calculate proportions by TE class
  res_summary <- res %>%
    dplyr::count(.data$TE_expression, .data$TE_class) %>%
    dplyr::group_by(.data$TE_class) %>%
    dplyr::reframe(
      TE_class = .data$TE_class,
      TE_expression = .data$TE_expression,
      prop = .data$n / sum(.data$n) * 100
    )
  
  res_summary$TE_expression <- factor(res_summary$TE_expression)
  nTE_class <- length(unique(res$TE_class))
  
  p <- ggplot2::ggplot(
    res_summary,
    ggplot2::aes(x = "", y = .data$prop, fill = .data$TE_expression)
  ) + 
    ggplot2::geom_bar(width = 1, stat = "identity", color = "white") +
    ggplot2::coord_polar("y", start = 0) +
    ggplot2::theme_void(base_size = 24) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$TE_class), 
      ncol = nTE_class
    )  +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%0.1f%%", .data$prop)), 
                       position = ggplot2::position_stack(vjust=0.5),
                       size=6) +
    ggplot2::scale_fill_manual(
      labels = labels,
      values = colors
    ) + 
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 20)
    ) 

  # Add title if provided
  if (!is.null(plot.title)) {
    p <- p + ggplot2::ggtitle(plot.title)
  }
  
  p
}


create_grouped_stack_bar <- function(res,
                             plot.title = NULL) {
 
  # Calculate proportions by TE class
  res_summary <- res %>%
    dplyr::count(.data$TE_expression, 
                 .data$expression_type, 
                 .data$TE_class) %>%
    dplyr::group_by(.data$TE_class) %>%
    dplyr::reframe(
      TE_class = .data$TE_class,
      TE_expression = .data$TE_expression,
      expression_type = .data$expression_type,
      prop = .data$n / sum(.data$n) * 100
    )
  
  res_summary$TE_expression <- factor(res_summary$TE_expression)
  res_summary$expression_type <- factor(res_summary$expression_type)
  levels(res_summary$TE_expression)[levels(res_summary$TE_expression) == "dependent"] <- "dep."
  nTE_class <- length(unique(res$TE_class))
  
  p <- ggplot2::ggplot(
    res_summary,
    ggplot2::aes(x = .data$TE_expression, 
                 y = .data$prop, 
                 fill = .data$expression_type,
                 ymin = 0,
                 ymax = max(.data$prop + 5, 100)
    )
  ) + 
    ggplot2::geom_bar(width = 1, stat = "identity", color = "white") +
    ggplot2::facet_grid(
      cols = ggplot2::vars(.data$TE_class)
  )   
  
  
  # Add percentage labels
  p <- p +
    ggfittext::geom_fit_text(
      ggplot2::aes(label = sprintf("%0.1f%%", .data$prop)),
      min.size = 8,
      size = 16,
      position = ggplot2::position_stack(vjust = 0.5),
      contrast = TRUE,
      show.legend = FALSE
    )
  
  # Apply theme
  p <- p +
    ggplot2::theme(
      # Grid
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),

      # Background
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      
      # Axis
      axis.text = ggplot2::element_text(colour = "grey30", size = 12),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(
        size = 20,
        colour = "grey30",
        hjust = 0.5,
        margin = ggplot2::margin(r = 10)
      ),
      axis.ticks.y = ggplot2::element_line(color = "grey30", linewidth = 0.5),
      axis.ticks.length.y = ggplot2::unit(0.25, "cm"),
      
      # Legend
      legend.text = ggplot2::element_text(size = 12),
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 24, color = "red")
    ) +
    ggplot2::labs(
      y = "% TE loci"
    )
  
  # Format y-axis as percentages
  p <- p +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(x, " %"),
      expand = c(0, 1)
    )
  
  
  # Add title if provided
  if (!is.null(plot.title)) {
    p <- p + ggplot2::ggtitle(plot.title)
  }
  
  p 
}

create_pie_donut <- function(res,
                             TE_feature,
                             plot.title = NULL,
                             output_folder = ".",
                             height = 7,
                             width = 7,
                             prefix = NULL) {
  
  # Validate input
  
  if (!dir.exists(output_folder)) {
    tryCatch(
      dir.create(output_folder, recursive = TRUE),
      error = function(e) {
        stop(
          "Cannot create output directory '", output_folder, "': ",
          e$message,
          call. = FALSE
        )
      }
    )
  }
  
  # Get TE_feature type
  TE_type <- .determine_broad_type(res, TE_feature)
  
  # Filter by TE_feature
  res_summary <- res %>% 
    dplyr::filter(.data[[TE_type]] == TE_feature)
  names(res_summary)[names(res_summary) == "TE_expression"] <- TE_feature
  
  res_summary$expression_type <- droplevels(res_summary$expression_type)
  
  
  # Create pie plots
  if(!is.null(prefix) & !startsWith(prefix, "_")) {
    prefix <- paste0("_", prefix)
  }
  
  filename <- file.path(output_folder, 
                        paste0("donut_pie_", TE_feature, prefix, ".png"))
  
  grDevices::png(filename, width=width, height=height, res=300, units = "in")
  webr::PieDonut(res_summary,
                 ggplot2::aes(pies=!!TE_feature,
                              donuts=!!quote(expression_type)), 
           showRatioThreshold = 0.001, 
           labelposition=0, 
           pieLabelSize = 8, donutLabelSize = 5,  showPieName=TRUE,
           ratioByGroup=FALSE, titlesize = 10)
  
  grDevices::dev.off()
}
