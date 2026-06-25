#' Create Stacked Bar Plot for Region or TE Class Distribution
#'
#' Internal function to create stacked bar plots showing percentage distribution
#' of genomic regions or TE classes across samples.
#'
#' @param aCR Data frame with columns: region (or class), percent, clon
#' @param plot.title Character string. Plot title.
#' @param analysis Character vector. Labels for x-axis (sample names)
#' @param plot.type Character. Either "region" or "TEclass" to determine
#'   categories and color scheme.
#'
#' @importFrom rlang .data
#' @return A ggplot object.
#'
#' @keywords internal
#' @noRd
stackBarPlot <- function(aCR,
                         plot.title,
                         analysis,
                         plot.type = c("region", "TEclass")) {
  
  plot.type <- match.arg(plot.type)
  
  # Input validation
  if (!is.data.frame(aCR)) {
    stop("'aCR' must be a data frame.", call. = FALSE)
  }
  
  required_cols <- c("region", "percent", "clon")
  if (!all(required_cols %in% colnames(aCR))) {
    stop(
      "'aCR' must contain columns: ",
      paste(required_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  # Define levels and colors based on plot type (shared palettes from
  # R/graph_theme.R keep colours consistent across all figures).
  if (plot.type == "region") {
    my.levels <- names(.region_palette)
    palette <- .region_palette
  } else {  # TEclass
    my.levels <- names(.te_class_palette)
    palette <- .te_class_palette

    # Collapse non-standard classes to "Other"
    aCR$region[!(aCR$region %in% my.levels[1:4])] <- "Other"
  }
  
  # Set factor levels
  aCR$region <- factor(aCR$region, levels = my.levels)
  
  # Create plot using internal helper
  p <- .create_stack_bar(
    data = aCR,
    x = "clon",
    y = "percent",
    fill = "region",
    palette = palette,
    x_labels = analysis,
    title = plot.title
  )
  
  p
}

#' Internal Helper for Creating Stack Bar Plots
#'
#' @keywords internal
#' @noRd
.create_stack_bar <- function(data, x, y, fill, palette, x_labels, title) {
  
  # Create base plot
  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(fill = .data[[fill]], y = .data[[y]], x = .data[[x]])
  ) +
    ggplot2::geom_bar(position = "fill", stat = "identity", width = 0.95) +
    ggplot2::scale_x_discrete(
      breaks = levels(data[[x]]),
      labels = x_labels,
      expand = c(0, 0)
    ) +
    ggplot2::scale_fill_manual(values = palette)
  
  # Apply shared publication theme
  p <- p +
    .theme_texpress(base_size = 14) +
    ggplot2::theme(
      aspect.ratio = 1,
      axis.title.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = title,
      y = "% TE loci"
    )
  
  # Format y-axis as percentages
  p <- p +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(),
      expand = c(0, 0)
    )
  
  # Add percentage labels
  p <- p +
    ggfittext::geom_fit_text(
      ggplot2::aes(label = sprintf("%0.1f%%", .data[[y]])),
      stat = "identity",
      min.size = 8,
      size = 10,
      position = ggplot2::position_fill(vjust = 0.5),
      contrast = TRUE,
      show.legend = FALSE
    )
  
  p
}
