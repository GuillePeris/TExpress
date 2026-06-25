# Shared, publication-oriented theme and colour palettes for all TExpress plots.
# Centralising these keeps the figures visually consistent (fonts, sizes, title
# styling, panel look, colours) so they form a cohesive, journal-ready set.

#' Publication theme for TExpress plots
#'
#' A clean, volcano-plot-style ggplot2 theme: white background, a black box
#' (panel border) around the plotting area, no grid lines, black ticks and
#' black, bold, centred title/subtitle. Used by every plotting function in the
#' package so the figures look consistent.
#'
#' @param base_size Base font size in points.
#' @param aspect Optional numeric aspect ratio for the panel (e.g. 1 for the
#'   square volcano/MA plots). If NULL (default) the aspect ratio is left free.
#'
#' @return A ggplot2 theme object.
#' @keywords internal
#' @noRd
.theme_texpress <- function(base_size = 13, aspect = NULL) {

  th <- ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      # Background + box around the plot, no grid
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA,
                                           linewidth = 0.6),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),

      # Axes
      axis.text  = ggplot2::element_text(colour = "grey20", size = base_size),
      axis.title = ggplot2::element_text(colour = "grey20", size = base_size + 2),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = ggplot2::unit(0.2, "cm"),

      # Titles — black, centred, never the old red
      plot.title = ggplot2::element_text(
        colour = "black", face = "bold", hjust = 0.5,
        size = base_size + 3,
        margin = ggplot2::margin(b = 6)
      ),
      plot.subtitle = ggplot2::element_text(
        colour = "grey25", hjust = 0.5,
        size = base_size,
        margin = ggplot2::margin(b = 8)
      ),

      # Legend
      legend.title    = ggplot2::element_blank(),
      legend.text     = ggplot2::element_text(size = base_size - 1),
      legend.key      = ggplot2::element_blank(),
      legend.position = "right"
    )

  if (!is.null(aspect)) {
    th <- th + ggplot2::theme(aspect.ratio = aspect)
  }

  th
}

# ---- Shared colour palettes -------------------------------------------------

# Genomic regions (sequential earth tones, kept from the original design).
.region_palette <- c(
  "Promoter"   = "#005F73",
  "5' UTR"     = "#0A9396",
  "Exon"       = "#94D2BD",
  "Intron"     = "#E9D8A6",
  "3' UTR"     = "#EE9B00",
  "Downstream" = "#BB3E03",
  "Intergenic" = "#870000"
)

# TE classes (kept from the original design; reused across stacked bars and
# the self-expressed class barplot for cross-figure consistency).
.te_class_palette <- c(
  "LINE"  = "#e64b35",
  "LTR"   = "#4dbbd5",
  "SINE"  = "#00a087",
  "DNA"   = "#3c5488",
  "Other" = "#f39b7f"
)

# Up-/down-regulation accents used by volcano, MA and violin rains.
.expr_palette <- c(
  "up"   = "#990f02",
  "down" = "#0f4392",
  "ns"   = "grey60"
)

# Single accent colour for barplots that do not encode colour meaningfully.
.bar_accent <- "#0A9396"
