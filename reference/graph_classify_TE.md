# Create Classification Plots for TE Expression

Generates multiple types of plots showing the proportion of
self-expressed vs. gene-dependent transposable elements across different
TE classes (LINE, SINE, LTR, DNA, etc.).

## Usage

``` r
graph_classify_TE(
  res,
  plot.title = NULL,
  device = "png",
  width.pie = 14,
  height.pie = 7,
  width.bar = 7,
  height.bar = 7,
  output_folder = ".",
  colors = c(dependent = "#FF1F5B", self = "#009ADE"),
  labels = c(dependent = "Gene-dependent TEs", self = "Self-expressed TEs"),
  save = "all"
)
```

## Arguments

- res:

  Data frame containing TE classification results. Must include columns:
  `TE_expression` (classification: "dependent" or "self"), and
  `expression_type` (detailed classification for stacked bar plots).

- plot.title:

  Character string. Title for the plot. If NULL (default), no title is
  displayed.

- device:

  Character vector. File format(s) for output plots. Supported: "pdf",
  "svg", "eps", "png", "tiff", "jpeg". Default is "png".

- width.pie:

  Numeric. Plot width in inches for pie graph. Default is 14 (suitable
  for 4-6 TE classes).

- height.pie:

  Numeric. Plot height in inches for pie graph. Default is 7.

- width.bar:

  Numeric. Plot width in inches for stack. Default is 7 (suitable for
  4-6 TE classes).

- height.bar:

  Numeric. Plot height in inches for stack bar. Default is 7.

- output_folder:

  Character string. Directory where plots will be saved. Default is
  current directory (".").

- colors:

  Named character vector. Colors for TE expression types. Default is
  c("dependent" = "#FF1F5B", "self" = "#009ADE"). Names must match
  values in `TE_expression` column.

- labels:

  Named character vector. Labels for TE expression types in legend.
  Default is c("dependent" = "Gene-dependent TEs", "self" =
  "Self-expressed TEs").

- save:

  Character string. Which TEs to save in output file:

  - "all": All expressed TEs (default)

  - "dys": Only significantly dysregulated TEs (up or down)

  - "up": Only significantly upregulated TEs

  - "down": Only significantly downregulated TEs

  The function creates three types of plots:

  **1. Pie Charts:** One pie chart per TE class showing the overall
  proportion of:

  - **Self-expressed TEs**: Transcribed from their own promoter

  - **Gene-dependent TEs**: Transcribed as part of gene runthrough

  **2. Stacked Bar Charts:** Grouped bars showing detailed
  classification breakdown:

  - Exon, Intron (expressed), Intron (not expressed), etc.

  - Grouped by self-expressed vs. gene-dependent

  - One panel per TE class

  **3. Bar plot:** One bar per TE class showing the overall proportion
  and absolute number of TE elements en each class.

## Value

Invisible NULL. Called for side effects (creating plot files).

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming 'classified_TEs' has TE_expression and TE_class columns
graph_classify_TE(
  res = TE_results$res.TEs,
  plot.title = "TE Expression Classification",
  device = c("pdf", "png"),
  output_folder = "results/classification"
)
} # }

  
```
