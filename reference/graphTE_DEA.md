# Plot Differential Expression Analysis Results for Transposable Elements

Creates volcano and MA plots from DESeq2 differential expression results
and saves them in multiple file formats.

## Usage

``` r
graphTE_DEA(
  res,
  maxpadj = 0.05,
  minlfc = 1,
  device = "png",
  output_folder = ".",
  width = 7,
  height = 7,
  plot.title = NULL
)
```

## Arguments

- res:

  Data frame or DESeqResults object containing differential expression
  results. Must include columns: `log2FoldChange`, `padj`, and
  `baseMean`.

- maxpadj:

  Numeric. Adjusted p-value threshold for significance. Features with
  padj \< maxpadj are considered significant. Default is 0.05.

- minlfc:

  Numeric. Minimum absolute log2 fold change threshold for differential
  expression. Features with \|log2FC\| \> minlfc are considered
  differentially expressed. Default is 1.

- device:

  Character vector. Output file format(s). Supported formats: "svg",
  "eps", "png", "tiff", "jpeg". Multiple formats can be specified as a
  vector (e.g., `c("tiff", "png")`). Default is "png".

- output_folder:

  Character string. Directory where plots will be saved. Must exist or
  be creatable. Default is current directory (".").

- width:

  Numeric. Plot width in inches. Default is 7.

- height:

  Numeric. Plot height in inches. Default is 7.

- plot.title:

  Character string. Title for both plots. If NULL (default), uses
  generic titles "Volcano Plot" and "MA Plot".

## Value

Does not return any object

## Examples

``` r
if (FALSE) { # \dontrun{
# Assume 'results' is a DESeq2 results object
library(DESeq2)

# Basic usage - save as PNG in current directory
graphTE_DEA(
  res = results,
  maxpadj = 0.05,
  minlfc = 1
)

# Save in multiple formats with custom title
graphTE_DEA(
  res = results,
  maxpadj = 0.01,
  minlfc = 1.5,
  device = c("png", "svg"),
  output_folder = "results/plots",
  width = 10,
  height = 8,
  plot.title = "KO vs WT"
)

# Custom dimensions for publication
graphTE_DEA(
  res = results,
  maxpadj = 0.05,
  minlfc = 1,
  device = "jpeg",
  output_folder = "figures",
  width = 6,
  height = 5,
  plot.title = "Differential TE Expression"
)
} # }
```
