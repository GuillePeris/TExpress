# Create violin plot for specific TE list

Generates a violin plot showing log2 fold changes for a specified list
of transposable elements (TEs).

## Usage

``` r
violinPlotByTEList(
  res.TEs,
  TE_list,
  broad_type = NULL,
  minlfc = 1,
  maxpadj = 0.05,
  min.N = 10,
  width = 7,
  height = 7,
  device = "png",
  output_folder = ".",
  output_file = NULL,
  plot.title = "Violin plot",
  onlySignif = FALSE
)
```

## Arguments

- res.TEs:

  Data frame containing TE differential expression results with columns:
  log2FoldChange, padj, TE_element, TE_name, TE_family, TE_class

- TE_list:

  Character vector of TE identifiers to plot. All elements should belong
  to the same hierarchical level

- broad_type:

  Character string specifying the TE classification level (e.g.,
  "TE_class", "TE_family"). If NULL, automatically detected from TE_list

- minlfc:

  Numeric minimum log2 fold change threshold for significance (default:
  1)

- maxpadj:

  Numeric maximum adjusted p-value threshold for significance (default:
  0.05)

- min.N:

  Integer minimum number of elements required per group (default: 10)

- width:

  Numeric plot width in inches (default: 7)

- height:

  Numeric plot height in inches (default: 7)

- device:

  Character string specifying output device (default: "png")

- output_folder:

  Character string path to output directory (default: ".")

- output_file:

  Character string with the output file name, without extension (added
  from `device`). If NULL, uses "violinPlot\_\<broad_type\>" (default:
  NULL)

- plot.title:

  Character string for plot title (default: "Violin plot")

- onlySignif:

  Logical. If TRUE, violins are drawn only from features with padj \<=
  maxpadj, "\_signif" is appended to the file name and the threshold is
  added to the subtitle; if FALSE, from all features. Up- and
  down-regulated points are plotted in both cases (default: FALSE)

## Value

Invisibly returns NULL. Saves plot to file as side effect.

## Examples

``` r
if (FALSE) { # \dontrun{
violinPlotByTEList(res.TEs, TE_list = c("LINE", "SINE", "LTR", "DNA"),
                   minlfc = 1.5, maxpadj = 0.01)

violinPlotByTEList(res.TEs, TE_list = c("Alu", "L1", "ERV1", "SVA"),
                   minlfc = 1.5, maxpadj = 0.01)
} # }
```
