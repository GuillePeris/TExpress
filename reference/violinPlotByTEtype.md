# Create violin plot for specific TE type

Generates a violin plot showing the most dysregulated transposable
elements within a specific TE type, organized by a finer classification
level.

## Usage

``` r
violinPlotByTEtype(
  res.TEs,
  TE_type,
  specific_type = "TE_name",
  broad_type = NULL,
  order = "up",
  minlfc = 1,
  maxpadj = 0.05,
  nTop = 6,
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

- TE_type:

  Character string specifying the TE type to analyze (example: "LTR",
  "Alu")

- specific_type:

  Character string for finer classification level: "TE_family" or
  "TE_name" (default: "TE_name")

- broad_type:

  Character string for broader classification level: "TE_class" or
  "TE_family". If NULL, automatically detected from TE_type (default:
  NULL)

- order:

  Character string specifying expression direction: "up", "down", or
  "all" (default: "up")

- minlfc:

  Numeric minimum log2 fold change threshold for significance (default:
  1)

- maxpadj:

  Numeric maximum adjusted p-value threshold for significance (default:
  0.05)

- nTop:

  Integer number of top dysregulated elements to display (default: 6)

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
  from `device`). If NULL, uses
  "violinPlot\_\<TE_type\>\_\<specific_type\>" (default: NULL)

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
violinPlotByTEtype(res.TEs, TE_type = "LINE",
                   specific_type = "TE_family", nTop = 10) 

violinPlotByTEtype(res.TEs, TE_type = "LINE",
                   specific_type = "TE_name", nTop = 10)
} # }
```
