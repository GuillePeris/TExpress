# Generate Genomic Region and TE Class Distribution Plots

Creates stacked bar plots showing the distribution of transposable
elements across genomic regions (promoters, exons, introns, etc.) and TE
classes (LINE, SINE, LTR, DNA) for different sample groups and
differential expression categories.

## Usage

``` r
graphTEregion(
  TE_results,
  device = "png",
  output_folder = ".",
  width = 10,
  height = 7,
  minCounts = 10,
  minlfc = 1,
  maxpadj = 0.05,
  plot.title = NULL
)
```

## Arguments

- TE_results:

  List object returned by
  [`TE_DEA`](https://guilleperis.github.io/TExpress/reference/TE_DEA.md)
  and annotated by
  [`annotate_TE_regions`](https://guilleperis.github.io/TExpress/reference/annotate_TE_regions.md),
  containing:

  res.TEs

  :   Data frame of DESeq2 results with genomic annotations. Must
      include columns: `log2FoldChange`, `padj`, `annotation`,
      `TE_class`

  TE.count

  :   Data frame of normalized counts with genomic annotations. Columns
      should be sample names plus annotation columns

  gene.count

  :   Data frame of normalized counts for genes

  metadata

  :   Sample metadata with columns: `Sample`, `Condition`

- device:

  Character vector. File format(s) for output plots. Supported: "svg",
  "eps", "png", "tiff", "jpeg". Can specify multiple formats. Default is
  "png".

- output_folder:

  Character string. Directory where plots will be saved. Default is
  current directory (".").

- width:

  Numeric. Plot width in inches. Default is 10.

- height:

  Numeric. Plot height in inches. Default is 7.

- minCounts:

  Numeric. Minimum total normalized counts threshold for considering a
  TE "expressed" in a condition. TEs below this threshold across all
  samples in a condition are excluded. Default is 10.

- minlfc:

  Numeric. Minimum absolute log2 fold change for defining
  up/down-regulated TEs. Default is 1.

- maxpadj:

  Numeric. Maximum adjusted p-value for defining significant
  differential expression. Default is 0.05.

- plot.title:

  Character string. Title for the plots. If NULL, no title is added.
  Default is NULL.

## Value

Invisible NULL. Called for side effects (creating plot files).

## Details

The function generates two types of stacked bar plots:

1.  Genomic Region Distribution: Shows percentage of TEs in each genomic
    context (Promoter, 5' UTR, Exon, Intron, 3' UTR, Downstream,
    Intergenic)

2.  TE Class Distribution: Shows percentage of different TE classes
    (LINE, SINE, LTR, DNA, Other)

Four categories are plotted for each:

- Control: TEs expressed (total counts \> minCounts) in control samples

- Treatment: TEs expressed in treatment samples

- Down: Significantly down-regulated TEs (log2FC \< -minlfc, padj \<
  maxpadj)

- Up: Significantly up-regulated TEs (log2FC \> minlfc, padj \< maxpadj)

Sample group names are derived from the common prefix of sample names in
each condition (e.g., "WT1, WT2, WT3" becomes "WT").

## Output Files

Two plot files are created in the specified format(s):

- `stackBar_region.<format>`: Genomic region distribution

- `stackBar_TEclass.<format>`: TE class distribution

## Examples

``` r
if (FALSE) { # \dontrun{
# After running TE_DEA and TE_regionAnnot
TE_results <- TE_DEA(
  metafile = "metadata.txt",
  folder = "counts",
  gtf.TE.file = "TEs.gtf"
)

TE_results <- TE_regionAnnot(
  TE_results = TE_results,
  gtf.genes.file = "genes.gtf"
)

# Generate region distribution plots
graphTEregion(
  TE_results = TE_results,
  device = c("tiff", "png"),
  output_folder = "results/plots",
  width = 10,
  height = 7,
  minCounts = 10,
  minlfc = 1,
  maxpadj = 0.05,
  plot.title = "TE Distribution"
)
} # }
```
