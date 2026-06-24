# Annotate TEs with Genomic Regions and Generate Region Distribution Plots

Annotates transposable elements (TEs) with their genomic context
relative to protein-coding genes and generates stacked bar plots showing
the distribution of TEs across genomic regions (promoters, exons,
introns, etc.). This is a wrapper function that combines genomic
annotation with visualization.

## Usage

``` r
annotate_TE_regions(
  TE_results,
  gtf.genes.file,
  gtf.format = "gtf",
  output_folder = ".",
  device = "png",
  gene_biotypes = "protein_coding",
  plot.title = NULL,
  maxpadj = 0.05,
  minlfc = 1,
  width = 10,
  height = 7,
  minCounts = 10,
  TSSminus = -1000,
  TSSplus = 1000,
  downstream = 10000,
  is_ext_3UTR = FALSE,
  ext_3UTR_file = NULL
)
```

## Arguments

- TE_results:

  List object returned by
  [`TE_DEA`](https://guilleperis.github.io/TExpress/reference/TE_DEA.md),
  containing:

  res.TEs

  :   Data frame of DESeq2 results for TEs

  TE.count

  :   Data frame of normalized counts for TEs

  gene.count

  :   Data frame of normalized counts for genes

  metadata

  :   Sample metadata data frame

- gtf.genes.file:

  Character string. Full path to GTF file containing gene annotations.
  Should include protein-coding genes with transcript information. Must
  have columns: `type`, `transcript_biotype`, `gene_biotype`, `gene_id`,
  `gene_name`.

- gtf.format:

  Character string. Format of gene GTF file: "gtf" (default), "gff3· or
  "gff".

- output_folder:

  Character string. Path for output directory. A subdirectory
  "TEs_annotated" will be created inside this folder. Default is current
  directory (".").

- device:

  Character vector. File format(s) for output plots. Supported: "svg",
  "eps", "png", "tiff", "jpeg". Can specify multiple formats. Default is
  "png".

- gene_biotypes:

  Character vector. Features with exons. Default is "protein_coding",

- plot.title:

  Character string. Title for the region distribution plots. If NULL
  (default), a generic title will be used.

- maxpadj:

  Numeric. Adjusted p-value threshold for significance. Features with
  padj \< maxpadj are considered significant. Default is 0.05.

- minlfc:

  Numeric. Minimum absolute log2 fold change threshold for differential
  expression. Default is 1.

- width:

  Numeric. Plot width in inches. Default is 10.

- height:

  Numeric. Plot height in inches. Default is 7.

- minCounts:

  Numeric. Minimum total normalized counts threshold for considering a
  TE loci as expressed. Used for region distribution plots. Default is
  10.

- TSSminus:

  Integer. Number of bases upstream of TSS to define promoter region.
  Should be negative. Default is -1000 (1kb upstream).

- TSSplus:

  Integer. Number of bases downstream of TSS to define promoter region.
  Should be positive. Default is 1000 (1kb downstream).

- downstream:

  Integer. Maximum distance downstream of gene end to consider for
  "Downstream" annotation. Default is 10000 (10kb).

- is_ext_3UTR:

  Boolean. TRUE if extended 3'UTR analysis is to be performed. Defaults
  to FALSE. To be implemented

- ext_3UTR_file:

  Character string. Filename for 3'UTR analysis result, in gff format.
  To be implemented.

## Value

A list with four updated elements for downstream analysis:

- res.TEs:

  Data frame of DESeq2 results with added genomic annotation columns
  (annotation, geneChr, geneStart, geneEnd, geneId, gene_name, etc.)

- TE.count:

  Data frame of normalized TE counts with added genomic annotation
  columns

- gene.count:

  Data frame of normalized gene counts

- metadata:

  Original sample metadata (unchanged)

## Details

The function performs the following steps:

1.  Imports and filters gene GTF to protein-coding genes only

2.  Creates a TxDb object for ChIPseeker annotation

3.  Extracts transcript features for nearest gene assignment

4.  Annotates TEs with genomic context using
    [`TE_genomic_context`](https://guilleperis.github.io/TExpress/reference/TE_genomic_context.md)

5.  Saves annotated results to TSV file

6.  Generates region distribution plots using
    [`graphTEregion`](https://guilleperis.github.io/TExpress/reference/graphTEregion.md)

Only protein-coding genes are considered for annotation to focus on
biologically relevant gene-TE associations. Non-coding genes and
pseudogenes are excluded.

## Output Files

The function creates a subdirectory structure:


    output_folder/
      └── TEs_annotated/
          ├── DESeq2_TE_results_annotated.tsv
          └── [region distribution plots in specified format(s)]

## Examples

``` r
if (FALSE) { # \dontrun{
# Run differential expression analysis first
TE_results <- TE_DEA(
  metafile = "metadata.txt",
  folder = "counts",
  gtf.TE.file = "TEs.gtf"
)

# Annotate TEs and generate region plots
TE_results_annotated <- annotate_TE_regions(
  TE_results = TE_results,
  gtf.genes.file = "genes.gtf",
  output_folder = "results",
  device = c("jpeg", "png"),
  plot.title = "TE Distribution by Genomic Region",
  minCounts = 10
)

} # }
```
