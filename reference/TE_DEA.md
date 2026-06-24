# Perform Differential Expression Analysis for Transposable Elements

Main pipeline function for analyzing differential expression of
transposable elements (TEs) and genes using DESeq2. Reads count data,
performs statistical analysis, and generates visualization plots.

## Usage

``` r
TE_DEA(
  metafile,
  folder,
  output = ".",
  maxpadj = 0.05,
  minlfc = 1,
  gtf.TE.file,
  device = "png",
  plot.title = "",
  useCtrlGenes = FALSE,
  shrinklog2FC = FALSE,
  saveNorm = TRUE
)
```

## Arguments

- metafile:

  Character string. Full path to the metadata file. Should be a
  tab-separated file readable by
  [`read_metadata`](https://guilleperis.github.io/TExpress/reference/read_metadata.md).

- folder:

  Character string. Full path to the directory containing count files
  listed in the metadata.

- output:

  Character string. Path for output directory. Will be created if it
  doesn't exist. Default is current directory (".").

- maxpadj:

  Numeric. Adjusted p-value threshold for significance. Features with
  padj \< maxpadj are considered significant. Default is 0.05.

- minlfc:

  Numeric. Minimum absolute log2 fold change threshold for differential
  expression. Default is 1.

- gtf.TE.file:

  Character string. Full path to GTF file containing TE annotations.
  Must be the same GTF file used for upstream TE counting with TElocal
  from TEtranscript package.

- device:

  Character vector. File format(s) for output plots. Supported: "svg",
  "eps", "png", "tiff", "jpeg". Can specify multiple formats as a vector
  (e.g., `c("jpeg", "png")`). Default is "png".

- plot.title:

  Character string. Title for volcano and MA plots. If empty string
  (default), generic titles will be used.

- useCtrlGenes:

  Logical. If TRUE, uses only genes (not TEs) for estimating DESeq2 size
  factors. This can improve normalization when TEs have very different
  expression patterns than genes. Default is FALSE.

- shrinklog2FC:

  Logical. If TRUE, applies apeglm log2 fold change shrinkage in DESeq2
  for more accurate effect size estimates. Recommended for ranking and
  visualization. Default is FALSE.

- saveNorm:

  Logical. If TRUE, saves normalized count matrices for both genes and
  TEs. If FALSE, only saves DESeq2 results. Default is TRUE.

## Value

A list with four elements:

- res.TEs:

  Data frame of DESeq2 results for TEs, including genomic coordinates
  and TE annotations

- TE.count:

  Data frame of normalized counts for TEs across samples

- gene.count:

  Data frame of normalized counts for genes across samples

- metadata:

  Original metadata data frame used for the analysis

## Details

The pipeline performs the following steps:

1.  Reads sample metadata and count files

2.  Imports TE annotations from GTF file

3.  Adds genomic coordinates to TE features

4.  Filters non-standard chromosomes

5.  Performs DESeq2 differential expression analysis

6.  Extracts normalized counts

7.  Separates results for genes vs. TEs

8.  Saves results to separate directories

9.  Generates volcano and MA plots for both genes and TEs

The function distinguishes TEs from genes by the presence of a colon (:)
in the feature name, which is the standard format for TElocal TE loci
identifiers: `TE_element:TE_name:TE_family:TE_class` (e.g.,
"L1PA2_dup501:L1PA2:L1:LINE"). Gene identifiers contain no colon (e.g.,
"ENSMUSG00000000001").

## Output Structure

The function creates the following directory structure:


    output/
      ├── genes_DEA/
      │   ├── DESeq2_gene_results.tsv
      │   ├── gene_normalizedCounts.tsv (if saveNorm = TRUE)
      │   ├── volcanoPlot.<format>
      │   └── maPlot.<format>
      └── TEs_DEA/
          ├── DESeq2_TE_results.tsv
          ├── TE_normalizedCounts.tsv (if saveNorm = TRUE)
          ├── volcanoPlot.<format>
          └── maPlot.<format>

TE results include additional columns for genomic coordinates (seqnames,
start, end, strand, width) and TE class/family information.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with default parameters
results <- TE_DEA(
  metafile = "metadata.txt",
  folder = "counts",
  gtf.TE.file = "TEs.gtf"
)

# Full analysis with custom parameters
results <- TE_DEA(
  metafile = "samples.txt",
  folder = "count_data",
  output = "results/DE_analysis",
  maxpadj = 0.01,
  minlfc = 1.5,
  gtf.TE.file = "annotations/TE_annotation.gtf",
  device = c("tiff", "png"),
  plot.title = "KO vs WT",
  useCtrlGenes = TRUE,
  shrinklog2FC = TRUE,
  saveNorm = TRUE
)

# Access TE results
head(results$res.TEs)
dim(results$TE.count)

# Find significantly upregulated TEs
sig_up_TEs <- results$res.TEs[
  results$res.TEs$padj < 0.05 &
  results$res.TEs$log2FoldChange > 1,
]
} # }
```
