# Perform Differential Expression Analysis Using DESeq2

Performs differential expression analysis on count data using DESeq2.
Optionally uses only gene features (non-TE) for size factor estimation.

## Usage

``` r
call_deseq2(countData, metadata, useCtrlGenes)
```

## Arguments

- countData:

  Data frame with count matrix as returned by
  [`readTEcounts`](https://guilleperis.github.io/TExpress/reference/readTEcounts.md).
  Features in rows, samples in columns. Row names identify features; TE
  features are expected to contain ":" in their names.

- metadata:

  Data frame with sample metadata as returned by
  [`read_metadata`](https://guilleperis.github.io/TExpress/reference/read_metadata.md).
  Must contain columns: File, Sample, Group, and Condition. Sample names
  must match column names in countData.

- useCtrlGenes:

  Logical. If TRUE, only gene features (those without ":" in their
  names) are used to estimate size factors. If FALSE, all features are
  used. Default is FALSE.

## Value

A DESeqDataSet object after running the DESeq2 analysis pipeline.

## Details

The function creates a DESeqDataSet with a design formula `~ condition`,
where condition is derived from the Condition column in metadata
(Control vs Treat).

When `useCtrlGenes = TRUE`, size factors are estimated using only gene
features (rownames without ":"), which can be more appropriate when
analyzing datasets containing both genes and transposable elements. This
approach assumes that TEs may have more variable expression and could
bias normalization.

The function assumes TE features contain ":" in their rownames (e.g.,
"transcript_id:TE_name:TE_family:TE_class").

## Examples

``` r
if (FALSE) { # \dontrun{
# Read data
metadata <- read_metadata("metadata.txt")
countData <- readTEcounts(metadata, "counts")

# Run DESeq2 using all features for normalization
dds <- call_deseq2(countData, metadata, useCtrlGenes = FALSE)

# Run DESeq2 using only genes for normalization
dds_ctrl <- call_deseq2(countData, metadata, useCtrlGenes = TRUE)

# Extract results
res <- results_deseq2(dds)
} # }
```
