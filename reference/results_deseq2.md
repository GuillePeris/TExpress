# Get Results After DESeq2 Analysis

Extracts differential expression results from a DESeqDataSet object,
with optional log2 fold change shrinkage for improved effect size
estimates.

## Usage

``` r
results_deseq2(dds, shrinklog2FC = FALSE)
```

## Arguments

- dds:

  A DESeqDataSet object as returned by
  [`call_deseq2`](https://guilleperis.github.io/TExpress/reference/call_deseq2.md)
  or [`DESeq2::DESeq()`](https://rdrr.io/pkg/DESeq2/man/DESeq.html).
  Must have completed the DESeq2 analysis pipeline.

- shrinklog2FC:

  Logical. If TRUE, applies log2 fold change shrinkage using the apeglm
  method to produce more accurate effect size estimates. If FALSE,
  returns standard DESeq2 results. Default is FALSE. Note: Shrinkage
  requires the `apeglm` package to be installed.

## Value

A data frame with DESeq2 results with rows sorted by adjusted p-value
(most significant first). Row names correspond to feature names from the
count data.

## Details

This function extracts differential expression results comparing the
"Treat" condition against the "Control" condition (reference level).

When `shrinklog2FC = TRUE`, log2 fold changes are shrunk using the
apeglm method
([`DESeq2::lfcShrink`](https://rdrr.io/pkg/DESeq2/man/lfcShrink.html))

Rows with NA values (typically low-count genes filtered by DESeq2's
independent filtering) are removed, and results are sorted by adjusted
p-value.

## Examples

``` r
if (FALSE) { # \dontrun{
# Run DESeq2 analysis
dds <- call_deseq2(countData, metadata, useCtrlGenes = FALSE)

# Get results without shrinkage
res <- results_deseq2(dds, shrinklog2FC = FALSE)

# Get results with log2FC shrinkage (requires apeglm package)
res_shrink <- results_deseq2(dds, shrinklog2FC = TRUE)

# View top results
head(res)
} # }
```
