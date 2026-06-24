# Get Normalized Counts from DESeq2 Object

Extracts normalized count data from a DESeqDataSet object. Counts are
normalized by size factors to account for differences in sequencing
depth between samples.

## Usage

``` r
norm_counts(dds)
```

## Arguments

- dds:

  A DESeqDataSet object as returned by
  [`call_deseq2`](https://guilleperis.github.io/TExpress/reference/call_deseq2.md).
  Size factors must have been calculated (automatically done by
  [`DESeq2::DESeq()`](https://rdrr.io/pkg/DESeq2/man/DESeq.html)).

## Value

A numeric matrix with normalized counts. Features in rows, samples in
columns. Row and column names are preserved from the DESeqDataSet.

## Details

This is a convenience wrapper around
`DESeq2::counts(dds, normalized = TRUE)`.

Normalized counts are calculated by dividing raw counts by size factors.
These counts are appropriate for visualization and clustering, but
should NOT be used as input for differential expression analysis (use
raw counts for that).

## Examples

``` r
if (FALSE) { # \dontrun{
# Run DESeq2 analysis
dds <- call_deseq2(countData, metadata)

# Extract normalized counts
norm_counts_matrix <- norm_counts(dds)

# Convert to data frame if needed
norm_counts_df <- as.data.frame(norm_counts_matrix)
} # }
```
