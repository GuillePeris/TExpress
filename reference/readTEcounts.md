# Read and Merge TE Count Data by Sample

Reads transposable element (TE) count files for multiple samples and
merges them into a single count matrix. Count files from TElocal
(TEtranscripts) should be tab-separated with two columns: feature names
and counts.

## Usage

``` r
readTEcounts(metadata, folder)
```

## Arguments

- metadata:

  Data frame with sample metadata, as returned by
  [`read_metadata`](https://guilleperis.github.io/TExpress/reference/read_metadata.md).
  Must contain columns: File, Sample, Group, and Condition.

- folder:

  Character string. Path to the directory containing the count files
  listed in the first column of `metadata`. Can be absolute or relative
  path.

## Value

A data frame with features as row names and samples as columns. Column
names correspond to the Sample column in `metadata`. Column order
matches the row order in `metadata`.

## Details

Each count file must:

- Be tab-separated

- Have a header row

- Contain exactly 2 columns: feature name and count value

- Use the same feature names across all samples (order can vary)

The function attempts to efficiently merge count data:

- If all files have identical features in the same order, uses fast
  column binding

- If feature order differs, performs full outer join to align features

- Missing features in some samples are filled with NA (should not occur
  in well-formatted data)

## Examples

``` r
if (FALSE) { # \dontrun{
# Read metadata
metadata <- read_metadata("metadata.txt")

# Read count files from directory
countData <- readTEcounts(metadata, folder = "counts")

# Verify dimensions
dim(countData)
head(countData)
} # }
```
