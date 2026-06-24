# Read Metadata File with Sample Information

Reads and validates a tab-separated metadata file containing sample
information for differential expression analysis. The file must contain
specific columns and follow formatting requirements for DESeq2
compatibility.

## Usage

``` r
read_metadata(datafile)
```

## Arguments

- datafile:

  Character string. Path to a tab-separated file (.txt, .tsv, or .csv)
  containing metadata. Must have exactly 4 columns with headers: File,
  Sample, Group, and Condition.

## Value

A data frame with four columns (File, Sample, Group, Condition) and
validated metadata. The Condition column is converted to a factor with
"Control" as the reference level.

## Details

The metadata file must follow these requirements:

- Exactly 4 columns with headers: File, Sample, Group, Condition

- At least one row of data (beyond the header)

- The "Condition" column must contain exactly two unique values:
  "Control" and "Treat"

- No duplicate values in the "Sample" or "File" columns

- Tab-separated format (other delimiters not supported)

The "Control" and "Treat" labels in the Condition column are required
for proper DESeq2 analysis, where "Control" will be set as the reference
level.

## Examples

``` r
if (FALSE) { # \dontrun{
# Read metadata from file
metadata <- read_metadata("path/to/metadata.txt")

# Example valid file structure:
# File            Sample  Group   Condition
# WT1.cntTable    WT1     WT      Control
# WT2.cntTable    WT2     WT      Control
# KO1.cntTable    KO1     KO      Treat
# KO2.cntTable    KO2     KO      Treat
} # }
```
