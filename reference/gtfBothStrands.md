# Create Antisense TE GTF with Both Strands

Generates a GTF file containing both sense and antisense annotations by
duplicating all features and inverting the strand for antisense copies.
Useful for analyzing antisense transcription or bidirectional promoters.

## Usage

``` r
gtfBothStrands(gtf.file, output.file = NULL, add_suffix = "_AS")
```

## Arguments

- gtf.file:

  Character string. Path to input TE GTF file. The file will be read,
  duplicated with strand inversion, and exported to a new file with
  "\_AS" suffix.

- output.file:

  Character string. Optional path for output TE GTF file. If NULL
  (default), creates a file in the same directory as input with "\_AS"
  added before the extension (e.g., "hg38_rmsk_TE.gtf" becomes
  "hg38_rmsk_TE_AS.gtf").

- add_suffix:

  Character string. Suffix to add to gene_id and transcript_id for
  antisense features. Default is "\_AS".

## Value

Character string (invisibly). Path to the created output GTF file. Also
prints a message with file location and feature counts.

## Details

The function performs the following steps:

1.  Imports the input TE GTF file

2.  Creates antisense copies by inverting strand information

3.  Appends suffix to gene_id and transcript_id in antisense copies

4.  Merges sense and antisense features

5.  Sorts by genomic coordinates

6.  Exports to new TE GTF file

For each feature on the "+" strand, an antisense copy on the "-" strand
is created, and vice versa. Gene and transcript IDs are modified to
distinguish antisense features (e.g., "ENSG00000000001" becomes
"ENSG00000000001_AS").

## Examples

``` r
if (FALSE) { # \dontrun{
# Create antisense GTF with default naming
gtfBothStrands("genes.gtf")
# Output: genes_AS.gtf

# Specify custom output file
gtfBothStrands(
  gtf.file = "hg38_rmsk_TE.gtf",
  output.file = "hg38_rmsk_TE_AS.gtf"
)

# Custom suffix for antisense features
gtfBothStrands(
  gtf.file = "hg38_rmsk_TE.gtf",
  add_suffix = "_antisense"
)
} # }
```
