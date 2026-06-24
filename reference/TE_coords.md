# Add Genomic Coordinates to TE Count Data

Extracts transposable element (TE) information from count data rownames
and adds genomic coordinates by joining with TE annotation data.

## Usage

``` r
TE_coords(countData.TEs, TE_annot.df)
```

## Arguments

- countData.TEs:

  Data frame with TE count data. Rownames must follow the format:
  "TE_element:TE_name:TE_family:TE_class" (colon-separated). Each column
  should contain count data for a sample.

- TE_annot.df:

  Data frame with TE annotations, typically derived from a GTF file.
  Must contain columns: seqnames, start, end, strand, and transcript_id.
  The transcript_id column is used to match with the first field in
  countData.TEs rownames.

## Value

Data frame containing the original count data plus additional columns:

- TE_element: Transcript ID extracted from rownames

- TE_name: TE name extracted from rownames

- TE_family: TE family extracted from rownames

- TE_class: TE class extracted from rownames

- seqnames: Chromosome/sequence name from annotation

- start: Start position from annotation

- end: End position from annotation

- strand: Strand information from annotation

Original rownames are preserved.

@importFrom rlang .data

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming countData.TEs has rownames like:
# "L1PA2_dup501:L1PA2:L1:LINE"
# "AluY_dup68349:AluY:Alu:SINE"

# And TE_annot.df has columns: transcript_id, seqnames, start, end, strand

countData_with_coords <- TE_coords(countData.TEs, TE_annot.df)

# Check results
head(countData_with_coords)
} # }
```
