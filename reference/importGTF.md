# Import and Standardize GTF/GFF Annotation File

Imports a GTF or GFF annotation file and standardizes it by keeping only
standard chromosomes and converting to UCSC chromosome naming style
(e.g., "chr1" instead of "1").

## Usage

``` r
importGTF(
  gtf.file,
  format = "gtf",
  keep.standard.chroms = TRUE,
  seqlevels.style = "UCSC"
)
```

## Arguments

- gtf.file:

  Character string. Path to a GTF or GFF format annotation file. Can be
  gzip-compressed (.gz).

- format:

  Character string. Format of the annotation file. Must be either "gtf"
  or "gff" (also accepts "gff3"). Default is "gtf".

- keep.standard.chroms:

  Logical. If TRUE (default), removes non-standard chromosomes (keeps
  only autosomes, sex chromosomes, and mitochondria). Scaffolds,
  alternative haplotypes, and patches are removed.

- seqlevels.style:

  Character string. Chromosome naming style to use. Options include
  "UCSC" (chr1, chr2, ...), "NCBI" (1, 2, ...), "Ensembl" (1, 2, ...).
  Default is "UCSC". Use NULL to skip style conversion.

## Value

A `GRanges` object containing the genomic annotations with standardized
chromosome names and filtered to standard chromosomes.

## Details

This function performs three operations:

1.  Imports the GTF/GFF file using `rtracklayer::import`

2.  Optionally filters to standard chromosomes (autosomes, sex
    chromosomes, mitochondria) using
    [`GenomeInfoDb::keepStandardChromosomes`](https://rdrr.io/pkg/GenomeInfoDb/man/seqlevels-wrappers.html)

3.  Optionally converts chromosome naming to specified style using
    [`GenomeInfoDb::seqlevelsStyle`](https://rdrr.io/pkg/GenomeInfoDb/man/seqlevelsStyle.html)

## Examples

``` r
if (FALSE) { # \dontrun{
# Import GTF file with default settings
gtf <- importGTF("annotations.gtf")

# Import GFF3 file
gff <- importGTF("annotations.gff3", format = "gff3")

# Import without filtering chromosomes
gtf_all <- importGTF("annotations.gtf", keep.standard.chroms = FALSE)

# Import and keep NCBI chromosome style (1, 2, 3...)
gtf_ncbi <- importGTF("annotations.gtf", seqlevels.style = "NCBI")

# Import for a specific species
gtf_mouse <- importGTF(
  "mouse_annotations.gtf",
  species = "Mus_musculus"
)

# Import compressed file
gtf <- importGTF("annotations.gtf.gz")
} # }
```
