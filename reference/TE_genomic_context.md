# Annotate TEs with Genomic Context Relative to Genes

Annotates transposable element (TE) loci with their genomic context
relative to protein-coding genes. Uses ChIPseeker to determine if TEs
overlap promoters, exons, introns, UTRs, or intergenic regions. Assigns
the nearest gene to intergenic TEs.

## Usage

``` r
TE_genomic_context(
  df.TEs,
  gene.TxDb,
  gene_names,
  transcript.gr,
  TSSminus,
  TSSplus,
  downstream,
  is_ext_3UTR = FALSE,
  ext_3UTR_file = NULL
)
```

## Arguments

- df.TEs:

  Data frame containing TE loci with genomic coordinates. Must include
  columns: `seqnames`, `start`, `end`, and `strand`. Typically the
  output from differential expression analysis.

- gene.TxDb:

  A TxDb object containing gene/transcript annotations. Can be created
  from GTF files using `GenomicFeatures::makeTxDbFromGFF` or loaded from
  Bioconductor annotation packages.

- gene_names:

  Data frame mapping gene IDs to gene names. Must contain columns:
  `gene_id` (matching IDs in TxDb) and `gene_name` (human-readable gene
  symbols). Used to add gene names to final output.

- transcript.gr:

  GRanges object containing transcript annotations. Must include a
  `gene_id` column in metadata. Used for assigning the nearest gene to
  intergenic TEs. Typically obtained by filtering a gene GTF for entries
  where `type == "transcript"`.

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

Data frame combining original TE information with genomic annotation
columns:

- annotation:

  Simplified genomic context (Promoter, 5' UTR, Exon, Intron, 3' UTR,
  Downstream, Intergenic)

- geneChr, geneStart, geneEnd, geneLength, geneStrand:

  Genomic coordinates of the associated gene

- geneId:

  Gene identifier from the TxDb

- gene_name:

  Human-readable gene symbol (if available in gene_names)

- distanceToTSS:

  Distance to transcription start site (from ChIPseeker)

## Details

The function performs genomic annotation in several steps:

1.  Converts TE data frame to GRanges object

2.  Annotates TEs using ChIPseeker with priority to coding regions

3.  Simplifies annotation categories (e.g., "Exon (uc057wby.1/267097,
    exon 2 of 7)" becomes "Exon")

4.  For intergenic TEs, assigns the nearest transcript from
    `transcript.gr`

5.  Adds gene names from the provided mapping table

Annotation priorities (from highest to lowest): 5' UTR → 3' UTR → Exon →
Promoter → Intron → Downstream → Intergenic

Note: ChIPseeker may internally convert UCSC chromosome naming (chr1) to
NCBI style (1). This function converts them back to UCSC format for
consistency.

## Examples

``` r
if (FALSE) { # \dontrun{
# Annotate TEs (assuming TE_results from TE_DEA)
TE_annotated <- TE_genomic_context(
  df.TEs = TE_results$res.TEs,
  gene.TxDb = txdb,
  gene_names = gene_names,
  transcript.gr = transcript.gr,
  TSSminus = -3000,
  TSSplus = 3000
)
} # }
```
