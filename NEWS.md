# TExpress 1.0

First public release.

* Differential expression analysis of transposable elements from TElocal count
  tables via `TE_DEA()` (DESeq2-based, Treat-vs-Control contrast).
* Genomic-context annotation of TEs relative to a gene GTF with
  `annotate_TE_regions()` (Promoter / 5' UTR / Exon / Intron / 3' UTR /
  Downstream / Intergenic).
* Transcriptional classification of TEs as self-expressed or gene-dependent
  with `classify_TE_transcription()`.
* GTF utilities: `filter_GTF()`, `filter_TE_annot()`, and `gtfBothStrands()`.
* Publication-ready plots: volcano/MA plots, violin plots, and stacked bar
  plots of genomic-region and TE-class distributions.
* `downloadTestData()` to fetch a small chr22 example dataset for testing and
  the package examples.
