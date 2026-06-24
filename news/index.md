# Changelog

## TExpress 1.0

First public release.

- Differential expression analysis of transposable elements from TElocal
  count tables via
  [`TE_DEA()`](https://guilleperis.github.io/TExpress/reference/TE_DEA.md)
  (DESeq2-based, Treat-vs-Control contrast).
- Genomic-context annotation of TEs relative to a gene GTF with
  [`annotate_TE_regions()`](https://guilleperis.github.io/TExpress/reference/annotate_TE_regions.md)
  (Promoter / 5’ UTR / Exon / Intron / 3’ UTR / Downstream /
  Intergenic).
- Transcriptional classification of TEs as self-expressed or
  gene-dependent with
  [`classify_TE_transcription()`](https://guilleperis.github.io/TExpress/reference/classify_TE_transcription.md).
- GTF utilities: `filter_GTF()`, `filter_TE_annot()`, and
  [`gtfBothStrands()`](https://guilleperis.github.io/TExpress/reference/gtfBothStrands.md).
- Publication-ready plots: volcano/MA plots, violin plots, and stacked
  bar plots of genomic-region and TE-class distributions.
- [`downloadTestData()`](https://guilleperis.github.io/TExpress/reference/downloadTestData.md)
  to fetch a small chr22 example dataset for testing and the package
  examples.
