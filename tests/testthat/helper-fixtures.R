# Shared fixture builders for the TExpress test suite.
# testthat auto-sources helper-*.R before running the tests.

# A valid metadata data frame: 4 columns, Control + Treat conditions.
make_metadata_df <- function() {
  data.frame(
    File      = c("WT1.cntTable", "WT2.cntTable", "KO1.cntTable", "KO2.cntTable"),
    Sample    = c("WT1", "WT2", "KO1", "KO2"),
    Group     = c("WT", "WT", "KO", "KO"),
    Condition = c("Control", "Control", "Treat", "Treat"),
    stringsAsFactors = FALSE
  )
}

# Write a tab-separated metadata file; returns the path.
write_metadata_file <- function(df = make_metadata_df(),
                                path = tempfile(fileext = ".txt")) {
  utils::write.table(df, path, sep = "\t", quote = FALSE,
                     row.names = FALSE, col.names = TRUE)
  path
}

# Write one 2-column TElocal-style count file per sample into `dir`.
# `features` is a character vector of feature names.
# `sample_counts` is a named list (names = sample) of integer count vectors,
# each the same length and order as `features`.
# When `shuffle = TRUE`, each file's rows are written in a different random
# order to exercise readTEcounts()'s row-name merge (slow) path.
write_count_files <- function(features, sample_counts,
                              dir = tempfile(), shuffle = FALSE) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  samples <- names(sample_counts)
  for (i in seq_along(samples)) {
    s <- samples[i]
    idx <- seq_along(features)
    if (shuffle) {
      # Deterministic per-file reordering without relying on set.seed.
      idx <- (idx + i) %% length(features) + 1L
      idx <- unique(c(idx, seq_along(features)))[seq_along(features)]
    }
    df <- data.frame(feature = features[idx],
                     count   = sample_counts[[s]][idx],
                     stringsAsFactors = FALSE)
    colnames(df) <- c("feature/sample", s)
    utils::write.table(df, file.path(dir, paste0(s, ".cntTable")),
                       sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = TRUE)
  }
  dir
}

# A count matrix whose rownames follow transcript_id:TE_name:TE_family:TE_class.
make_counts_TE_df <- function() {
  df <- data.frame(
    WT1 = c(10, 20, 30),
    KO1 = c(15, 25, 35),
    stringsAsFactors = FALSE
  )
  rownames(df) <- c(
    "L1PA2_dup501:L1PA2:L1:LINE",
    "AluY_dup68349:AluY:Alu:SINE",
    "MIRb_dup12:MIRb:MIR:SINE"
  )
  df
}

# Annotation data frame keyed by transcript_id (matches TE_element field).
make_TE_annot_df <- function() {
  data.frame(
    seqnames      = c("chr1", "chr2", "chr3"),
    start         = c(100L, 200L, 300L),
    end           = c(150L, 250L, 350L),
    strand        = c("+", "-", "+"),
    transcript_id = c("L1PA2_dup501", "AluY_dup68349", "MIRb_dup12"),
    stringsAsFactors = FALSE
  )
}

# A DEA-results-like data frame for the violinplot / region helpers.
make_res_TEs_df <- function() {
  data.frame(
    log2FoldChange = c( 2.5, -3.0,  0.1,  1.8, -2.2,  0.0),
    padj           = c(0.001, 0.002, 0.5, 0.01, 0.03, 0.9),
    TE_element     = paste0("dup", 1:6),
    TE_name        = c("L1PA2", "AluY", "MIRb", "L1PA3", "AluSx", "MIRb"),
    TE_family      = c("L1", "Alu", "MIR", "L1", "Alu", "MIR"),
    TE_class       = c("LINE", "SINE", "SINE", "LINE", "SINE", "SINE"),
    stringsAsFactors = FALSE
  )
}

# Build a tiny GRanges and export it as a GTF; returns the path.
# Includes both strands and the mcols the filter functions require.
write_small_gtf <- function(path = tempfile(fileext = ".gtf")) {
  gr <- GenomicRanges::GRanges(
    seqnames = c("1", "1", "2"),
    ranges   = IRanges::IRanges(start = c(100, 500, 900),
                                end   = c(200, 600, 1000)),
    strand   = c("+", "-", "+")
  )
  GenomicRanges::mcols(gr) <- S4Vectors::DataFrame(
    type               = c("gene", "gene", "transcript"),
    gene_id            = c("g1", "g2", "g3"),
    transcript_id      = c("t1", "t2", "t3"),
    gene_biotype       = c("protein_coding", "lncRNA", "protein_coding"),
    transcript_biotype = c("protein_coding", "lncRNA", "protein_coding"),
    class_id           = c("LINE", "SINE", "DNA")
  )
  rtracklayer::export(gr, path, format = "gtf")
  path
}

# Export a GRanges missing the biotype/class mcols the filter functions require
# (and lacking transcript_id, to exercise gtfBothStrands' warning path).
write_incomplete_gtf <- function(path = tempfile(fileext = ".gtf")) {
  gr <- GenomicRanges::GRanges(
    seqnames = c("1", "2"),
    ranges   = IRanges::IRanges(start = c(100, 500), end = c(200, 600)),
    strand   = c("+", "-")
  )
  GenomicRanges::mcols(gr) <- S4Vectors::DataFrame(
    type    = c("gene", "gene"),
    gene_id = c("g1", "g2")
  )
  rtracklayer::export(gr, path, format = "gtf")
  path
}

# Number of features in a GTF on disk.
count_gtf_features <- function(path) {
  length(rtracklayer::import(path, format = "gtf"))
}

# Creates a small counts list with genes and TEs
make_counts_list <- function() {
  list(
    features = c("geneA", "geneB", "L1PA2:L1PA2:L1:LINE", "AluY:AluY:Alu:SINE"),
    counts   = list(
      WT1 = c(1L, 2L, 3L, 4L),
      WT2 = c(5L, 6L, 7L, 8L),
      KO1 = c(9L, 10L, 11L, 12L),
      KO2 = c(13L, 14L, 15L, 16L)
    )
  )
}