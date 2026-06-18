test_that("importGTF errors on missing argument", {
  expect_error(importGTF(), "missing with no default")
})

test_that("importGTF errors when the file does not exist", {
  expect_error(importGTF("no_such_file_54321.gtf"), "not found")
})

test_that("importGTF imports a GRanges and applies UCSC seqlevels style", {
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomeInfoDb")

  path <- write_small_gtf()
  on.exit(unlink(path))

  gtf <- importGTF(path)

  expect_s4_class(gtf, "GRanges")
  # seqnames written as "1"/"2" should become "chr1"/"chr2" under UCSC style.
  expect_true(all(grepl("^chr", as.character(GenomicRanges::seqnames(gtf)))))
})

