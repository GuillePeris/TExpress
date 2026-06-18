# ---- filter_TE_annot -----------------------------------------------------

test_that("filter_TE_annot errors on missing arg / non-existent file", {
  expect_error(filter_TE_annot(), "missing with no default")
  expect_error(filter_TE_annot("no_such_file.gtf"), "not found")
})

test_that("filter_TE_annot errors when class_id is absent", {
  skip_if_not_installed("rtracklayer")
  path <- write_incomplete_gtf()
  on.exit(unlink(path))
  expect_error(filter_TE_annot(path), "missing required column")
})

test_that("filter_TE_annot keeps only requested TE classes", {
  skip_if_not_installed("rtracklayer")
  path <- write_small_gtf()
  on.exit(unlink(path), add = TRUE)
  
  # Fixture class_id values are LINE, SINE, DNA; keep only LINE.
  filter_TE_annot(path, TE_class = "LINE", suffix = "filtered")
  
  out <- sub("\\.gtf$", "_filtered.gtf", path)
  on.exit(unlink(out), add = TRUE)
  expect_true(file.exists(out))
  expect_equal(count_gtf_features(out), 1L)
})
