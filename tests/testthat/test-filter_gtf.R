test_that("filter_GTF errors on missing arg / non-existent file", {
  expect_error(filter_GTF(), "missing with no default")
  expect_error(filter_GTF("no_such_file.gtf"), "not found")
})

test_that("filter_GTF errors when required columns are absent", {
  skip_if_not_installed("rtracklayer")
  path <- write_incomplete_gtf()
  on.exit(unlink(path))
  expect_error(filter_GTF(path), "missing required column")
})

test_that("filter_GTF keeps only requested biotypes and names output", {
  skip_if_not_installed("rtracklayer")
  path <- write_small_gtf()
  on.exit(unlink(path), add = TRUE)

  filter_GTF(path, features = "protein_coding", suffix = "filtered")

  out <- sub("\\.gtf$", "_filtered.gtf", path)
  on.exit(unlink(out), add = TRUE)
  expect_true(file.exists(out))
  # Two of the three fixture features are protein_coding.
  expect_equal(count_gtf_features(out), 2L)
})

