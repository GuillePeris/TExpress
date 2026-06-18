test_that("readTEcounts errors on missing arguments", {
  expect_error(readTEcounts(), "missing with no default")
  expect_error(readTEcounts(make_metadata_df()), "missing with no default")
})

test_that("readTEcounts errors when the folder does not exist", {
  expect_error(
    readTEcounts(make_metadata_df(), "no_such_folder_98765"),
    "does not exist"
  )
})

test_that("readTEcounts errors when a listed count file is absent", {
  spec <- make_counts_list()
  dir <- write_count_files(spec$features, spec$counts)
  on.exit(unlink(dir, recursive = TRUE))
  unlink(file.path(dir, "KO2.cntTable"))

  expect_error(readTEcounts(make_metadata_df(), dir), "Missing count file")
})

test_that("readTEcounts merges identically-ordered files (fast path)", {
  spec <- make_counts_list()
  dir <- write_count_files(spec$features, spec$counts)
  on.exit(unlink(dir, recursive = TRUE))

  md <- make_metadata_df()
  res <- readTEcounts(md, dir)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), length(spec$features))
  expect_identical(colnames(res), md$Sample)
  expect_setequal(rownames(res), spec$features)
  # Spot-check a value: geneA for WT1 should be 1.
  expect_equal(res["geneA", "WT1"], 1)
})

test_that("readTEcounts aligns features when row order differs (slow path)", {
  spec <- make_counts_list()
  dir <- write_count_files(spec$features, spec$counts, shuffle = TRUE)
  on.exit(unlink(dir, recursive = TRUE))

  md <- make_metadata_df()
  res <- readTEcounts(md, dir)

  expect_equal(nrow(res), length(spec$features))
  expect_setequal(rownames(res), spec$features)
  # Regardless of on-disk order, each feature keeps its per-sample value.
  expect_equal(res["geneB", "KO1"], 10)
  expect_equal(res["AluY:AluY:Alu:SINE", "WT2"], 8)
})
