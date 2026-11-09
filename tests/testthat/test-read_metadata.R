test_that("read_metadata errors when datafile is missing", {
  expect_error(read_metadata(), "missing with no default")
})

test_that("read_metadata errors when the file does not exist", {
  expect_error(read_metadata("does_not_exist_12345.txt"), "not found")
})

test_that("read_metadata errors when column count is wrong", {
  bad <- make_metadata_df()[, 1:3]
  path <- write_metadata_file(bad)
  on.exit(unlink(path))
  expect_error(read_metadata(path), "columns")
})

test_that("read_metadata errors when Condition lacks Control and Treat", {
  bad <- make_metadata_df()
  bad$Condition <- c("A", "A", "B", "B")
  path <- write_metadata_file(bad)
  on.exit(unlink(path))
  expect_error(read_metadata(path), "Control.*Treat")
})

test_that("read_metadata returns a validated data frame on a valid file", {
  path <- write_metadata_file()
  on.exit(unlink(path))

  md <- read_metadata(path)

  expect_s3_class(md, "data.frame")
  expect_identical(colnames(md), c("File", "Sample", "Group", "Condition"))
  expect_equal(nrow(md), 4L)
  expect_setequal(unique(md$Condition), c("Control", "Treat"))
})

test_that("read_metadata validates the batch argument", {
  path <- write_metadata_file()
  on.exit(unlink(path))
  expect_error(read_metadata(path, batch = NA), "'batch' must be TRUE or FALSE")
  expect_error(read_metadata(path, batch = "yes"), "'batch' must be TRUE or FALSE")
})

test_that("read_metadata errors when batch = TRUE and Batch column is absent", {
  path <- write_metadata_file()
  on.exit(unlink(path))
  expect_error(read_metadata(path, batch = TRUE), "Batch")
})

test_that("read_metadata accepts a Batch column with batch = TRUE or FALSE", {
  df <- make_metadata_df()
  df$Batch <- c(1, 2, 1, 2)
  path <- write_metadata_file(df)
  on.exit(unlink(path))

  for (b in c(TRUE, FALSE)) {
    md <- read_metadata(path, batch = b)
    expect_identical(colnames(md),
                     c("File", "Sample", "Group", "Condition", "Batch"))
  }
})

test_that("read_metadata rejects a fifth column not named Batch", {
  df <- make_metadata_df()
  df$Other <- 1:4
  path <- write_metadata_file(df)
  on.exit(unlink(path))
  expect_error(read_metadata(path), "columns")
  expect_error(read_metadata(path, batch = TRUE), "Batch")
})

test_that("read_metadata errors when batch = TRUE and a batch lacks a condition", {
  df <- make_metadata_df()
  df$Batch <- c("b1", "b2", "b2", "b3")
  path <- write_metadata_file(df)
  on.exit(unlink(path))
  expect_error(read_metadata(path, batch = TRUE),
               "both conditions.*b1 \\(Control=1, Treat=0\\); b3 \\(Control=0, Treat=1\\)")
  expect_no_error(read_metadata(path, batch = FALSE))
})
