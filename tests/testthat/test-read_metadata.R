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
