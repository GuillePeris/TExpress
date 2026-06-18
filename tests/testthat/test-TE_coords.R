test_that("TE_coords errors on missing arguments", {
  expect_error(TE_coords(), "missing with no default")
  expect_error(TE_coords(make_counts_TE_df()), "missing with no default")
})

test_that("TE_coords errors when TE_annot.df lacks required columns", {
  annot <- make_TE_annot_df()[, c("seqnames", "start", "transcript_id")]
  expect_error(
    TE_coords(make_counts_TE_df(), annot),
    "missing required column"
  )
})

test_that("TE_coords errors when rownames have no colon", {
  counts <- make_counts_TE_df()
  rownames(counts) <- c("a", "b", "c")
  expect_error(TE_coords(counts, make_TE_annot_df()), "colon-separated")
})

test_that("TE_coords splits rownames and joins coordinates", {
  counts <- make_counts_TE_df()
  annot <- make_TE_annot_df()

  res <- TE_coords(counts, annot)

  expect_true(all(c("TE_element", "TE_name", "TE_family", "TE_class",
                    "seqnames", "start", "end", "strand") %in% colnames(res)))
  # Rownames preserved.
  expect_identical(rownames(res), rownames(counts))
  # First feature parsed and joined correctly.
  expect_equal(res["L1PA2_dup501:L1PA2:L1:LINE", "TE_name"], "L1PA2")
  expect_equal(res["L1PA2_dup501:L1PA2:L1:LINE", "TE_class"], "LINE")
  expect_equal(res["L1PA2_dup501:L1PA2:L1:LINE", "seqnames"], "chr1")
  expect_equal(res["L1PA2_dup501:L1PA2:L1:LINE", "start"], 100L)
})

test_that("addTEColumns splits rownames and joins on TE_element", {
  counts <- make_counts_TE_df()
  annot <- make_TE_annot_df()
  # addTEColumns joins on a TE_element column rather than transcript_id.
  colnames(annot)[colnames(annot) == "transcript_id"] <- "TE_element"

  res <- TExpress:::addTEColumns(counts, annot)

  expect_identical(rownames(res), rownames(counts))
  expect_equal(res["AluY_dup68349:AluY:Alu:SINE", "TE_family"], "Alu")
  expect_equal(res["AluY_dup68349:AluY:Alu:SINE", "strand"], "-")
})
