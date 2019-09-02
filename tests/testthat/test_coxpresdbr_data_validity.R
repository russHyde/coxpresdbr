###############################################################################

context("Tests for data-structure validity in `coxpresdbr` package")

###############################################################################

test_that(".is_gene_statistics_df", {
  # should be a non-null dataframe
  expect_error(
    object = .is_gene_statistics_df(),
    info = "x should be defined in .is_gene_statistics_df"
  )
  expect_false(
    object = .is_gene_statistics_df("Not a data-frame"),
    info = "x should be a data-frame in .is_gene_statistics_df"
  )
  expect_false(
    object = .is_gene_statistics_df(tibble::tibble()),
    info = "data frame x should have some contents in .is_gene_statistics"
  )
  expect_false(
    object = .is_gene_statistics_df(tibble::tibble(
      gene_id = character(0), p_value = numeric(0), direction = numeric(0)
    ))
  )
  # should have a gene_id column
  expect_false(
    object = .is_gene_statistics_df(tibble::tibble(
      GENE_ID = "abc", p_value = 1, direction = -1
    )),
    info = "should have a gene_id column"
  )
  # should have a p_value column
  expect_false(
    object = .is_gene_statistics_df(tibble::tibble(
      gene_id = "abc", P_VALUE = 0.5, direction = 1
    )),
    info = "should have a p_value column"
  )
  # should have a direction column
  expect_false(
    object = .is_gene_statistics_df(tibble::tibble(
      gene_id = "abc", p_value = 0.5, DIRECTION = 1
    )),
    info = "should have a direction column"
  )
  # should be at most one copy of each gene in the gene_id column
  expect_false(
    object = .is_gene_statistics_df(tibble::tibble(
      gene_id = rep("abc", 2), p_value = 0:1, direction = rep(1, 2)
    )),
    info = "should be no repeats of a gene_id"
  )
})
