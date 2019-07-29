###############################################################################

context("Tests for `CoxpresDbAccessor` objects in `coxpresdbr` package")

###############################################################################

test_data_genes <- as.character(c(
  2538791, 2539499, 2540594, 2541560, 2541907,
  2542210, 2542294, 2543492, 3361219, 3361512
))

test_gene_statistics <- tibble::tibble(
  gene_id = test_data_genes,
  p_value = c(0.1, 0.3, 0.5, 0.2, 1.0, 0.0, 0.0000001, 0.5, 0.5, 0.5),
  direction = sample(c(-1, 1), size = 10, replace = TRUE)
)

###############################################################################

test_that("CoxpresDbAccessor: class definition", {
  # Can't construct a CoxpresDbAccessor since it is a virtual class
})

###############################################################################

test_that("CoxpresDbArchiveAccessor: class definition", {
  archive_accessor <- methods::new("CoxpresDbArchiveAccessor")

  expect_is(
    object = archive_accessor,
    class = "CoxpresDbArchiveAccessor",
    info = "new `CoxpresDbArchiveAccessor` is a `CoxpresDbArchiveAccessor`"
  )

  expect_is(
    object = archive_accessor,
    class = "CoxpresDbAccessor",
    info = "new `CoxpresDbArchiveAccessor` is a subclass of `CoxpresDbAccessor`"
  )

  # TODO: object creation using the function `CoxpresDbAccessor`
})

###############################################################################

test_that("CoxpresDbDataframeAccessor: class definition", {
  # Note that CoxpresDbDataframeAccessor 'new' constructor is not exported
  # so the user does not have access to the following code.

  archive_accessor <- methods::new("CoxpresDbDataframeAccessor")

  expect_is(
    object = archive_accessor,
    class = "CoxpresDbDataframeAccessor",
    info = "new `CoxpresDbDataframeAccessor` is a `CoxpresDbDataframeAccessor`"
  )

  expect_is(
    object = archive_accessor,
    class = "CoxpresDbAccessor",
    info = paste(
      "new `CoxpresDbDataframeAccessor` is a subclass of `CoxpresDbAccessor`"
    )
  )

  df <- tibble::tibble(
    source_id = "123",
    target_id = "456",
    mutual_rank = 1.357
  )

  expect_is(
    object = CoxpresDbAccessor(df),
    class = "CoxpresDbDataframeAccessor",
    info = paste(
      "User can make a `CoxpresDbDataframeAccessor` by passing a `data.frame`",
      "to the `CoxpresDbAccessor` function"
    )
  )

  # If the user provides a dataframe for in-memory access to CoxpresDb, then
  # the columns `source_id`, `target_id`, `mutual_rank` should all be present
  expect_error(
    object = CoxpresDbAccessor(dplyr::rename(df, not_source = source_id)),
    info = paste(
      "If user makes a `CoxpresDbAccessor` from a `dataframe`, then",
      "`source_id` should be present in columns"
    )
  )
  expect_error(
    object = CoxpresDbAccessor(dplyr::rename(df, not_target = target_id)),
    info = paste(
      "If user makes a `CoxpresDbAccessor` from a `dataframe`, then",
      "`target_id` should be present in columns"
    )
  )
  expect_error(
    object = CoxpresDbAccessor(dplyr::rename(df, not_rank = mutual_rank)),
    info = paste(
      "If user makes a `CoxpresDbAccessor` from a `dataframe`, then",
      "`mutual_rank` should be present in columns"
    )
  )
})
