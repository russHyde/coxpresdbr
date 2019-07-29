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
})

###############################################################################
