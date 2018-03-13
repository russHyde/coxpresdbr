###############################################################################

context("tests for class definitions in `coxpresdbr` package")

###############################################################################

test_that("CoxpresDbImporter: class definition", {
  expect_is(
    object = new("CoxpresDbImporter"),
    class = "CoxpresDbImporter",
    info = "new `CoxpresDbImporter` object has appropriate class-name"
  )
})

###############################################################################

test_that("CoxpresDbPartners: class definition", {
  expect_is(
    object = new("CoxpresDbPartners"),
    class = "CoxpresDbPartners",
    info = "new `CoxpresDbPartners` object has appropriate class-name"
  )

  expect_error(
    object = new("CoxpresDbPartners", cluster_graph = "Not an igraph"),
    info = paste(
      "If not NULL, `cluster_graph` should be an `igraph` in",
      "`CoxpresDbPartners`"
    )
  )
})

###############################################################################
