###############################################################################

context("Tests for class definitions in `coxpresdbr` package")

###############################################################################

test_data_genes <- as.character(c(
  2538791, 2539499, 2540594, 2541560, 2541907,
  2542210, 2542294, 2543492, 3361219, 3361512
))

test_gene_statistics <- tibble::data_frame(
  gene_id = test_data_genes,
  p_value = runif(10),
  direction = sample(c(-1, 1), size = 10, replace = TRUE)
)

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

test_that("CoxpresDbPartners: field matches", {
  expect_equal(
    object = new(
      "CoxpresDbPartners",
      gene_statistics = test_gene_statistics
    )@gene_statistics,
    expected = test_gene_statistics,
    info = paste(
      "The input to the gene_statistics field of CoxpresDbPartners",
      "should match it's output"
    )
  )

  expect_equal(
    object = get_gene_universe(
      new(
        "CoxpresDbPartners",
        gene_statistics = test_gene_statistics
      )
    ),
    expected = test_gene_statistics$gene_id,
    info = paste(
      "A CoxpresDbPartners made from just a gene_statistics",
      "should return a gene_universe that matches the gene set in the",
      "gene_statistics used to construct it"
    )
  )
})

###############################################################################
