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

  expect_error(
    object = new(
      "CoxpresDbPartners",
      gene_statistics = tibble::data_frame(
        GENE_ID = "NOT gene_id", P_VALUE = -1, DIRECTION = Inf
      )
    ),
    info = paste(
      "If not NULL/empty, `gene_statistics` should have `gene_id`, `p_value`",
      "and `direction` columns"
    )
  )

  expect_error(
    object = new(
      "CoxpresDbPartners",
      partners = tibble::data_frame(
        SOURCE_ID = "NOT A VALID `source_id` COLUMN",
        target_id = "valid_col_name"
      )
    ),
    info = "If not NULL/empty, partners should have a `source_id` column"
  )

  expect_error(
    object = new(
      "CoxpresDbPartners",
      partners = tibble::data_frame(
        source_id = "valid_col_name",
        TaRget_Id = "NOT A VALID target_id COLUMN"
      )
    ),
    info = "If not NULL/empty, partners should have a `target_id` column"
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
})

test_that("CoxpresDbPartners: method checks", {
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

  test_partners_df <- tibble::data_frame(
    source_id = rep(letters[1:3], each = 2),
    target_id = rep(LETTERS[6:4], times = 2)
  )

  expect_equal(
    object = get_source_genes(
      new("CoxpresDbPartners", partners = test_partners_df)
    ),
    expected = letters[1:3],
    info = paste(
      "A CoxpresDbPartners made from a partners data-frame should return",
      "the entries in the source_id column when get_source_genes() is called"
    )
  )
  expect_equal(
    object = get_source_genes(
      new("CoxpresDbPartners")
    ),
    expected = character(0),
    info = paste(
      "A `CoxpresDbPartners` with no `partners` field should return NULL",
      "when get_source_genes() is called"
    )
  )
})

###############################################################################
