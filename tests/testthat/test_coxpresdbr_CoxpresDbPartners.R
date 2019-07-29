###############################################################################

context("Tests for `CoxpresDbPartners` objects in `coxpresdbr` package")

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
      gene_statistics = tibble::tibble(
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
      partners = tibble::tibble(
        SOURCE_ID = "NOT A VALID `source_id` COLUMN",
        target_id = "valid_col_name"
      )
    ),
    info = "If not NULL/empty, partners should have a `source_id` column"
  )

  expect_error(
    object = new(
      "CoxpresDbPartners",
      partners = tibble::tibble(
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
      "The input to the `gene_statistics` field of `CoxpresDbPartners`",
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
      "A `CoxpresDbPartners` made from just a `gene_statistics`",
      "should return a `gene_universe` that matches the gene set in the",
      "`gene_statistics` used to construct it"
    )
  )

  test_partners_df <- tibble::tibble(
    source_id = rep(letters[1:3], each = 2),
    target_id = rep(LETTERS[6:4], times = 2)
  )

  expect_equal(
    object = get_source_genes(
      new("CoxpresDbPartners", partners = test_partners_df)
    ),
    expected = letters[1:3],
    info = paste(
      "A `CoxpresDbPartners` made from a partners data-frame should return",
      "the entries in the `source_id` column when `get_source_genes()` is",
      "called"
    )
  )
  expect_equal(
    object = get_source_genes(
      new("CoxpresDbPartners")
    ),
    expected = character(0),
    info = paste(
      "A `CoxpresDbPartners` with no `partners` field should return `NULL`",
      "when `get_source_genes()` is called"
    )
  )
})

###############################################################################

test_that("compute z scores from 2-tailed p and direction", {
  expect_error(
    object = .compute_z_scores(),
    info = "No input to `.compute_z_scores`"
  )

  p_vals <- c(0.5, 0.5, 0, 0, 1, 1, 1, 0.1, 0.1)
  dirs <- c(1, -1, 1, -1, 1, -1, 0, 1, -1)
  z_scores <- -1 * dirs * qnorm(p_vals / 2)

  expect_equal(
    object = .compute_z_scores(p_values = 1, directions = 0),
    expected = 0,
    info = "middling z-score"
  )
  expect_error(
    object = .compute_z_scores(p_values = p_vals),
    info = "No directions in `.compute_z_scores`"
  )
  expect_error(
    object = .compute_z_scores(directions = dirs),
    info = "No p-values in `.compute_z_scores`"
  )
  expect_error(
    object = .compute_z_scores(p_values = -1, directions = 1),
    info = "p-values should not be < 0"
  )
  expect_error(
    object = .compute_z_scores(p_values = 10, directions = 1),
    info = "p-values should not be > 1"
  )

  expect_equal(
    object = sign(
      .compute_z_scores(p_values = p_vals, directions = dirs)
    ),
    expected = sign(dirs) * (p_vals != 1),
    info = "signs of the z-scores should match the directions (for p < 1)"
  )
  expect_equal(
    object = -1 * .compute_z_scores(p_values = p_vals, directions = dirs),
    expected = .compute_z_scores(p_values = p_vals, directions = -1 * dirs),
    info = "negative-z-score is z-score of negative direction"
  )

  expect_equal(
    object = .compute_z_scores(p_vals, dirs),
    expected = z_scores,
    info = "hand-cranked z-score values"
  )
})

test_that(
  "add z-scores to the gene-statistics of an existing `CoxpresDbPartners`", {
    coex_partners <- new(
      "CoxpresDbPartners",
      gene_statistics = test_gene_statistics
    )

    expect_is(
      object = .add_z_scores(coex_partners),
      "CoxpresDbPartners",
      info = "`.add_z_scores` should return a `CoxpresDbPartners` object"
    )

    expect_equal(
      object = .add_z_scores(coex_partners)@gene_statistics$z_score,
      expected = .compute_z_scores(
        coex_partners@gene_statistics$p_value,
        coex_partners@gene_statistics$direction
      ),
      info = "add z-scores to a `CoxpresDbPartners`"
    )
  }
)

###############################################################################
