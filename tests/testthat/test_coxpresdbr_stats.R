###############################################################################
# Tests for statistics functions in `coxpresdbr` package
#
###############################################################################

context("Tests for geneset-combination statistics in `coxpresdbr` package")

###############################################################################

test_that("evaluate_coex_partners-data.frame: invalid input", {
  pval_empty <- tibble::data_frame(
    gene_id = character(0),
    p_value = numeric(0),
    direction = integer(0)
  )

  coex_empty <- tibble::data_frame(
    source_id = character(0),
    target_id = character(0)
  )

  # `x` should inherit from data.frame
  # ... (or on generalising be a DGELRT/MArrayLM object)
  expect_error(
    object = evaluate_coex_partners(
      x = "NOT A DATA FRAME",
      coex_partners = coex_empty
    ),
    info = "`x` should be a data-frame in evaluate_coex_partners"
  )

  expect_error(
    object = evaluate_coex_partners(
      x = tibble::data_frame(
        GENE_ID = character(0),
        P_VALUE = numeric(0),
        DIRECTION = integer(0)
      )
    ),
    info = "`x` should have columns `gene_id`, `p_value` and `direction`"
  )

  # `coex_partners` should inherit from data.frame
  expect_error(
    object = evaluate_coex_partners(
      x = pval_empty,
      coex_partners = "Not a data-frame"
    ),
    info = "`coex_partners` should be a data-frame"
  )

  expect_error(
    object = evaluate_coex_partners(
      x = pval_empty,
      coex_partners = tibble::data_frame(
        SOURCE_ID = character(0),
        TARGET_ID = character(0)
      )
    ),
    info = "`coex_partners` should have columns `source_id` and `target_id`"
  )
})

###############################################################################

test_that("evaluate_coex_partners-data.frame: valid input", {
  # source gene has no partners: no-row data-frame returns
  pval <- data.frame(
    gene_id = letters[1:3],
    p_value = c(1e-12, 0.001, 0.9),
    direction = c(1L, -1L, 1L),
    stringsAsFactors = FALSE
  )
  coex_no_partners <- tibble::data_frame(
    source_id = "d",
    target_id = "e"
  )
  expect <- tibble::data_frame(
    gene_id = character(0),
    n_partners = integer(0),
    z_score = numeric(0),
    p_value = numeric(0)
  )
  expect_equal(
    object = evaluate_coex_partners(pval, coex_no_partners),
    expected = expect,
    info = "source gene has no partners => source gene is absent from output"
  )

  # source gene has a single partner: the returned p-value should match the
  # input p-value
  coex_one_partner <- tibble::data_frame(
    source_id = "d",
    target_id = "a"
  )
  expect <- tibble::data_frame(
    gene_id = "d",
    n_partners = 1L,
    z_score = 0.5 * (
      qnorm(5e-13, lower.tail = FALSE) -
        qnorm(1 - 5e-13, lower.tail = FALSE)
    ),
    p_value = 1e-12
  )
  expect_equal(
    object = as.data.frame(
      evaluate_coex_partners(pval, coex_one_partner)
    )[, 1:2],
    expected = as.data.frame(expect)[, 1:2],
    info = paste(
      "source gene has a single partner => gene-id/#partners for the result",
      "should match that of the target gene (non-numeric cols)"
    )
  )
  expect_equal(
    object = evaluate_coex_partners(pval, coex_one_partner)$z_score,
    expected = expect$z_score,
    info = paste(
      "source gene has a single partner => z-score should match that of the",
      "target gene (z-score col)"
    ),
    tolerance = 1e-4,
    scale = expect$z_score
  )
  expect_equal(
    object = evaluate_coex_partners(pval, coex_one_partner)$p_value,
    expected = expect$p_value,
    info = paste(
      "source gene has a single partner => p-value should match that of the",
      "target gene (p_value col)"
    ),
    tolerance = 1e-4,
    scale = expect$p_value
  )

  # p-value returned should be independent of `direction` and
  # z-scores returned should be reversed if `direction` is switched
  expect_equal(
    object = pval %>%
      dplyr::mutate(direction = -1 * direction) %>%
      evaluate_coex_partners(coex_one_partner) %>%
      as.data.frame(),

    expected = pval %>%
      evaluate_coex_partners(coex_one_partner) %>%
      dplyr::mutate(z_score = -1 * z_score) %>%
      as.data.frame(),

    info = paste(
      "if all directions change, the p-value should remain the same and the",
      "z-score should reverse"
    )
  )

  # TODO: input p-values are identical for each target gene: returned p-value
  # should match the input p-vaue / sqrt(num_targets)
})

###############################################################################

test_that("cluster_by_coex_partnership: invalid input", {
  # - Input should be a CoxpresDbPartners object
  # - The CoxpresDbPartners object should have a valid `gene_statistics`
  # - The CoxpresDbPartners object should have a valid `partners`
  # - `drop_disparities` should be Boolean

  expect_error(
    object = cluster_by_coex_partnership(
      coex_partners = "NOT A CoxpresDbPartners object",
      drop_disparities = TRUE
    ),
    info = "`coex_partners` should be a `CoxpresDbPartners` object"
  )

  expect_error(
    object = cluster_by_coex_partnership(
      coex_partners = new(
        "CoxpresDbPartners",
        partners = tibble::data_frame(source_id = "a", target_id = "b")
      ),
      drop_disparities = TRUE
    ),
    info = "`coex_partners` should have a valid `gene_statistics` entry"
  )

  expect_error(
    object = cluster_by_coex_partnership(
      coex_partners = new(
        "CoxpresDbPartners",
        gene_statistics = tibble::data_frame(
          gene_id = "abc", p_value = "0.2", direction = "1"
        )
      ),
      drop_disparities = TRUE
    ),
    info = "`coex_partners` should have a valid/non-empty `partners` entry"
  )

  test_coex_partners <- new(
    "CoxpresDbPartners",
    gene_statistics = tibble::data_frame(
      gene_id = letters[1:2], p_value = c(0.1, 0.2), direction = c(-1, 1)
    ),
    partners = tibble::data_frame(
      source_id = "a",
      target_id = "b"
    )
  )

  expect_error(
    object = cluster_by_coex_partnership(
      coex_partners = test_coex_partners,
      drop_disparities = "NOT A BOOLEAN"
    ),
    info = "`drop_disparities` should be Boolean"
  )
})

###############################################################################

test_that("cluster_by_coex_partnership: valid input", {

})

###############################################################################
