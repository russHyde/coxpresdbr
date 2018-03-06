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

test_that("evaluate_coex_partners-data.frame: valid input", {
  # test using both tibble and data.frame inputs

  # source gene has no partners: no-row data-frame returns
  pval <- data.frame(
    gene_id = letters[1:3],
    p_value = c(1e-12, 0.001, 0.9),
    direction = c(1L, -1L, 1L),
    stringsAsFactors = FALSE
  )
  coex_no_partners <- tibble::data_frame(
    source_id = "D",
    target_id = "E"
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
  expect_equal(
    object = evaluate_coex_partners(pval, coex_one_partner)$p_value,
    expected = 1e-12,
    info = paste(
      "source gene has a single partner => p-value should match that of the",
      "target gene"
    )
  )

  # p-value returned should be independent of `direction`
  # z-scores returned should be negatives if direction is switched

  # input p-values are identical for each target gene: returned p-value should
  # match the input p-vaue
})

###############################################################################
