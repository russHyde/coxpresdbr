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
  # `coex_partners` should have columns `source_id` and `target_id`
  expect_error(
    object = evaluate_coex_partners(
      x = pval_empty,
      coex_partners = tibble::data_frame(
        SOURCE_ID = character(0),
        TARGET_ID = character(0)
        ),
      info = "`coex_partners` should have columns `source_id` and `target_id`"
      )
    )
})

test_that("evaluate_coex_partners-data.frame: valid input", {
  # test using both tibble and data.frame inputs
})

###############################################################################
