###############################################################################

context("Tests for compute_cor")

###############################################################################

random_matrix <- function(
  n_rows = sample(1:50, 1),
  n_cols = sample(1:50, 1),
  row_names = paste0("R", seq_len(n_rows)),
  col_names = paste0("C", seq_len(n_cols))
) {
  matrix(
    rnorm(n = n_rows * n_cols),
    nrow = n_rows, ncol = n_cols, dimnames = list(row_names, col_names)
  )
}

###############################################################################

.z_transform <- function(r) {
  # different but equivalent to that used in the source code
  0.5 * (log(1 + r) - log(1 - r))
}

format_results <- function(
  source_id = character(0),
  target_id = character(0),
  partial_cor = numeric(max(length(source_id), length(target_id))),
  z_score = .z_transform(partial_cor)
) {
  tibble::tibble(
    source_id = source_id,
    target_id = target_id,
    partial_cor = partial_cor,
    z_score = z_score
  )
}

empty_result <- format_results()

# Test for
# - identity of the character columns, and
# - numeric identity of the numeric columns (allowing numeric-tolerance)
# .. returned by compute_cor()

expect_cor <- function(object, expected, ..., info = NULL) {
  expect_equal(
    object[c("source_id", "target_id")],
    expected[c("source_id", "target_id")],
    info = info
  )
  expect_equal(
    as.matrix(object[, c("partial_cor", "z_score")]),
    as.matrix(expected[, c("partial_cor", "z_score")]),
    ...,
    info = info
  )
}

###############################################################################

test_that(
  "compute_cor returns the (partial) correlation for a set of column pairs", {

  mat <- random_matrix(n_cols = 2, col_names = c("a", "b"))

  # Here we use a one-column binary control:
  # [1, 0, 1, 0, 1, 0, 1]^T
  #
  # The user could use
  # [1, 0, 1, 0, 1, 0, 1;
  #  0, 1, 0, 1, 0, 1, 0]^T
  # .. and would get the same partial-correlation values, but ppcor will
  # complain about the variance-covariance matrix being singular
  control <- matrix(
    rep(c(1, 0), length.out = nrow(mat)),
    nrow = nrow(mat),
    ncol = 1
  )

  # no source/target pairs
  expect_cor(
    compute_cor(
      x = mat, source_id = character(0), target_id = character(0),
      control = NULL
    ),
    empty_result
  )

  # a single source/target pair with no controlling variables
  s0 <- "a"
  t0 <- "b"
  cor0 <- cor(mat[, s0], mat[, t0])
  expect_cor(
    compute_cor(x = mat, source_id = s0, target_id = t0, control = NULL),
    format_results(s0, t0, cor0),
    info = paste(
      "Without a control matrix, the reported `partial_cor` is",
      "cor(source_column, target_column)"
    )
  )

  # a single source/target pair with a control matrix
  expect_cor(
    compute_cor(
      x = mat, source_id = s0, target_id = t0, control = control
    ),
    format_results(
      s0, t0, ppcor::pcor.test(mat[, s0], mat[, t0], control)$estimate
    ),
    info = paste(
      "With a control matrix, the reported `partial_cor` is computed using",
      "ppcor"
    )
  )

  expect_error(
    compute_cor(x = mat, source_id = "a", target_id = "a"),
    info = "source/target pairs should be distinct"
  )
})

test_that("compute_cor for multiple source/target pairs", {
  #
  mat <- random_matrix(n_col = 3)

  # note that we use cor(vec, matrix) to get the correlation between vec and
  # each column of matrix as a shorthand here

  s0 <- "C1"
  t0 <- c("C2", "C3")
  expected <- format_results(
    s0, t0, as.vector(cor(mat[, s0], mat[, t0]))
  )
  expect_cor(
    compute_cor(x = mat, source_id = s0, target_id = t0),
    expected,
    tolerance = 1, scale = 1,
    info = "compute correlation between one source-id and several target-ids"
  )

  s0 <- c("C1", "C3")
  t0 <- "C2"
  expected <- format_results(
    s0, t0, as.vector(cor(mat[, t0], mat[, s0]))
  )
  expect_cor(
    compute_cor(x = mat, source_id = s0, target_id = t0),
    expected,
    tolerance = 1, scale = 1,
    info = "compute correlation between several source-ids and one target-id"
  )

  expect_error(
    compute_cor(x = mat, source_id = "C1", target_id = c("C2", "C3", "C1")),
    info = "can't correlate source against itself (1 source, several targets)"
  )
  expect_error(
    compute_cor(x = mat, source_id = c("C3", "C2", "C1"), target_id = "C1"),
    info = "can't correlate source against itself (1 target, several sources)"
  )
})

###############################################################################
