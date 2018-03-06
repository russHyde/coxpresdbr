###############################################################################

if (requireNamespace("lintr", quietly = TRUE)) {
  context("Tests for lints in `coxpresdbr` package")
  test_that("Package Style", {
    lintr::expect_lint_free()
  })
}

###############################################################################
