###############################################################################

context("Tests for dataset I/O in `coxpresdbr` package")

###############################################################################
# - Test data contains files for 10 genes from the fission yeast coxpresdb
# dataset
# - Each file contains coexpression data for just those 10 genes
# - NOTE: testthat tests are ran with "<pkg>/tests/testthat/" as the working
# directory
test_data_file <- "spo_v14_subset.tar.bz2"
test_data_genes <- as.character(c(
  2538791, 2539499, 2540594, 2541560, 2541907,
  2542210, 2542294, 2543492, 3361219, 3361512
))

###############################################################################

test_that("import_coex_db: invalid input", {
  expect_error(
    object = import_coex_db(gene_id = "2538791", db_archive = "NOT A FILE"),
    info = "Attempt to load a missing file as coxpresdb archive"
  )

  expect_error(
    object = import_coex_db("NOT_A_GENE", db_archive = test_data_file),
    info = "Attempt to import a missing gene from a coxpresdb archive"
  )
})
