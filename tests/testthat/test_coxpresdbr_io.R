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

test_that(".is_coxpresdb_archive", {
  expect_equal(
    object = .is_coxpresdb_archive(test_data_file),
    expected = TRUE,
    info = "Checks a single valid coxpresdb archive is a valid archive"
  )

  expect_equal(
    object = .is_coxpresdb_archive("test_coxpresdbr_io.R"),
    expected = FALSE,
    info = "A valid coxpresdb archive should have a `.tar.bz2 extension`"
  )

  expect_equal(
    object = .is_coxpresdb_archive(rep(test_data_file, 2)),
    expected = FALSE,
    info = "User should only use a single coxpresdb archive at a time"
  )
})

###############################################################################

test_that(".get_coxpresdb_file_paths", {
  expect_equal(
    object = .get_coxpresdb_file_paths(test_data_file),
    expected = sort(file.path("spo_v14_subset", test_data_genes)),
    info = "File contents of a valid CoxpresDB.jp archive"
  )
})

###############################################################################

test_that("import_coex_db_universe: invalid input", {
  expect_error(
    object = import_coex_db_universe(db_archive = "NOT A FILE"),
    info = paste(
      "Attempt to load a gene-universe from a missing file as",
      "coxpresdb archive"
    )
  )

  expect_error(
    object = import_coex_db_universe(db_archive = "test_coxpresdbr_io.R"),
    info = paste(
      "Attempt to load a gene-universe from an existing file that",
      "is not a .tar.bz2 coxpresdb archive"
    )
  )
})

test_that("import_coex_db_universe: valid input", {
  expect_silent(
    object = import_coex_db_universe(db_archive = test_data_file)
  )

  expect_equal(
    object = import_coex_db_universe(db_archive = test_data_file),
    expected = test_data_genes,
    info = "Gene-set parsing from a coxpresdb archive"
  )
})

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

  expect_error(
    object = import_coex_db(test_data_genes[1:2], db_archive = test_data_file),
    info = paste(
      "User should only request the coexpression database for a",
      "single gene"
    )
  )
})

test_that("import_coex_db: valid input", {
  message(getwd())

  expect_silent(
    object = import_coex_db(gene_id = "2538791", db_archive = test_data_file)
  )

  coex_db_2538791 <- import_coex_db(
    gene_id = "2538791", db_archive = test_data_file
  )

  expect_is(
    object = coex_db_2538791,
    class = "tbl_df",
    info = "Coexpression database should be returned as a tibble::tbl_df"
  )

  expect_equal(
    object = colnames(coex_db_2538791),
    expected = c("source_id", "target_id", "mutual_rank", "correlation"),
    info = paste(
      "Colnames of a coexpression database should be `source_id`,",
      "`target_id`, `mutual_rank`, and `correlation`"
    )
  )

  expect_equal(
    object = unique(coex_db_2538791[["source_id"]]),
    expected = "2538791",
    info = paste(
      "A single source-gene should be present in a returned",
      "coexpression database"
    )
  )

  expect_false(
    object = "2538791" %in% coex_db_2538791[["target_id"]],
    info = paste(
      "The source-gene should be absent from the targets in its",
      "coexpression database"
    )
  )
  expect_equal(
    object = as.vector(sapply(coex_db_2538791, class)),
    expected = c("character", "character", "numeric", "numeric"),
    info = paste(
      "Coltypes should be character for source/target-id and numeric for",
      "mutual-rank/correlation"
    )
  )
})
