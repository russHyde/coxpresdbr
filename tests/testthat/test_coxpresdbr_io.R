###############################################################################

context("Tests for dataset I/O in `coxpresdbr` package")

###############################################################################

# - Test data contains files for 10 genes from the fission yeast coxpresdb
# dataset
# - Each file contains coexpression data for just those 10 genes
# - NOTE: testthat tests are ran with "<pkg>/tests/testthat/" as the working
# directory
test_data_bz2 <- "spo_v14_subset.tar.bz2"
test_data_tar <- "spo_v14_subset.tar"
test_data_zip <- "spo_v14_subset_two_cols.zip"

test_data_genes <- as.character(c(
  2538791, 2539499, 2540594, 2541560, 2541907,
  2542210, 2542294, 2543492, 3361219, 3361512
))

coex_data_2538791 <- tibble::tibble(
  # copied & reformatted from the file
  source_id = "2538791",
  target_id = c(
    "2543492", "2539499", "3361512", "2540594", "2541907", "2542294",
    "2542210", "2541560", "3361219"
  ),
  mutual_rank = c(
    631.52, 875.09, 1125.90, 1214.15, 1261.86, 2939.49,
    4110.55, 4184.63, 4250.63
  )
)

###############################################################################

test_df <- tibble::tibble(
  source_id = c("123", "123", "987", "123", "123"),
  target_id = c("123", "456", "123", "987", "568"),
  mutual_rank = c(1, 1.5, 2, 3.5, 2.2)
)

test_df_genes <- as.character(c(123, 456, 568, 987))

###############################################################################

test_that("CoxpresDbAccessor: Constructor", {

  # Ensure that overwrite_in_bunzip2 is true whenever running unit tests that
  # decompress an archive - otherwise, R.utils::bunzip2 will throw an exception
  # on re-extracting the same file.
  expect_is(
    object = CoxpresDbAccessor(
      db_archive = test_data_bz2,
      overwrite_in_bunzip2 = TRUE
    ),
    "CoxpresDbArchiveAccessor",
    info = paste(
      "Construction of a `CoxpresDbArchiveAccessor` from a bz2-compressed",
      "archive"
    )
  )

  expect_is(
    object = CoxpresDbAccessor(db_archive = test_data_tar),
    "CoxpresDbArchiveAccessor",
    info = paste(
      "Construction of a `CoxpresDbArchiveAccessor` from an uncompressed",
      "archive"
    )
  )

  expect_is(
    object = CoxpresDbAccessor(db_archive = test_data_zip),
    "CoxpresDbArchiveAccessor",
    info = paste(
      "Construction of a `CoxpresDbArchiveAccessor` from a zip-compressed",
      "archive"
    )
  )

  expect_error(
    object = CoxpresDbAccessor(db_archive = "NOT A FILE"),
    info = "A valid file must be passed as `db_archive` to `CoxpresDbAccessor`"
  )

  expect_error(
    object = CoxpresDbAccessor(db_archive = "test_coxpresdbr_io.R"),
    info = paste(
      "A valid (.tar.bz2, .tar or .zip) file must be passed as",
      "`db_archive` to `CoxpresDbAccessor`"
    )
  )
})

###############################################################################

test_that("CoxpresDbAccessor: accessors for filenames", {
  expect_equal(
    object = get_raw_archive(
      CoxpresDbAccessor(
        db_archive = test_data_bz2,
        overwrite_in_bunzip2 = TRUE
      )
    ),
    expected = test_data_bz2,
    info = paste(
      "raw_archive for a compressed archive should be the compressed",
      "archive's filename"
    )
  )

  expect_equal(
    object = get_raw_archive(
      CoxpresDbAccessor(db_archive = test_data_tar)
    ),
    expected = test_data_tar,
    info = paste(
      "raw_archive for a uncompressed .tar archive should be the uncompressed",
      "archive's filename"
    )
  )

  expect_equal(
    object = get_raw_archive(
      CoxpresDbAccessor(db_archive = test_data_zip)
    ),
    expected = test_data_zip,
    info = paste(
      "raw_archive for a .zip archive should be the .zip file itself"
    )
  )

  expect_equal(
    object = get_uncompressed_archive(
      CoxpresDbAccessor(
        db_archive = test_data_bz2,
        overwrite_in_bunzip2 = TRUE
      )
    ),
    expected = file.path(tempdir(), test_data_tar),
    info = paste(
      "uncompressed_archive for a compressed archive should be in the tempdir",
      "by default"
    )
  )

  expect_equal(
    object = get_uncompressed_archive(
      CoxpresDbAccessor(db_archive = test_data_tar)
    ),
    expected = test_data_tar,
    info = paste(
      "uncompressed_archive for a uncompressed archive should be the",
      "uncompressed archive's filename"
    )
  )

  expect_equal(
    object = get_uncompressed_archive(
      CoxpresDbAccessor(db_archive = test_data_zip)
    ),
    expected = test_data_zip,
    info = paste(
      "since we can access files from a .zip, we don't need to uncompress it",
      "so the uncompressed archive's filename should match the .zips filename"
    )
  )
})

###############################################################################

test_that(".is_coxpresdb_archive", {
  expect_true(
    object = .is_coxpresdb_archive(test_data_bz2),
    info = paste(
      "Checks a single, valid, coxpresdb archive (*.tar.bz2) is a valid",
      "archive"
    )
  )

  expect_true(
    object = .is_coxpresdb_archive(test_data_tar),
    info = paste(
      "Checks a single, valid, uncompressed coxpresdb (*.tar) is a valid",
      "archive"
    )
  )

  expect_true(
    .is_coxpresdb_archive(test_data_zip),
    info = paste(
      "A .zip archive is a valid coxpresdb archive"
    )
  )

  expect_false(
    object = .is_coxpresdb_archive("test_coxpresdbr_io.R"),
    info = paste(
      "A valid coxpresdb archive should have a `.tar`, `.tar.bz2` or `.zip`",
      "extension`"
    )
  )

  expect_false(
    object = .is_coxpresdb_archive(rep(test_data_bz2, 2)),
    info = "User should only use a single coxpresdb archive at a time"
  )
})

###############################################################################

test_that("get the file-paths for all genes in the archive", {
  expected_tar <- tibble::tibble(
    gene_id = test_data_genes,
    file_path = file.path("spo_v14_subset", test_data_genes)
  ) %>%
    dplyr::arrange(gene_id)

  expected_zip <- tibble::tibble(
    gene_id = test_data_genes,
    file_path = file.path("spo_v14_subset_two_cols", test_data_genes)
  ) %>%
    dplyr::arrange(gene_id)

  bz2_accessor <- CoxpresDbAccessor(test_data_bz2, overwrite_in_bunzip2 = TRUE)
  tar_accessor <- CoxpresDbAccessor(test_data_tar)
  zip_accessor <- CoxpresDbAccessor(test_data_zip)

  expect_equal(
    object = get_file_paths(bz2_accessor),
    expected = expected_tar,
    info = paste(
      "File paths for the gene-partner datafames in a compressed",
      "CoxpresDB archive"
    )
  )

  expect_equal(
    object = get_file_paths(tar_accessor),
    expected = expected_tar,
    info = paste(
      "File paths for the gene-partner datafames in an",
      "uncompressed CoxpresDB archive"
    )
  )

  expect_equal(
    object = get_file_paths(zip_accessor),
    expected = expected_zip,
    info = paste(
      "File paths for the gene-partner dataframes as obtained from a .zip",
      "CoxpresDB archive"
    )
  )
})

###############################################################################

test_that("get file path for a specific gene", {
  bz2_accessor <- CoxpresDbAccessor(test_data_bz2, overwrite_in_bunzip2 = TRUE)
  zip_accessor <- CoxpresDbAccessor(test_data_zip)

  expect_equal(
    object = get_file_path_for_gene("2538791", bz2_accessor),
    expected = file.path("spo_v14_subset", "2538791"),
    info = "file path for a specific gene in a CoxpresDb archive"
  )

  expect_equal(
    object = get_file_path_for_gene("2538791", zip_accessor),
    expected = file.path("spo_v14_subset_two_cols", "2538791"),
    info = "file path for a specific gene in a .zip Coxpresdb archive"
  )
})

###############################################################################

test_that(
  "get the genes that are defined in a CoxpresDB archive: valid input", {
    bz2_accessor <- CoxpresDbAccessor(
      db_archive = test_data_bz2, overwrite_in_bunzip2 = TRUE
    )
    tar_accessor <- CoxpresDbAccessor(
      test_data_tar, overwrite_in_bunzip2 = TRUE
    )
    zip_accessor <- CoxpresDbAccessor(test_data_zip)
    df_accessor <- CoxpresDbAccessor(test_df)

    expect_silent(
      object = get_gene_ids(bz2_accessor)
    )

    expect_equal(
      object = get_gene_ids(bz2_accessor),
      expected = test_data_genes,
      info = "Accessor test for gene-ids from a compressed CoxpresDB archive"
    )

    expect_equal(
      object = get_gene_ids(tar_accessor),
      expected = test_data_genes,
      info = "Accessor test for gene-ids from an uncompressed CoxpresDB archive"
    )

    expect_equal(
      object = get_gene_ids(zip_accessor),
      expected = test_data_genes,
      info = "Accessor test for gene-ids from a .zip CoxpresDb archive"
    )

    expect_equal(
      object = get_gene_ids(df_accessor),
      expected = test_df_genes,
      info = "Accessor test for gene-ids from a data-frame based CoxpresDb"
    )
  }
)

###############################################################################

test_that("get_all_coex_partners: invalid input", {
  bz2_accessor <- CoxpresDbAccessor(
    db_archive = test_data_bz2, overwrite_in_bunzip2 = TRUE
  )

  expect_error(
    object = get_all_coex_partners(
      gene_id = "2538791", importer = "NOT AN IMPORTER"
    ),
    info = paste(
      "Attempt to load gene-partners from a string, not an `CoxpresDbAccessor`"
    )
  )

  expect_error(
    object = get_all_coex_partners("NOT_A_GENE", importer = bz2_accessor),
    info = "Attempt to import a missing gene from a coxpresdb archive"
  )

  expect_error(
    object = get_all_coex_partners(
      test_data_genes[1:2], importer = bz2_accessor
    ),
    info = paste(
      "User should only request the coexpression database for a",
      "single gene"
    )
  )
})

###############################################################################

test_that("get_all_coex_partners: input from an archive", {
  bz2_accessor <- CoxpresDbAccessor(test_data_bz2, overwrite_in_bunzip2 = TRUE)
  zip_accessor <- CoxpresDbAccessor(test_data_zip)

  expect_equal(
    object = get_all_coex_partners(
      gene_id = "2538791", importer = bz2_accessor
    ),
    expected = coex_data_2538791,
    info = paste(
      "Given: a single gene ID and a bz2-based CoexpresDB accessor;",
      "When: the user requests all coexpression partners of that gene;",
      "Then: the returned data should have `source_id`, `target_id` and",
      "`mutual_rank` columns, the source gene should be absent from the",
      "target gene columns and the mutual-rank values should be presented in",
      "increasing order"
    )
  )

  expect_equal(
    object = get_all_coex_partners(
      gene_id = "2538791", importer = zip_accessor
    ),
    expected = coex_data_2538791,
    info = paste(
      "Given: a single gene ID and a zip-based CoexpresDB accessor;",
      "When: the user requests all coexpression partners of that gene;",
      "Then: the returned data should have `source_id`, `target_id` and",
      "`mutual_rank` columns, the source gene should be absent from the target",
      "gene columns and the mutual rank values should be presented in
      increasing order"
    )
  )
})

###############################################################################

# Rather than unit testing `get_all_coex_partners` for a data-frame-based
# CoxpresDB archive, recommend assessing the exported function
# `get_coex_partners`

###############################################################################
