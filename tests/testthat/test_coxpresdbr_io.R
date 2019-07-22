###############################################################################

context("Tests for dataset I/O in `coxpresdbr` package")

###############################################################################

# - Test data contains files for 10 genes from the fission yeast coxpresdb
# dataset
# - Each file contains coexpression data for just those 10 genes
# - NOTE: testthat tests are ran with "<pkg>/tests/testthat/" as the working
# directory
test_data_file <- "spo_v14_subset.tar.bz2"
test_data_uncompressed <- "spo_v14_subset.tar"
test_data_genes <- as.character(c(
  2538791, 2539499, 2540594, 2541560, 2541907,
  2542210, 2542294, 2543492, 3361219, 3361512
))

###############################################################################

test_that("CoxpresDbImporter: Constructor", {

  # Ensure that overwrite_in_bunzip2 is true whenever running unit tests that
  # decompress an archive - otherwise, R.utils::bunzip2 will throw an exception
  # on re-extracting the same file.
  expect_is(
    object = CoxpresDbImporter(
      db_archive = test_data_file,
      overwrite_in_bunzip2 = TRUE
    ),
    "CoxpresDbImporter",
    info = "Construction of a CoxpresDbImporter with a compressed archive"
  )

  expect_is(
    object = CoxpresDbImporter(db_archive = test_data_uncompressed),
    "CoxpresDbImporter",
    info = "Construction of a CoxpresDbImporter with an uncompressed archive"
  )

  expect_error(
    object = CoxpresDbImporter(db_archive = "NOT A FILE"),
    info = "A valid file must be passed as db_archive to CoxpresDbImporter"
  )

  expect_error(
    object = CoxpresDbImporter(db_archive = "test_coxpresdbr_io.R"),
    info = "A valid file must be passed as db_archive to CoxpresDbImporter"
  )
})

###############################################################################

test_that("CoxpresDbImporter: accessors for filenames", {
  expect_equal(
    object = get_raw_archive(
      CoxpresDbImporter(
        db_archive = test_data_file,
        overwrite_in_bunzip2 = TRUE
      )
    ),
    expected = test_data_file,
    info = paste(
      "raw_archive for a compressed archive should be the compressed",
      "archive's filename"
    )
  )

  expect_equal(
    object = get_raw_archive(
      CoxpresDbImporter(db_archive = test_data_uncompressed)
    ),
    expected = test_data_uncompressed,
    info = paste(
      "raw_archive for a uncompressed archive should be the uncompressed",
      "archive's filename"
    )
  )

  expect_equal(
    object = get_uncompressed_archive(
      CoxpresDbImporter(
        db_archive = test_data_file,
        overwrite_in_bunzip2 = TRUE
      )
    ),
    expected = file.path(tempdir(), test_data_uncompressed),
    info = paste(
      "uncompressed_archive for a compressed archive should be in the tempdir"
    )
  )

  expect_equal(
    object = get_uncompressed_archive(
      CoxpresDbImporter(db_archive = test_data_uncompressed)
    ),
    expected = test_data_uncompressed,
    info = paste(
      "uncompressed_archive for a uncompressed archive should be the",
      "uncompressed archive's filename"
    )
  )
})

###############################################################################

test_that(".is_coxpresdb_archive", {
  expect_true(
    object = .is_coxpresdb_archive(test_data_file),
    info = paste(
      "Checks a single, valid, coxpresdb archive (*.tar.bz2) is a valid",
      "archive"
    )
  )

  expect_true(
    object = .is_coxpresdb_archive(test_data_uncompressed),
    info = paste(
      "Checks a single, valid, uncompressed coxpresdb (*.tar) is a valid",
      "archive"
    )
  )

  expect_false(
    object = .is_coxpresdb_archive("test_coxpresdbr_io.R"),
    info = "A valid coxpresdb archive should have a `.tar.bz2 extension`"
  )

  expect_false(
    object = .is_coxpresdb_archive(rep(test_data_file, 2)),
    info = "User should only use a single coxpresdb archive at a time"
  )
})

###############################################################################

test_that("get the file-paths for all genes in the archive", {
  expected <- tibble::data_frame(
    gene_id = test_data_genes,
    file_path = file.path("spo_v14_subset", test_data_genes)
  ) %>%
    dplyr::arrange(gene_id)

  expect_equal(
    object = get_file_paths(
      CoxpresDbImporter(test_data_file, overwrite_in_bunzip2 = TRUE)
    ),
    expected = expected,
    info = paste(
      "File paths for the gene-partner datafames in a compressed",
      "CoxpresDB archive"
    )
  )

  expect_equal(
    object = get_file_paths(
      CoxpresDbImporter(test_data_uncompressed)
    ),
    expected = expected,
    info = paste(
      "File paths for the gene-partner datafames in an",
      "uncompressed CoxpresDB archive"
    )
  )
})

test_that("get file path for a specific gene", {
  expect_equal(
    object = get_file_path_for_gene(
      "2538791",
      CoxpresDbImporter(test_data_file, overwrite_in_bunzip2 = TRUE)
    ),
    expected = file.path("spo_v14_subset", "2538791"),
    info = "file path for a specific gene in a CoxpresDb archive"
  )
})

###############################################################################

test_that(
  "get the genes that are defined in a CoxpresDb archive: invalid input", {
    expect_error(
      object = get_gene_ids(
        CoxpresDbImporter(db_archive = "NOT A FILE")
      ),
      info = paste(
        "Attempt to load a gene-universe from a missing file as",
        "coxpresdb archive"
      )
    )

    expect_error(
      object = get_gene_ids(
        CoxpresDbImporter(db_archive = "test_coxpresdbr_io.R")
      ),
      info = paste(
        "Attempt to load a gene-universe from an existing file that",
        "is not a .tar.bz2 coxpresdb archive"
      )
    )
  }
)

###############################################################################

test_that(
  "get the genes that are defined in a CoxpresDB archive: valid input", {
    expect_silent(
      object = get_gene_ids(
        CoxpresDbImporter(
          db_archive = test_data_file, overwrite_in_bunzip2 = TRUE
        )
      )
    )

    expect_equal(
      object = get_gene_ids(
        CoxpresDbImporter(test_data_file, overwrite_in_bunzip2 = TRUE)
      ),
      expected = test_data_genes,
      info = "Accessor test for gene-ids from a compressed CoxpresDB archive"
    )

    expect_equal(
      object = get_gene_ids(
        CoxpresDbImporter(test_data_uncompressed, overwrite_in_bunzip2 = TRUE)
      ),
      expected = test_data_genes,
      info = "Accessor test for gene-ids from an uncompressed CoxpresDB archive"
    )
  }
)

###############################################################################

test_that("import_all_coex_partners: invalid input", {
  importer <- CoxpresDbImporter(
    db_archive = test_data_file, overwrite_in_bunzip2 = TRUE
  )

  expect_error(
    object = import_all_coex_partners(
      gene_id = "2538791", importer = "NOT AN IMPORTER"
    ),
    info = "Attempt to load gene-partners from a string, not an importer"
  )

  expect_error(
    object = import_all_coex_partners("NOT_A_GENE", importer = importer),
    info = "Attempt to import a missing gene from a coxpresdb archive"
  )

  expect_error(
    object = import_all_coex_partners(
      test_data_genes[1:2],
      importer = importer
    ),
    info = paste(
      "User should only request the coexpression database for a",
      "single gene"
    )
  )
})

test_that("import_all_coex_partners: valid input", {
  message(getwd())

  importer <- CoxpresDbImporter(test_data_file, overwrite_in_bunzip2 = TRUE)

  # expect_silent(
  #  object = import_all_coex_partners(
  #    gene_id = "2538791", importer = importer
  #  )
  # )

  coex_db_2538791 <- import_all_coex_partners(
    gene_id = "2538791", importer = importer
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
    object = vapply(
      coex_db_2538791, function(x) class(x)[1],
      FUN.VALUE = character(1),
      USE.NAMES = FALSE
    ),
    expected = c("character", "character", "numeric", "numeric"),
    info = paste(
      "Coltypes should be character for source/target-id and numeric for",
      "mutual-rank/correlation"
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
})
