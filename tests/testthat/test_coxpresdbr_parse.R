###############################################################################
# Tests for parsing functions for `coxpresdbr` package
#
###############################################################################

context("Tests for dataset parsing in `coxpresdbr` package")

###############################################################################

test_that(".filter_coex_partners", {
  df1 <- tibble::tibble(
    source_id = rep("A", 5),
    target_id = LETTERS[2:6],
    mutual_rank = 1:5
  )

  expect_equal(
    object = .filter_coex_partners(df1),
    expected = df1,
    info = "dataset should be unchanged if no filters are applied"
  )

  expect_equal(
    object = .filter_coex_partners(df1, gene_universe = c("A", "B", "D")),
    expected = df1[c(1, 3), ],
    info = "gene_universe imposes a filter over both source_id and target_id"
  )

  expect_equal(
    object = .filter_coex_partners(df1, gene_universe = c("B", "D")),
    expected = df1[integer(0), ],
    info = paste(
      "if none of the source_ids are in the gene_universe, an empty frame",
      "results"
    )
  )

  expect_equal(
    object = .filter_coex_partners(
      df1[sample(seq(nrow(df1))), ],
      n_partners = 2
    ),
    expected = df1[1:2, ],
    info = "`n_partners` imposes a filter on the number of returned entries"
  )

  expect_equal(
    object = .filter_coex_partners(
      df1[sample(seq(nrow(df1))), ],
      mr_threshold = 4
    ),
    expected = df1[1:4, ],
    info = "`mr_threshold` imposes a filter on the mutual_rank value"
  )
})

###############################################################################

test_that("get_coex_partners", {
  # TODO: tests for multiple genes and filtering

  test_data_file <- "spo_v14_subset.tar.bz2"
  test_data_genes <- as.character(c(
    2538791, 2539499, 2540594, 2541560, 2541907,
    2542210, 2542294, 2543492, 3361219, 3361512
  ))
  importer <- CoxpresDbImporter(test_data_file, overwrite_in_bunzip2 = TRUE)

  expect_equal(
    object = get_coex_partners(
      gene_ids = "2538791", importer = importer
    ),
    expected = import_all_coex_partners(
      gene_id = "2538791", importer = importer
    ),
    info = paste(
      "for a single gene and no filters, get_coex_partners should match",
      "import_all_coex_partners"
    )
  )

  expect_equal(
    object = get_coex_partners(
      gene_ids = "2538791", importer = importer, n_partners = 3
    ),
    expected = import_all_coex_partners(
      gene_id = "2538791", importer = importer
    )[1:3, ],
    info = paste(
      "a single gene with a `n_partners` filter"
    )
  )

  expect_equal(
    object = get_coex_partners(
      gene_ids = c("2538791", "2539499"), importer = importer
    ),
    expected = bind_rows(
      import_all_coex_partners("2538791", importer = importer),
      import_all_coex_partners("2539499", importer = importer)
    ),
    info = paste(
      "a pair of source genes, without any filters"
    )
  )
})
