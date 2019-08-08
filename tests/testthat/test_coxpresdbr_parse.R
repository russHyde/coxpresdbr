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
  # TODO: convert tests for `get_all_coex_partners` into tests on
  #   `get_coex_partners` since `get_all_coex_partners` is not exported

  test_data_file <- "spo_v14_subset.tar.bz2"
  test_data_genes <- as.character(c(
    2538791, 2539499, 2540594, 2541560, 2541907,
    2542210, 2542294, 2543492, 3361219, 3361512
  ))
  importer <- CoxpresDbAccessor(test_data_file, overwrite_in_bunzip2 = TRUE)

  expect_equal(
    object = get_coex_partners(
      gene_ids = "2538791", importer = importer
    ),
    expected = get_all_coex_partners(
      gene_id = "2538791", importer = importer
    ),
    info = paste(
      "for a single gene and no filters, get_coex_partners should match",
      "get_all_coex_partners"
    )
  )

  expect_equal(
    object = get_coex_partners(
      gene_ids = "2538791", importer = importer, n_partners = 3
    ),
    expected = get_all_coex_partners(
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
      get_all_coex_partners("2538791", importer = importer),
      get_all_coex_partners("2539499", importer = importer)
    ),
    info = paste(
      "a pair of source genes, without any filters"
    )
  )
})

###############################################################################

test_that("get_coex_partners from a dataframe", {
  # --- test data
  test_df <- tibble::tibble(
    source_id = c("123", "123", "987", "123", "123", "987"),
    target_id = c("123", "456", "123", "987", "568", "456"),
    mutual_rank = c(1, 1.5, 2, 3.5, 2.2, 1.2)
  )

  test_df_123_partners <- tibble::tibble(
    source_id = "123",
    target_id = c("456", "568", "987"),
    # note this is sorted by mutual rank, even though the input (eg, for
    # source = 123) isn't sorted;
    # also note source = 987 is absent from the source id column
    # also note target = 123 is absent from the target ids
    mutual_rank = c(1.5, 2.2, 3.5)
  )

  test_df_987_partners <- tibble::tibble(
    source_id = "987",
    target_id = c("456", "123"),
    mutual_rank = c(1.2, 2)
  )

  df_accessor <- CoxpresDbAccessor(test_df)

  # --- tests

  expect_equal(
    object = get_coex_partners(
      gene_ids = "123", importer = df_accessor
    ),
    expected = test_df_123_partners,
    info = paste(
      "Given: a dataframe-based CoxpresDbAccessor and a single gene ID;",
      "When: the user requests the coexpression partners of that gene;",
      "Then: the coexpression partners are returned in order of increasing",
      "mutual rank and the source gene is not a partner of itself"
    )
  )

  expect_error(
    object = get_coex_partners(
      gene_ids = "NOT PRESENT", importer = df_accessor
    ),
    info = paste(
      "Given: a dataframe-based CoxpresDbAccessor and a gene that is not",
      "present in that dataframe;",
      "When: the user requests the coexpression partners of that gene;",
      "Then: an error is thrown"
    )
  )

  expect_error(
    object = get_coex_partners(
      gene_ids = "456", importer = df_accessor
    ),
    info = paste(
      "Given: a dataframe-based CoxpresDbAccessor and a gene that is not",
      "a 'source' gene (but is a 'target' gene);",
      "When: the user requests the coexpression partners of that gene;",
      "Then: an error is thrown"
    )
  )

  expect_equal(
    object = get_coex_partners(
      gene_ids = c("987", "123"), importer = df_accessor
    ),
    expected = dplyr::bind_rows(
      test_df_987_partners, test_df_123_partners
    ),
    info = paste(
      "Given: a dataframe-based CoxpresDbAccessor and two genes that are",
      "'source' genes in that dataframe;",

      "When: the user requests the coexpression partners of those genes;",

      "Then: the coexpression partners are returned in source-gene associated",
      "blocks, and the source genes are ordered as in the input"
    )
  )

  expect_equal(
    object = get_coex_partners(
      gene_ids = c("987", "123"), importer = df_accessor, n_partners = 1
    ),
    expected = dplyr::bind_rows(
      test_df_987_partners[1, ], test_df_123_partners[1, ]
    ),
    info = paste(
      "Given: a dataframe-based CoxpresDbAccessor, and two genes that are",
      "'source' genes in that dataframe;",

      "When: the user requests the top coexpression partner of each gene;",

      "Then: a single coexpression partner is returned in a dataframe ordered",
      "by the input order of the source-genes"
    )
  )

  expect_equal(
    object = get_coex_partners(
      gene_ids = c("987", "123"), importer = df_accessor, mr_threshold = 2.1
    ),
    expected = dplyr::bind_rows(
      dplyr::filter(test_df_987_partners, mutual_rank <= 2.1),
      dplyr::filter(test_df_123_partners, mutual_rank <= 2.1)
    ),
    info = paste(
      "Given: a dataframe-based CoxpresDbAccessor, and two genes that are",
      "'source' genes in that dataframe;",

      "When: the user requests coexpression partners up to a mutual-rank",
      "threshold;",

      "Then: all coexpression partners returned should have mutual-rank <=",
      "that threshold"
    )
  )
})
