###############################################################################

context("Tests for workflows in `coxpresdbr` package")

###############################################################################

set.seed(1)

###############################################################################

test_data_file <- "spo_v14_subset.tar.bz2"

test_data_uncompressed <- "spo_v14_subset.tar"

test_data_genes <- as.character(c(
  2538791, 2539499, 2540594, 2541560, 2541907,
  2542210, 2542294, 2543492, 3361219, 3361512
))

test_gene_statistics <- tibble::data_frame(
  gene_id = test_data_genes,
  p_value = runif(10),
  direction = sample(c(-1, 1), size = 10, replace = TRUE)
)

test_importer <- CoxpresDbImporter(db_archive = test_data_uncompressed)

###############################################################################

test_that("run_coex_partner_workflow: invalid input", {

  # - fail on NULL gene_ids
  expect_error(
    object = run_coex_partner_workflow(),
    info = "gene_ids should not be missing"
  )
  expect_error(
    object = run_coex_partner_workflow(gene_ids = NULL),
    info = "gene_ids should not be NULL"
  )

  # - fail on NULL gene_statistics
  expect_error(
    object = run_coex_partner_workflow(
      gene_ids = test_data_genes
    ),
    info = "gene_statistics should not be missing"
  )

  expect_error(
    object = run_coex_partner_workflow(
      gene_ids = test_data_genes,
      gene_statistics = NULL
    ),
    info = "gene_statistics should not be NULL"
  )
  # - fail on non-CoxpresDbImporter as importer
  expect_error(
    object = run_coex_partner_workflow(
      gene_ids = test_data_genes,
      gene_statistics = test_gene_statistics
    ),
    info = "importer should not be missing"
  )
  expect_error(
    object = run_coex_partner_workflow(
      gene_ids = test_data_genes,
      gene_statistics = test_gene_statistics,
      importer = "NOT AN IMPORTER"
    ),
    info = "importer should be a CoxpresDbImporter"
  )
  # - fail on gene_statistics that doesn't have gene_id, p_value, direction
  # columns
  expect_error(
    object = run_coex_partner_workflow(
      gene_ids = test_data_genes,
      gene_statistics = tibble::data_frame(
        Gene_id = "a", P_VALUE = 1, direktion = -1
      ),
      importer = test_importer
    ),
    info = "gene_statistics must pass .is_gene_statistics_df"
  )
  # - fail on gene_statistics that is a subset of gene_universe
  expect_error(
    object = run_coex_partner_workflow(
      gene_ids = test_data_genes,
      gene_statistics = test_gene_statistics[-1, ],
      importer = test_importer,
      gene_universe = test_data_genes
    ),
    info = paste(
      "if gene_universe is provided, all it's elements should be present in",
      "gene_statistics$gene_id"
    )
  )

  # - fail on importer for which it's gene-ids are a subset of gene_universe
  expect_error(
    object = run_coex_partner_workflow(
      gene_ids = test_data_genes[1:2],
      gene_statistics = test_gene_statistics,
      importer = test_importer,
      gene_universe = c(test_data_genes, "not_a_gene")
    ),
    info = paste(
      "if gene_universe is provided, all it's elements should have data",
      "inside the coexpression database that is being mined"
    )
  )

  # - fail on gene_universe that is a subset of gene_ids
  expect_error(
    object = run_coex_partner_workflow(
      gene_ids = test_data_genes,
      gene_statistics = test_gene_statistics,
      importer = test_importer,
      gene_universe = test_data_genes[-1]
    ),
    info = "if gene_universe is provided, it should contain all of `gene_id`s"
  )
})

test_that("run_coex_partner_workflow: valid input", {

  # output should be a CoxpresDbPartners object
  expect_is(
    object = run_coex_partner_workflow(
      gene_ids = test_data_genes,
      gene_statistics = test_gene_statistics,
      importer = test_importer
    ),
    class = "CoxpresDbPartners",
    info = paste(
      "the output from `get_gene_universe` should be a `CoxpresDbPartners`",
      "object"
    )
  )

  # when gene_universe is NULL, it should equal the intersection between
  # the gene-ids present in CoxpresDB and in `gene_statistics`
  non_overlap_gene_statistics <- bind_rows(
    test_gene_statistics[-c(1:2), ],
    tibble::data_frame(
      gene_id = "1234567", p_value = 0.4, direction = -1
    )
  )

  expect_equal(
    object = get_gene_universe(
      run_coex_partner_workflow(
        gene_ids = test_data_genes[3:6],
        gene_statistics = non_overlap_gene_statistics,
        importer = test_importer,
        gene_universe = NULL
      )
    ),
    expected = intersect(
      non_overlap_gene_statistics$gene_id,
      get_gene_ids(test_importer)
    ),
    info = paste(
      "If `gene_universe` is NULL, it should be defined to contain only those",
      "genes that are in both `gene_statistics` and accessible through",
      "`importer`"
    )
  )

  # `gene_statistics` field should be the restriction of `gene_statistics`
  # argument to the `gene_universe`
  expect_equal(
    object = run_coex_partner_workflow(
      gene_ids = test_data_genes[1:2],
      gene_statistics = test_gene_statistics,
      importer = test_importer,
      gene_universe = test_data_genes[1:4]
    )@gene_statistics,
    expected = dplyr::filter_(
      test_gene_statistics,
      ~ gene_id %in% test_data_genes[1:4]
    ),
    info = paste(
      "The `gene statistics` field stored by `run_coex_partner_workflow`",
      "should be the restriction of the input `gene_statistics` arg to the",
      "selected `gene_universe`"
    )
  )

  # `partners` entry should be the same as calling get_coex_partners() with the
  # same input arguments

  # `partner_summaries` should be the same as calling evaluate_coex_partners()
  # with the same input

  # `cluster_graph`
})
