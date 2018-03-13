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
  # TODO - pass this test:
  # - fail on gene_statistics that is a subset of gene_universe
  expect_error(
    object = run_coex_partner_workflow(
      gene_ids = test_data_genes,
      gene_statistics = test_gene_statistics,
      importer = test_importer,
      gene_universe = test_data_genes[-1]
    ),
    info = paste(
      "if gene_universe is provided, all genes in gene_statistics should be",
      "present in gene_universe"
    )
  )
  # - fail on gene_universe that is a subset of gene_ids
  #
})

test_that("run_coex_partner_workflow: valid input", {

  # output should be a CoxpresDbPartners object

  # `gene_stats` entry should be the restriction of `gene_statistics` input to
  # `gene_universe`

  # `partners` entry should be the same as calling get_coex_partners() with the
  # same input arguments

  # `partner_summaries` should be the same as calling evaluate_coex_partners()
  # with the same input

  # `cluster_graph`
})
